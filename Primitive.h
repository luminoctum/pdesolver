#ifndef PRIMITIVE
#define PRIMITIVE
#include <iostream>
#include "Variable.h"
#include "utils.h"
#include <vector>
#include "Stencil.h"
#include "MicroPhysics.h"
using namespace std;

template<int O, class T = double>
class Primitive{
protected:
    Stencil<CenterDifferenceX> cdx;
    ForwardIntegralY finty;
    BackwardIntegralY binty;
    double eps1, eps2, grav, T0, cp, f;
    vector<VariableBase*> var;
    Array<2, T> mass, rdist, t_ov_tc, tv0;
    Interval<2> cij;
public:
    Variable<O, T,
        Boundary<Reflective, ConstExtrap, LinearExtrap, ConstExtrap>
            > uwind, vwind;
    Variable<O, T,
        Boundary<Mirror, FixedConst<T>, ConstExtrap, ConstExtrap>
            > theta;
    Variable<O, Vector<2, T>,
        Boundary<Mirror, FixedConst<Vector<2, T> >, LinearExtrap, ConstExtrap>
            > mixr;
    Variable<O, Vector<2, T>,
        Boundary<Mirror, FixedConst<Vector<2, T> >, ConstExtrap, ConstExtrap>
            > svp, rh, liq;
    Variable<O, T,
        Boundary<Mirror, ConstExtrap, ConstExtrap, FixedConst<T> >
            > wwind;
    Variable<O, T,
        Boundary<Mirror, FixedConst<T>, LinearExtrap, ConstExtrap>
            > phi;
    Variable<O, T,
        Boundary<Mirror, FixedConst<T>, LinearExtrap, ConstExtrap>
            > temp, tempv;
    Variable<O, T,
        Boundary<Mirror, FixedConst<T>, LinearExtrap, ConstExtrap>
            > pdry;

    Primitive() :
        uwind("uwind"), vwind("vwind"), theta("tc"), mixr((char*[]){"xH2O", "xNH3"}), 
        wwind("wwind"), phi("phi"), temp("temp"), tempv("tempv"), pdry("pdry"),
        svp((char*[]){"svpH2O", "svpNH3"}),
        rh((char*[]){"hH2O", "hNH3"}), liq((char*[]){"qH2O", "qNH3"})
    {
        // put scalar here, but should be determined from attributes in ncfile
        eps1 = 8.21; eps2 = 7.75;
        grav = 10.44;
        T0 = 134.8;
        cp = 11455.;
        f = 2.128e-4;
        cij = setups::cij;
        mass.initialize(cij); mass = setups::ncvar["mass"];
        rdist.initialize(cij); rdist = setups::ncvar["rdist"];
        t_ov_tc.initialize(cij); t_ov_tc = setups::ncvar["t_ov_tc"];
        tv0.initialize(cij); tv0 = setups::ncvar["tv0"];

        var = {
            &uwind, &vwind, &theta, &mixr, 
            &wwind, &phi, &temp, &tempv, &pdry, 
            &svp, &rh, &liq
        };
        for (int i = 0; i < var.size(); i++){
            var[i]->fixBoundary();
            var[i]->printInfo(1);
        }
    };
    void updateDiagnostics() {
        wwind.cell(cij) = - binty(cdx(uwind.cell, cij), cij);
        temp.cell(cij) = t_ov_tc * theta.cell(cij);
        tempv.cell = temp.cell *
            (1. + mixr.cell.comp(0) + mixr.cell.comp(1)) 
            / (1. + eps1 * mixr.cell.comp(0) + eps2 * mixr.cell.comp(1)); 
        // it should use the integ stencil but the compiler does not let it go
        for (int j = cij[1].first() + 1; j <= cij[1].last(); j++)
            phi.cell(cij[0], j) = phi.cell(cij[0], j - 1) 
                + 0.5 * finty.dy * grav / T0 
                * (tempv.cell(cij[0], j - 1) + tempv.cell(cij[0], j)
                - tv0(cij[0], j - 1) - tv0(cij[0], j));
    }
    /*
    void updateBodyforce(){
        uwind.cell_t(cij) += - mass * (cdx(phi.cell, cij) + vwind.cell(cij) * (f + vwind.cell(cij) / rdist));
        vwind.cell_t(cij) += - uwind.cell(cij) / mass * (f + vwind.cell(cij) / rdist);
        theta.cell_t(cij) += theta.cell(cij) / (cp * temp.cell(cij)) * eSum(speci.Lv * speci.eps * liq.cell(cij));
    }
    void updateMicrophysics(){
        // calculate saturation vapor pressure
        // should define a virtual function that takes care of svp_from_t 
        // test condensing
        //mixr.cell.comp(0) = max(mixr.cell.comp(0));
        //mixr.cell.comp(1) = max(mixr.cell.comp(1));

        // disable microphysics
        moist.cell.comp(0) = H2O.svp_from_t(temp.cell.comp(0));
        moist.cell.comp(1) = NH3.svp_from_t(temp.cell.comp(0));
        // should define a binary function min, max, etc...
        for (int s = 0; s < 2; s++){
            // calculate relative humidity: h = x * pdry / svp
            moist.cell(cij).comp(s + 2) = ncvar["pdry"] * mixr.cell(cij).comp(s) / moist.cell(cij).comp(s);
            for (int i = cij[0].first(); i <= cij[0].last(); i++)
                for (int j = cij[1].first(); j <= cij[1].last(); j++)
                    if (moist.cell(i, j)(s + 2) > 1.2){
                        // calculate condenstates: q = (h - 1.) * svp / pdry
                        moist.cell(i, j)(s + 4) = (moist.cell(i, j)(s + 2) - 1.) * moist.cell(i, j)(s) / ncvar["pdry"](i, j);
                        moist.cell(i, j)(s + 2) = 1.;
                    }
        }
        pdry.cell(cij) = ncvar["ptol"] / (1. + mixr.cell(cij).comp(0) + mixr.cell(cij).comp(1));
        // calculate mixing ratio: x = h * svp / pdry
        for (int s = 0; s < 2; s++)
            mixr.cell(cij).comp(s) = moist.cell(cij).comp(s + 2) * moist.cell(cij).comp(s) / pdry.cell(cij);
    }*/
};

#endif
