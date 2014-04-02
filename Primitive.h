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
    Stencil<ElementSum> esum;
    Stencil<CenterDifferenceX> cdx;
    ForwardIntegralY finty;
    BackwardIntegralY binty;
    double grav, T0, cp, f; Vector<2> eps;
    //std::vector<VariableBase*> var;
    //VariableList vlist;
    Array<2, T> rdist, t_ov_tc, tv0, ptol;
    Interval<2> cij;
    Water H2O; Ammonia NH3;
    CondensateList<2> species;
public:
    Variable<O, T,
        Boundary<Reflective, ConstExtrap, LinearExtrap, ConstExtrap>
            > mass, uwind, vwind;
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

    Primitive() : mass("mass"),
        uwind("uwind"), vwind("vwind"), theta("tc"), mixr((char*[]){"xH2O", "xNH3"}), 
        wwind("wwind"), phi("phi"), temp("temp"), tempv("tempv"), pdry("pdry"),
        svp((char*[]){"svpH2O", "svpNH3"}), rh((char*[]){"hH2O", "hNH3"}), liq((char*[]){"qH2O", "qNH3"}),
        H2O(1.3E-2), NH3(4E-4), species(H2O, NH3)
    {
        // put scalar here, but should be determined from attributes in ncfile
        eps = species.mu / 2.2E-3;
        grav = 10.44;
        T0 = 134.8;
        cp = 11455.;
        f = 2.128e-4;
        cij = setups::cij;
        CenterDifferenceX::dx = setups::dx;
        ForwardIntegralY::dy = setups::dy;
        BackwardIntegralY::dy = setups::dy;

        rdist.initialize(cij); rdist = setups::ncvar["rdist"];
        t_ov_tc.initialize(cij); t_ov_tc = setups::ncvar["t_ov_tc"];
        tv0.initialize(cij); tv0 = setups::ncvar["tv0"];
        ptol.initialize(cij); ptol = setups::ncvar["ptol"];

        /*vlist.mass = &mass;
        vlist.mvar = {
            &uwind, &wwind
        };
        vlist.uvar = {
            &vwind, &theta, &mixr, &phi, &temp, &tempv, &pdry, &svp, &rh, &liq
        };
        vlist.pvar = {
            &uwind, &vwind, &theta, &mixr
        };
        vlist.printInfo()*/
    };
    void updateDiagnostics(){
        wwind.cell(cij) = - binty(cdx(uwind.cell, cij), cij);
        temp.cell(cij) = t_ov_tc * theta.cell(cij);
        tempv.cell = temp.cell * (1. + esum(mixr.cell))/(1. + esum(eps * mixr.cell));
        phi.cell(cij) = finty(grav / T0 * (tempv.cell(cij) - tv0), cij);
    }
    void updateBodyforce(){
        uwind.cell_t(cij) += - mass.cell(cij) * (cdx(phi.cell, cij) + vwind.cell(cij) * (f + vwind.cell(cij) / rdist));
        vwind.cell_t(cij) += - uwind.cell(cij) / mass.cell(cij) * (f + vwind.cell(cij) / rdist);
        theta.cell_t(cij) += theta.cell(cij) / (cp * temp.cell(cij)) * esum(species.Lv * eps * liq.cell(cij));
    }
    void updateMicrophysics(){
        // calculate saturation vapor pressure
        // should define a virtual function that takes care of svp_from_t 
        // test condensing
        //mixr.cell.comp(0) = max(mixr.cell.comp(0));
        //mixr.cell.comp(1) = max(mixr.cell.comp(1));

        svp.cell = species.svp_from_t(temp.cell);
        // relative humidity: h = x * pdry / svp
        rh.cell = mixr.cell * pdry.cell / svp.cell;
        for (int s = 0; s < 2; s++){
            for (int i = cij[0].first(); i <= cij[0].last(); i++)
                for (int j = cij[1].first(); j <= cij[1].last(); j++)
                    if (rh.cell(i, j)(s) > 1.2){
                        // calculate condenstates: q = (h - 1.) * svp / pdry
                        liq.cell(i, j)(s) = (rh.cell(i, j)(s) - 1.) * svp.cell(i, j)(s) / pdry.cell(i, j);
                        rh.cell(i, j)(s) = 1.;
                    }
        }
        pdry.cell(cij) = ptol / (1. + esum(mixr.cell, cij));
        // mixing ratio: x = h * svp / pdry
        mixr.cell(cij) = rh.cell(cij) * svp.cell(cij) / pdry.cell(cij);
    }
};

#endif
