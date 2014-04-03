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
    Stencil<ForwardDifferenceX> fdx;
    Stencil<ForwardDifferenceY> fdy;
    Stencil<GodunovX<T> > flux_x;
    Stencil<GodunovY<T> > flux_y;
    ForwardIntegralY finty;
    BackwardIntegralY binty;
    double grav, T0, cp, f; Vector<2> eps;
    VariableList<T> vlist;
    Array<2, T> rdist, t_ov_tc, tv0, ptol;
    Interval<2> cij, sicj, cisj;
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
        eps     = species.mu / 2.198E-3;
        grav    = 10.44;
        T0      = 134.8;
        cp      = 11455.;
        f       = 2.128e-4;
        cij     = setups::cij;
        sicj    = setups::sicj;
        cisj    = setups::cisj;

        CenterDifferenceX::dx = setups::dx;
        ForwardDifferenceX::dx = setups::dx;
        ForwardDifferenceY::dy = setups::dy;
        ForwardIntegralY::dy = setups::dy;
        BackwardIntegralY::dy = setups::dy;
        flux_x.function().xwind = &uwind.wallx;
        flux_y.function().ywind = &wwind.wally;

        rdist.initialize(cij); rdist = setups::ncvar["rdist"];
        t_ov_tc.initialize(cij); t_ov_tc = setups::ncvar["t_ov_tc"];
        tv0.initialize(cij); tv0 = setups::ncvar["tv0"];
        ptol.initialize(cij); ptol = setups::ncvar["ptol"];

        vlist.mass = &mass;
        vlist.mvar = {
            &uwind, &wwind
        };
        vlist.uvar = {
            &vwind, &theta, &mixr, &phi, &temp, &tempv, &pdry, &svp, &rh, &liq
        };
        vlist.pvar = {
            &uwind, &vwind, &theta, &mixr
        };
        vlist.fixBoundary();
        //vlist.printInfo(1);
    };
    void updateDiagnostics(){
        wwind.cell(cij) = - binty(cdx(uwind.cell, cij), cij);
        temp.cell(cij) = t_ov_tc * theta.cell(cij);
        // stencil output is zero based, need to specify cij explicitly
        tempv.cell(cij) = temp.cell(cij) * (1. + esum(mixr.cell, cij))/(1. + esum(eps * mixr.cell, cij));
        phi.cell(cij) = finty(grav / T0 * (tempv.cell(cij) - tv0), cij);
    }
    void updateBodyforce(){
        uwind.cell_t += - mass.cell(cij) * (cdx(phi.cell, cij) + vwind.cell(cij) * (f + vwind.cell(cij) / rdist));
        vwind.cell_t += - uwind.cell(cij) / mass.cell(cij) * (f + vwind.cell(cij) / rdist);
        theta.cell_t += theta.cell(cij) / (cp * temp.cell(cij)) * esum(species.Lv * eps * liq.cell(cij));
    }
    void updateMicrophysics(){
        // test condensing
        //mixr.cell.comp(0) = max(mixr.cell.comp(0));
        //mixr.cell.comp(1) = max(mixr.cell.comp(1));

        // saturation vapor pressure
        svp.cell.comp(0) = H2O.svp_from_t(temp.cell);
        svp.cell.comp(1) = NH3.svp_from_t(temp.cell);
        // relative humidity: h = x * pdry / svp
        rh.cell = mixr.cell * pdry.cell / svp.cell;
        // condensing liquid
        for (int s = 0; s < 2; s++){
            liq.cell.comp(s) = where(rh.cell.comp(s) > 1.2, (rh.cell.comp(s) - 1.) * svp.cell.comp(s) / pdry.cell, 0.);
            rh.cell.comp(s) = where(liq.cell.comp(s) > 0., 1., rh.cell.comp(s));
        }
        pdry.cell(cij) = ptol / (1. + esum(mixr.cell, cij));
        // mixing ratio: x = h * svp / pdry
        mixr.cell(cij) = rh.cell(cij) * svp.cell(cij) / pdry.cell(cij);
        liq.cell /= setups::dt;
    }
    void advection(){
        wwind.wallConstruct();
        vlist.wallConstruct();
        uwind.cell_t += fdx(flux_x(uwind.wallx(sicj))) + fdy(flux_y(uwind.wally(cisj)));
        vwind.cell_t += (fdx(flux_x(vwind.wallx(sicj))) + fdy(flux_y(vwind.wally(cisj)))) / mass.cell(cij);
        theta.cell_t += (fdx(flux_x(theta.wallx(sicj))) + fdy(flux_y(theta.wally(cisj)))) / mass.cell(cij);
        mixr.cell_t  += (fdx(flux_x(mixr.wallx(sicj))) + fdy(flux_y(mixr.wally(cisj)))) / mass.cell(cij);
    }
    void observe(double _time){
        vlist.ncwrite(_time);
    }
    void forward(){
        advection();

        updateBodyforce();

        updateMicrophysics();

        vlist.updateTendency();

        updateDiagnostics();

        vlist.fixBoundary();

        //vlist.printInfo(1);
        observe(1.);
    }
};

#endif
