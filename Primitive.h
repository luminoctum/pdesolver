#ifndef PRIMITIVE
#define PRIMITIVE
#include <iostream>
#include "Variable.h"
#include "utils.h"
#include <vector>
#include "Stencil.h"
#include "MicroPhysics.h"

class ModelBase{
protected:
    Stencil<ElementSum> esum;
    Stencil<ElementAverage> eavg;
    Stencil<FiniteDifference<D1, O2, DimX> > cdx;
    Stencil<FiniteDifference<D1, O2, DimY> > cdy;
    Stencil<FiniteDifference<D1, O1, DimX> > fdx;
    Stencil<FiniteDifference<D1, O1, DimY> > fdy;
    Integrate<Forward, Trapz, DimY> finty;
    Integrate<Backward, Staggered, DimY> binty;
    Stencil<Laplace> laplace;
    Array<2, double, ConstantFunction> ones;
    Array<2, Vector<2, double>, ConstantFunction> onesx, onesy;
    Interval<2> cij, sicj, cisj;
    struct RunTimeInfo{
        std::vector<std::string> header_fmt;
        std::vector<std::string> header;
        std::vector<std::string> value_fmt;
        std::vector<double> value;
        void printHeader(){
            for (int i = 0; i < header.size(); i++)
                printf(header_fmt[i].c_str(), header[i].c_str());
            std::cout << std::endl;
        }
    } runInfo;
    ModelBase(){
        cij     = setups::cij;
        sicj    = setups::sicj;
        cisj    = setups::cisj;
        FiniteDifferenceBase::dx = setups::dx;
        FiniteDifferenceBase::dy = setups::dy;
        ones.initialize(cij); ones.engine().setConstant(1.);
        onesx.initialize(sicj); onesx.engine().setConstant(Vector<2>(1., 1.));
        onesy.initialize(cisj); onesy.engine().setConstant(Vector<2>(1., 1.));
    }
};

template<int O, class T = double>
class Primitive : public ModelBase{
protected:
    Stencil<Godunov<T, DimX> > fux;
    Stencil<Godunov<T, DimY> > fuy;

    double grav, T0, cp, f; Vector<2> eps;
    VariableList<T> vlist;
    Array<2, T> rdist, t_ov_tc, tv0, ptol, imassx, imassy;
    Water H2O; Ammonia NH3;
    CondensateList<2> species;
public:
    Variable<O, T,
        Boundary<Reflective, ConstExtrap, ConstExtrap, ConstExtrap>
        //Boundary<FixedConst<T>, ConstExtrap, ConstExtrap, ConstExtrap>
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
        Boundary<Mirror, Dependent, ConstExtrap, FixedConst<T> >
            > wwind;
    Variable<O, T,
        Boundary<Mirror, FixedConst<T>, LinearExtrap, ConstExtrap>
        //Boundary<ConstExtrap, FixedConst<T>, LinearExtrap, ConstExtrap>
            > phi;
    Variable<O, T,
        Boundary<Mirror, FixedConst<T>, LinearExtrap, ConstExtrap>
        //Boundary<ConstExtrap, FixedConst<T>, LinearExtrap, ConstExtrap>
            > temp, tempv;
    Variable<O, T,
        Boundary<Mirror, FixedConst<T>, LinearExtrap, ConstExtrap>
            > pdry;

    Primitive() : mass("mass"),
        uwind("uwind"), vwind("vwind"), theta("tc"), mixr({"xH2O", "xNH3"}), 
        wwind("wwind"), phi("phi"), temp("temp"), tempv("tempv"), pdry("pdry"),
        svp({"svpH2O", "svpNH3"}), rh({"hH2O", "hNH3"}), liq({"qH2O", "qNH3"}),
        H2O(1.3E-2), NH3(4E-4), species(H2O, NH3)
    {
        // put scalar here, but should be determined from attributes in ncfile
        eps     = species.mu / 2.198E-3;
        grav    = 10.44;
        T0      = 134.8;
        cp      = 11455.;
        f       = 2.128e-4;

        fux.function().xwind = &uwind.wallx;
        fuy.function().ywind = &wwind.wally;

        _ncInitialize6_(rdist, t_ov_tc, tv0, ptol, imassx, imassy);

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
        runInfo.header_fmt = {"%-6s", "%-8s", "%-8s", "%-12s", "%-12s", "%-12s", "%-12s"};
        runInfo.header     = {"No.", "time", "uwind", "--theta--", "xH2O", "xNH3", "mass"};
        runInfo.printHeader();
    };
    void updateDiagnostics(){
        wwind.cell(cij) = eavg(wwind.wally, cij);
        temp.cell(cij) = t_ov_tc * theta.cell(cij);
        // stencil output is zero based, need to specify cij explicitly
        tempv.cell(cij) = temp.cell(cij) * (1. + esum(mixr.cell, cij))/(1. + esum(eps * mixr.cell, cij));
        finty(grav / T0 * (tempv.cell(cij) - tv0), phi.cell, cij);
    }
    void updateBodyforce(){
        uwind.cell_t += - mass.cell(cij) * (cdx(phi.cell, cij) + vwind.cell(cij) * (f + vwind.cell(cij) / rdist));
        vwind.cell_t += - uwind.cell(cij) / mass.cell(cij) * (f + vwind.cell(cij) / rdist);
        //theta.cell_t += theta.cell(cij) / (cp * temp.cell(cij)) * esum(species.Lv * eps * liq.cell(cij)) / setups::dt;
    }
    void updateMicrophysics(){
        // saturation vapor pressure
        svp.cell.comp(0) = H2O.svp_from_t(temp.cell);
        svp.cell.comp(1) = NH3.svp_from_t(temp.cell);
        // relative humidity: h = x * pdry / svp
        rh.cell = mixr.cell * pdry.cell / svp.cell;
        // condensing liquid
        for (int s = 0; s < 2; s++){
            rh.cell.comp(s) = where(rh.cell.comp(s) < 0., 0., rh.cell.comp(s));
            liq.cell.comp(s) = where(rh.cell.comp(s) > 1.2, (rh.cell.comp(s) - 1.) * svp.cell.comp(s) / pdry.cell, 0.);
            rh.cell.comp(s) = where(liq.cell.comp(s) > 0., 1., rh.cell.comp(s));
        }
        pdry.cell(cij) = ptol / (1. + esum(mixr.cell, cij));
        // mixing ratio: x = h * svp / pdry
        mixr.cell(cij) = rh.cell(cij) * svp.cell(cij) / pdry.cell(cij);
    }
    void advection(){
        vlist.wallConstruct();
        binty( - fdx(fux(onesx)), wwind.wally, cij);
        uwind.cell_t -= fdx(fux(uwind.wallx) * imassx) + fdy(fuy(uwind.wally) * imassy);
        vwind.cell_t -= (fdx(fux(vwind.wallx)) + fdy(fuy(vwind.wally))) / mass.cell(cij);
        theta.cell_t -= (fdx(fux(theta.wallx)) + fdy(fuy(theta.wally))) / mass.cell(cij);
        mixr.cell_t  -= (fdx(fux(mixr.wallx)) + fdy(fuy(mixr.wally))) / mass.cell(cij);
    }
    void dissipation(){
        uwind.cell_t += 0.03 / setups::dt * laplace(uwind.cell, cij);
        vwind.cell_t += 0.03 / setups::dt * laplace(vwind.cell, cij);
    }
    void observe(int i, double _time){
        printf("%-6d%-8.1f%-8.1f%-6.0f%-6.0f%-12.2E%-12.2E%-12.2E\n", i, _time, 
                max(uwind.cell(cij) / mass.cell(cij)), 
                min(theta.cell(cij)), max(theta.cell(cij)), 
                max(mixr.cell(cij).comp(0)), max(mixr.cell(cij).comp(1)),
                max(cdy(wwind.cell, cij) + cdx(uwind.cell, cij))
                );
        vlist.ncwrite(_time);
    }
    void forward(){
        advection();
        dissipation();
        updateBodyforce();
        updateMicrophysics();
        vlist.updateTendency();
        vlist.fixBoundary();
        updateDiagnostics();
        //vlist.printInfo(1);
    }
};

#endif
