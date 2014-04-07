#ifndef PRIMITIVE
#define PRIMITIVE
#include <vector>
#include "Variable.h"
#include "Stencil.h"
#include "MicroPhysics.h"

class ModelBase{
protected:
    Stencil<ElementSum> esum;
    Stencil<ElementAverage> eavg;
    Stencil<Average> avg;
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

class Primitive : public ModelBase{
    typedef _AtomicType_ T;
protected:
    Stencil<_FluxCalculator_<T, DimX> > fux;
    Stencil<_FluxCalculator_<T, DimY> > fuy;

    double grav, T0, cp, f; Vector<2, T> eps;
    Array<2, T> mass, rdist, t_ov_tc, tv0, ptol, massx, massy;

    VariableList vlist;
    Water H2O; Ammonia NH3;
    CondensateList<2> species;
public:
    Variable<T, Boundary<Reflective, ConstExtrap, ConstExtrap, ConstExtrap>
    //Variable<T, Boundary<FixedConst<T>, ConstExtrap, ConstExtrap, ConstExtrap>
        , MassView> uwind;
    Variable<T, Boundary<Reflective, ConstExtrap, ConstExtrap, ConstExtrap> 
        > vwind;
    Variable<T, Boundary<Mirror, FixedConst<T>, ConstExtrap, FixedConst<T> >
        > theta;
    Variable<Vector<2, T>, Boundary<Mirror, FixedConst<Vector<2, T> >, LinearExtrap, ConstExtrap>
        > mixr;
    Variable<Vector<2, T>, Boundary<Mirror, FixedConst<Vector<2, T> >, ConstExtrap, ConstExtrap>
        , AverageView> svp, rh, liq;
    Variable<T, Boundary<Mirror, Dependent, ConstExtrap, FixedConst<T> >
        , MassView> wwind;
    Variable<T, Boundary<LinearExtrap, FixedConst<T>, LinearExtrap, ConstExtrap>
        > phi, temp, tempv, pdry;

    Primitive() : uwind("uwind"), vwind("vwind"), theta("tc"), mixr({"xH2O", "xNH3"}), 
        wwind("wwind"), phi("phi"), temp("temp"), tempv("tempv"), pdry("pdry"),
        svp({"svpH2O", "svpNH3"}), rh({"hH2O", "hNH3"}), liq({"qH2O", "qNH3"}),
        H2O(1.3E-2), NH3(4E-4) ,species(H2O, NH3)
    {
        // put scalar here, but should be determined from attributes in ncfile
        eps     = species.mu / 2.198E-3;
        grav    = 10.44;
        T0      = 134.8;
        cp      = 11455.;
        f       = 2.128e-4;

        fux.function().xwind = &uwind.wallx;
        fuy.function().ywind = &wwind.wally;

        _ncInitialize7_(mass, rdist, t_ov_tc, tv0, ptol, massx, massy);

        vlist.mass = &mass;
        vlist.var = {
            &uwind, &vwind, &wwind, &theta, &mixr, &phi, &temp, &tempv, &pdry, &svp, &rh, &liq
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
        uwind.tendency += - fdx(avg(phi.cell(sicj)) * massx, cij) + mass * phi.cell(cij) / rdist
            + mass * vwind.cell(cij) * (f + vwind.cell(cij) / rdist);
        vwind.tendency += - uwind.cell(cij) / mass * (f + vwind.cell(cij) / rdist);
        theta.tendency += theta.cell(cij) / (cp * temp.cell(cij)) * esum(species.Lv * eps * liq.cell(cij)) / setups::dt;
    }
    void updateMicrophysics(){
        // saturation vapor pressure
        svp.cell.comp(0) = H2O.svp_from_t(temp.cell);
        svp.cell.comp(1) = NH3.svp_from_t(temp.cell);
        // relative humidity: h = x * pdry / svp
        rh.cell = mixr.cell * pdry.cell / svp.cell;
        // condensing liquid
        for (int s = 0; s < 2; s++){
            //rh.cell.comp(s) = where(rh.cell.comp(s) < 0., 0., rh.cell.comp(s));
            liq.cell.comp(s) = where(rh.cell.comp(s) > 1.2, (rh.cell.comp(s) - 1.) * svp.cell.comp(s) / pdry.cell, 0.);
            rh.cell.comp(s) = where(liq.cell.comp(s) > 0., 1., rh.cell.comp(s));
        }
        pdry.cell(cij) = ptol / (1. + esum(mixr.cell, cij));
        // mixing ratio: x = h * svp / pdry
        mixr.cell(cij) = rh.cell(cij) * svp.cell(cij) / pdry.cell(cij);
    }
    void advection(){
        vlist.wallConstruct();
        /* uwind at the left boundary must be zero */
        uwind.wallx(0, AllDomain<1>()).comp(0) = 0;
        uwind.wallx(-1, AllDomain<1>()).comp(1) = 0;
    
        binty(- fdx(fux(onesx)), wwind.wally, cij);
        uwind.tendency -= fdx(fux(uwind.wallx) / massx) + fdy(fuy(uwind.wally) / massy);
        vwind.tendency -= (fdx(fux(vwind.wallx)) + fdy(fuy(vwind.wally))) / mass;
        theta.tendency -= (fdx(fux(theta.wallx)) + fdy(fuy(theta.wally))) / mass;
        mixr.tendency  -= (fdx(fux(mixr.wallx)) + fdy(fuy(mixr.wally))) / mass;
    }
    void dissipation(){
        uwind.tendency += 0.03 / setups::dt * laplace(uwind.cell, cij);
        vwind.tendency += 0.03 / setups::dt * laplace(vwind.cell, cij);
    }
    void observe(int i, double _time){
        printf("%-6d%-8.1f%-8.1f%-6.0f%-6.0f%-12.2E%-12.2E%-12.2E\n", i, _time, 
                max(uwind.cell(cij) / mass), 
                min(theta.cell(cij)), max(theta.cell(cij)), 
                max(mixr.cell(cij).comp(0)), max(mixr.cell(cij).comp(1)),
                max(cdy(wwind.cell, cij) + cdx(uwind.cell, cij))
                );
        vlist.ncwrite(_time);
    }
    void forward(){
        vlist.step_count++;
        advection();
        updateBodyforce();
        dissipation();
        updateMicrophysics();
        vlist.updateTendency();
        updateDiagnostics();
        vlist.fixBoundary();
        vlist.updateView();
    }
};

#endif
