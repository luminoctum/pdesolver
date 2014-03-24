#ifndef PRIMITIVE
#define PRIMITIVE
#include "utils.h"
#include "Stencil.h"
#include "Boundary.h"
#include "MicroPhysics.h"

template<int DIM, class ET>
class Primitive{
    typedef ET Element_t;
    enum {dim = DIM};
};
//explicit specializations have to be at namespace scope

template<int N> struct Variable{};

template<>
struct Variable<0>{
    typedef double Element_t;
    typedef Array<2, Element_t> Array_t;
    Variable(){};
    Variable(Interval<2> cxy, Interval<2> sXy, Interval<2> sxY){
        cell.initialize(cxy);
        wallXm.initialize(sXy);
        wallXp.initialize(sXy);
        wallYm.initialize(sxY);
        wallYp.initialize(sxY);
        fluxX.initialize(sXy);
        fluxY.initialize(sxY);
    }
    static std::string name;
    static Array_t cell, wallXm, wallXp, wallYm, wallYp, fluxX, fluxY;
    typedef Int2Type<0> ctag;
    typedef Int2Type<1> ptag;
    typedef Boundary<Reflective, ConstExtrap, LinearExtrap, ConstExtrap> bd;
};
std::string Variable<0>::name = "uwind";

template<>
struct Variable<1>{
    typedef double Element_t;
    typedef Array<2, Element_t> Array_t;
    Variable(){};
    Variable(Interval<2> cxy, Interval<2> sXy, Interval<2> sxY){
        cell.initialize(cxy);
        wallXm.initialize(sXy);
        wallXp.initialize(sXy);
        wallYm.initialize(sxY);
        wallYp.initialize(sxY);
        fluxX.initialize(sXy);
        fluxY.initialize(sxY);
    }
    static std::string name;
    static Array_t cell, wallXm, wallXp, wallYm, wallYp, fluxX, fluxY;
    typedef Int2Type<0> ctag;
    typedef Int2Type<1> ptag;
    typedef Boundary<Reflective, ConstExtrap, LinearExtrap, ConstExtrap> bd;
};
std::string Variable<1>::name = "vwind";

template<>
struct Variable<2>{
    typedef double Element_t;
    typedef Array<2, Element_t> Array_t;
    Variable(){};
    Variable(Interval<2> cxy, Interval<2> sXy, Interval<2> sxY){
        cell.initialize(cxy);
        wallXm.initialize(sXy);
        wallXp.initialize(sXy);
        wallYm.initialize(sxY);
        wallYp.initialize(sxY);
        fluxX.initialize(sXy);
        fluxY.initialize(sxY);
    }
    static std::string name;
    static Array_t cell, wallXm, wallXp, wallYm, wallYp, fluxX, fluxY;
    typedef Int2Type<0> ctag;
    typedef Int2Type<1> ptag;
    typedef Boundary<Mirror, FixedConst<Array<1, Element_t>, 2>, ConstExtrap, ConstExtrap> bd;
};
std::string Variable<2>::name = "tc";

template<>
struct Variable<3>{
    typedef Vector<2> Element_t;
    typedef Array<2, Element_t> Array_t;
    Variable(){};
    Variable(Interval<2> cxy, Interval<2> sXy, Interval<2> sxY){
        cell.initialize(cxy);
        wallXm.initialize(sXy);
        wallXp.initialize(sXy);
        wallYm.initialize(sxY);
        wallYp.initialize(sxY);
        fluxX.initialize(sXy);
        fluxY.initialize(sxY);
    }
    static std::string name[2];
    static Array_t cell, wallXm, wallXp, wallYm, wallYp, fluxX, fluxY;
    typedef Int2Type<1> ctag;
    typedef Int2Type<1> ptag;
    typedef Boundary<Mirror, FixedConst<Array<1, Element_t>, 3>, LinearExtrap, ConstExtrap> bd;
};
std::string Variable<3>::name[2] = {"xH2O", "xNH3"};

template<>
struct Variable<4>{
    typedef double Element_t;
    typedef Array<2, Element_t> Array_t;
    Variable(){};
    Variable(Interval<2> cxy, Interval<2> sXy, Interval<2> sxY){
        cell.initialize(cxy);
    }
    static std::string name;
    static Array_t cell; 
    typedef Int2Type<0> ptag;
    typedef Boundary<Reflective, FixedConst<Array<1, Element_t>, 4>, LinearExtrap, ConstExtrap> bd;
};
std::string Variable<4>::name = "phi";

template<>
struct Variable<5>{
    typedef double Element_t;
    typedef Array<2, Element_t> Array_t;
    Variable(){};
    Variable(Interval<2> cxy, Interval<2> sXy, Interval<2> sxY){
        cell.initialize(cxy);
    }
    static std::string name;
    static Array_t cell;
    typedef Int2Type<0> ptag;
    typedef Boundary<Reflective, ConstExtrap, ConstExtrap, FixedConst<Array<1, Element_t>, 5> > bd;
};
std::string Variable<5>::name = "wwind";

template<>
struct Variable<6>{
    typedef Vector<2> Element_t;
    typedef Array<2, Element_t> Array_t;
    Variable(){};
    Variable(Interval<2> cxy, Interval<2> sXy, Interval<2> sxY){
        cell.initialize(cxy);
    }
    static std::string name[2];
    static Array_t cell;
    typedef Int2Type<0> ptag;
    typedef Boundary<Reflective, FixedConst<Array<1, Element_t>, 6>, LinearExtrap, ConstExtrap> bd;
};
std::string Variable<6>::name[2] = {"temp", "tempv"};

template<>
struct Variable<7>{
    typedef Vector<6> Element_t;
    typedef Array<2, Element_t> Array_t;
    Variable(){};
    Variable(Interval<2> cxy, Interval<2> sXy, Interval<2> sxY){
        cell.initialize(cxy);
    }
    static std::string name[6];
    static Array_t cell;
    typedef Int2Type<0> ptag;
    typedef Boundary<Reflective, FixedConst<Array<1, Element_t>, 7>, ConstExtrap, ConstExtrap> bd;
};
std::string Variable<7>::name[6] = {"svpH2O", "svpNH3", "hH2O", "hNH3", "qH2O", "qNH3"};

template<>
struct Variable<8>{
    typedef double Element_t;
    typedef Array<2, Element_t> Array_t;
    Variable(){};
    Variable(Interval<2> cxy, Interval<2> sXy, Interval<2> sxY){
        cell.initialize(cxy);
    }
    static std::string name;
    static Array_t cell; 
    typedef Int2Type<0> ptag;
    typedef Boundary<Reflective, FixedConst<Array<1, Element_t>, 8>, LinearExtrap, ConstExtrap> bd;
};
std::string Variable<8>::name = "pdry";

/* explicit definition for static variables */
// 0 
Variable<0>::Array_t Variable<0>::cell;
Variable<0>::Array_t Variable<0>::wallXm;
Variable<0>::Array_t Variable<0>::wallXp;
Variable<0>::Array_t Variable<0>::wallYm;
Variable<0>::Array_t Variable<0>::wallYp;
Variable<0>::Array_t Variable<0>::fluxX;
Variable<0>::Array_t Variable<0>::fluxY;

// 1
Variable<1>::Array_t Variable<1>::cell;
Variable<1>::Array_t Variable<1>::wallXm;
Variable<1>::Array_t Variable<1>::wallXp;
Variable<1>::Array_t Variable<1>::wallYm;
Variable<1>::Array_t Variable<1>::wallYp;
Variable<1>::Array_t Variable<1>::fluxX;
Variable<1>::Array_t Variable<1>::fluxY;

// 2 
Variable<2>::Array_t Variable<2>::cell;
Variable<2>::Array_t Variable<2>::wallXm;
Variable<2>::Array_t Variable<2>::wallXp;
Variable<2>::Array_t Variable<2>::wallYm;
Variable<2>::Array_t Variable<2>::wallYp;
Variable<2>::Array_t Variable<2>::fluxX;
Variable<2>::Array_t Variable<2>::fluxY;

// 3 
Variable<3>::Array_t Variable<3>::cell;
Variable<3>::Array_t Variable<3>::wallXm;
Variable<3>::Array_t Variable<3>::wallXp;
Variable<3>::Array_t Variable<3>::wallYm;
Variable<3>::Array_t Variable<3>::wallYp;
Variable<3>::Array_t Variable<3>::fluxX;
Variable<3>::Array_t Variable<3>::fluxY;

//4,5,6,7,8
Variable<4>::Array_t Variable<4>::cell;
Variable<5>::Array_t Variable<5>::cell;
Variable<6>::Array_t Variable<6>::cell;
Variable<7>::Array_t Variable<7>::cell;
Variable<8>::Array_t Variable<8>::cell;

/* Axisymmetric Geostrophic Adjustment */
template<class ET>
class Primitive<2, ET>{
    enum {dim = 2};
    typedef ET Element_t;
    typedef Interval<dim> Interval_t;
protected:
    Stencil<CenterDifferenceX> cdx;
    ForwardIntegralY finty;
    BackwardIntegralY binty;
    double eps1, eps2, grav, T0, cp, f;
public:
    Variable<0> uwind;
    Variable<1> vwind;
    Variable<2> theta;
    Variable<3> mixr;
    Variable<4> phi;
    Variable<5> wwind;
    Variable<6> temp;
    Variable<7> moist; //0:svpH2O, 1:svpNH3, 2:hH2O, 3:hNH3, 4:qH2O, 5:qNH3
    Variable<8> pdry; 
    Water H2O; Ammonia NH3;
    Array<2, double> svpH2O, svpNH3;
    Primitive(){
        // put scalar here, but should be determined from attributes in ncfile
        eps1 = 8.21;
        eps2 = 7.75;
        grav = 10.44;
        T0 = 134.8;
        cp = 11455.;
        f = 2.128e-4;
        H2O.mmr = 1.3E-2;
        NH3.mmr = 4.E-4;
    };
    template<class A>
    void updateDiagnostics(A& ncvar, Interval<2> cij) {
        wwind.cell(cij) = - binty(cdx(uwind.cell, cij), cij);
        temp.cell(cij).comp(0) = ncvar["t_ov_tc"] * theta.cell(cij);
        temp.cell.comp(1) = temp.cell.comp(0) *
            (1. + mixr.cell.comp(0) + mixr.cell.comp(1)) 
            / (1. + eps1 * mixr.cell.comp(0) + eps2 * mixr.cell.comp(1)); 
        // it should use the integ stencil but the compiler does not let it go
        for (int j = cij[1].first() + 1; j <= cij[1].last(); j++)
            phi.cell(cij[0], j) = phi.cell(cij[0], j - 1) 
                + 0.5 * finty.dy * grav / T0 
                * (temp.cell(cij[0], j - 1).comp(1) + temp.cell(cij[0], j).comp(1) 
                - ncvar["tv0"](cij[0], j - 1) - ncvar["tv0"](cij[0], j));
        // calculate saturation vapor pressure
        // should define a virtual function that takes care of svp_from_t 
        // test condensing
        //mixr.cell.comp(0) = max(mixr.cell.comp(0));
        //mixr.cell.comp(1) = max(mixr.cell.comp(1));

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
    }
};

#endif
