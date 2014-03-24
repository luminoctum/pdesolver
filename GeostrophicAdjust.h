#ifndef TEST
#define TEST
#include "Boundary.h"

template<int DIM, class ET>
class GeostrophicAdjust{
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
    typedef Vector<4> Element_t;
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
    static std::string name[4];
    static Array_t cell, wallXm, wallXp, wallYm, wallYp, fluxX, fluxY;
    typedef Int2Type<1> ctag;
    typedef Int2Type<1> ptag;
    typedef Boundary<Mirror, FixedConst<Array<1, Element_t>, 3>, LinearExtrap, ConstExtrap> bd;
};
std::string Variable<3>::name[4] = {"xH2O", "xNH3", "qH2O", "qNH3"};

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
    typedef double Element_t;
    typedef Array<2, Element_t> Array_t;
    Variable(){};
    Variable(Interval<2> cxy, Interval<2> sXy, Interval<2> sxY){
        cell.initialize(cxy);
    }
    static std::string name;
    static Array_t cell;
    typedef Int2Type<0> ptag;
    typedef Boundary<Reflective, FixedConst<Array<1, Element_t>, 6>, LinearExtrap, ConstExtrap> bd;
};
std::string Variable<6>::name = "temp";

template<>
struct Variable<7>{
    typedef Vector<2> Element_t;
    typedef Array<2, Element_t> Array_t;
    Variable(){};
    Variable(Interval<2> cxy, Interval<2> sXy, Interval<2> sxY){
        cell.initialize(cxy);
    }
    static std::string name[2];
    static Array_t cell;
    typedef Int2Type<0> ptag;
    typedef Boundary<Reflective, FixedConst<Array<1, Element_t>, 7>, ConstExtrap, ConstExtrap> bd;
};
std::string Variable<7>::name[2] = {"hH2O", "hNH3"};

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

//4,5,6,7
Variable<4>::Array_t Variable<4>::cell;
Variable<5>::Array_t Variable<5>::cell;
Variable<6>::Array_t Variable<6>::cell;
Variable<7>::Array_t Variable<7>::cell;

/* Axisymmetric Geostrophic Adjustment */
template<class ET>
class GeostrophicAdjust<2, ET>{
    enum {dim = 2};
    typedef ET Element_t;
    typedef Interval<dim> Interval_t;
public:
    Variable<0> uwind;
    Variable<1> vwind;
    Variable<2> theta;
    Variable<3> mixr;
    Variable<4> phi;
    Variable<5> wwind;
    Variable<6> temp;
    Variable<7> rh;
};

#endif
