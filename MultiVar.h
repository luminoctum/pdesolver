#ifndef TEST
#define TEST
#include "Boundary.h"

template<int DIM, class ET>
class MultiVar{
    typedef ET Element_t;
    enum {dim = DIM};
};
//explicit specializations have to be at namespace scope

template<int N> struct Variable{};
template<>
struct Variable<0>{
    typedef double Element_t;
    typedef Array<1, Element_t> Array_t;
    Variable(){};
    Variable(Interval<1> cx, Interval<1> sx){
        cell.initialize(cx);
        wallm.initialize(sx);
        wallp.initialize(sx);
    }
    static std::string name;
    static Array_t cell;
    static Array_t wallm;
    static Array_t wallp;
    typedef Int2Type<0> ctag;
    //typedef Boundary<Reflective, LinearExtrap, ConstExtrap, ConstExtrap> bd;
    typedef LinearExtrap bd;
};
std::string Variable<0>::name = "uwind";
template<>
struct Variable<1>{
    typedef Vector<2> Element_t;
    typedef Array<1, Element_t> Array_t;
    Variable(){};
    Variable(Interval<1> cx, Interval<1> sx){
        cell.initialize(cx);
        wallm.initialize(sx);
        wallp.initialize(sx);
    }
    static std::string name;
    static Array_t cell;
    static Array_t wallm;
    static Array_t wallp;
    typedef Int2Type<1> ctag;
    //typedef Boundary<FixedConst<Element_t, 1>, Mirror, ConstExtrap, ConstExtrap> bd;
    typedef ConstExtrap bd;
};
std::string Variable<1>::name = "vwind";

/* explicit definition for static variables */
Variable<0>::Array_t Variable<0>::cell;
Variable<0>::Array_t Variable<0>::wallm;
Variable<0>::Array_t Variable<0>::wallp;
Variable<1>::Array_t Variable<1>::cell;
Variable<1>::Array_t Variable<1>::wallm;
Variable<1>::Array_t Variable<1>::wallp;

/* Test for a system of multi variables */
template<class ET>
class MultiVar<1, ET>{
    enum {dim = 1};
    typedef ET Element_t;
    typedef Interval<dim> Interval_t;
public:
    Variable<0> uwind;
    Variable<1> vwind;
};

#endif
