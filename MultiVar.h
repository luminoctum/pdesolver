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
    Variable(){};
    Variable(Interval<1> cx, Interval<1> si){
        cell.initialize(cx);
        wall.initialize(si);
    }
    static std::string name;
    static Array<1, Element_t> cell;
    static Array<1, Element_t> wall;
    typedef Int2Type<0> ctag;
    typedef Int2Type<1> ptag;
    //typedef Boundary<Reflective, LinearExtrap, ConstExtrap, ConstExtrap> bd;
    typedef LinearExtrap bd;
};
std::string Variable<0>::name = "uwind";
template<>
struct Variable<1>{
    typedef Vector<2> Element_t;
    Variable(){};
    Variable(Interval<1> cx, Interval<1> si){
        cell.initialize(cx);
        wall.initialize(si);
    }
    static std::string name;
    static Array<1, Element_t> cell;
    static Array<1, Element_t> wall;
    typedef Int2Type<1> ctag;
    typedef Int2Type<0> ptag;
    //typedef Boundary<FixedConst<Element_t, 1>, Mirror, ConstExtrap, ConstExtrap> bd;
    typedef ConstExtrap bd;
};
std::string Variable<1>::name = "vwind";

/* explicit definition for static variables */
Array<1, Variable<0>::Element_t> Variable<0>::cell;
Array<1, Variable<0>::Element_t> Variable<0>::wall;
Array<1, Variable<1>::Element_t> Variable<1>::cell;
Array<1, Variable<1>::Element_t> Variable<1>::wall;

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
