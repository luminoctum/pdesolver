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
    Variable(){
        name = "var_a";
        cell.initialize(10);
        wall.initialize(11);
    }
    static std::string name;
    static Array<1> cell;
    static Array<1> wall;
    typedef Int2Type<0> ctag;
    typedef Int2Type<1> ptag;
};

template<>
struct Variable<1>{
    Variable(){
        cell.initialize(10);
        wall.initialize(11);
    }
    static std::string name;
    static Array<1> cell;
    static Array<1> wall;
    static Int2Type<1> ctag;
    static Int2Type<0> ptag;
};

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

// explicit definition for static variables;
std::string Variable<0>::name;
Array<1> Variable<0>::cell;
Array<1> Variable<0>::wall;
std::string Variable<1>::name;
Array<1> Variable<1>::cell;
Array<1> Variable<1>::wall;

#endif
