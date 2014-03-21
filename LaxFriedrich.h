#ifndef LAXFRIEDRICH
#define LAXFRIEDRICH
#include "Pooma/Arrays.h"
template<int DIM, class ET, class SYS>
class LaxFriedrich{
    typedef ET Element_t;
    typedef Array<DIM, ET> Array_t;
    typedef Interval<DIM> Interval_t;
private:
    double cfl;
public:
    LaxFriedrich(){}
    LaxFriedrich(double _cfl) : cfl(_cfl){}
    template<class EG1> Array_t wallFlux(const Array<DIM, Element_t, EG1>& ua, const Interval_t& ci){
        Array_t result(ci.length() + 1);
        for (int i = ci.first(); i <= ci.last() + 1; i++){
            result(i) = 0.5 * (SYS::cellFlux(ua(i - 1)), SYS::cellFlux(ua(i))) 
                - 0.5 * cfl * (SYS::cellFlux(ua(i)) - SYS::cellFlux(ua(i - 1)));
        }
        return result;
    }
};
#endif
