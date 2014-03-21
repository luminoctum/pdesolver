#ifndef SHALLOWWATER
#define SHALLOWWATER
#include "Boundary.h"

template<int DIM, class ET>
class ShallowWater{
    typedef ET Element_t;
    enum {dim = DIM};
};

/* 1D shallow water equation */
template<class ET>
class ShallowWater<1, ET>{
    /* FIXME: may consider to add a takeoffcheck to ensure all parameters are properly set.
     * Also you may want to add a function to check values of a variable like cfl, fliud depth etc*/
    enum {dim = 1};
    typedef ET Element_t;
    typedef Vector<2, Element_t> Vector_t;
    typedef Array<dim, Vector_t> Array_t;
    typedef Interval<dim> Interval_t;
    typedef Tensor<2, Element_t> Matrix_t;
    typedef FixedConst<Vector_t> Bd;
protected:
    Array_t ua;
    Array<dim, Element_t> xa;
    inline void fixBoundary(const Array_t& ua, Interval_t ci){
        Bd::fix(ua, ci);
    }
public:
    static inline Vector_t cellFlux(Vector_t u){ 
        return Vector_t(u(1), u(1) * u(1) / u(0) + 0.5 * u(0) * u(0));
    }
    static inline Matrix_t dFluxdu(Vector_t ua){
        Matrix_t dfdu;
        dfdu(0, 0) = 0; dfdu(0, 1) = 1;
        dfdu(1, 0) = -pow(ua(1) / ua(0), 2) + ua(0);
        dfdu(1, 1) = 2. * ua(1) / ua(0);
        return dfdu;
    }
    static inline Matrix_t roeMatrix(Vector_t ua, Vector_t ub){
        Matrix_t roe;
        Element_t h1, h2, hv1, hv2, v;
        h1 = ua(0); hv1 = ua(1);
        h2 = ub(0); hv2 = ub(1);
        v = (hv1 / sqrt(h1) + hv2 / sqrt(h2))/(sqrt(h1) + sqrt(h2));
        roe(0, 0) = 0.; roe(0, 1) = 1.;
        roe(1, 0) = 0.5 * (h1 + h2) - v * v;
        roe(1, 1) = 2. * v;
        return roe;
    }
};


#endif
