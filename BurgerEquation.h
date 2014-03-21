#ifndef BURGEREQUATION
#define BURGEREQUATION
#include "utils.h"

template<int DIM, class ET>
struct BurgerEquation{
    typedef ET Element_t;
    enum {dim = DIM};
};

template<class ET>
struct BurgerEquation<1, ET>{
    typedef ET Element_t;
    enum {dim = 1};
    static inline Element_t cellFlux(Element_t u){ 
        return 0.5 * u * u;
    }
    /* assuming a < b */
    static inline Element_t minFlux(Element_t a, Element_t b){
        return (a > 0 ? cellFlux(a) : (b < 0 ? cellFlux(b) : 0 ));
    }
    /* assuming a < b */
    static inline Element_t maxFlux(Element_t a, Element_t b){
        return (a > 0 ? cellFlux(b) : (b < 0 ? cellFlux(a) : _max_(cellFlux(a), cellFlux(b))));
    }
};

#endif
