#ifndef GODUNOV
#define GODUNOV
#include "Pooma/Arrays.h"
/* assume DIM = 1 */
template<int DIM, class SYS>
struct Godunov{
    template<int SIZE, class ET> Tensor<SIZE, ET>
    static inline entropyFix(const Vector<SIZE, ET>& ua, const Vector<SIZE, ET>& ub){
        Tensor<SIZE, ET> lp1, D1, rp1, lp2, D2, rp2, lp, D, rp;
        ET beta;
        EigenSystem<SIZE>::solve(SYS::dFluxdu(ua), lp1, D1, rp1);
        EigenSystem<SIZE>::solve(SYS::dFluxdu(ub), lp2, D2, rp2);
        EigenSystem<SIZE>::solve(SYS::roeMatrix(ua, ub), lp, D, rp);
        for (int s = 0; s < SIZE; s++)
            if ((D1(s, s) < 0.) and (D2(s, s) > 0.)){
                beta = (D2(s, s) - D(s, s))/(D2(s, s) - D1(s, s));
                D(s, s) = (1 - beta) * D2(s, s) - beta * D1(s, s);
            } else{
                D(s, s) = fabs(D(s, s));
            }
        return dot(dot(rp, D), lp);
    }

    template<int SIZE, class ET> Tensor<SIZE, ET>
    static inline absRoe(const Vector<SIZE, ET>& ua, const Vector<SIZE, ET>& ub){
        Tensor<SIZE, ET> lp, D, rp;
        EigenSystem<SIZE>::solve(SYS::roeMatrix(ua, ub), lp, D, rp);
        return dot(dot(rp, fabs(D)), lp);
    }

    template<class ET, class EG> Array<DIM, ET>
    static inline wallFlux(const Array<DIM, ET, EG>& ua, const Interval<DIM>& ci){
        Interval<DIM> si(ci.first(), ci.last() + 1);
        return wallFlux(ua(si - 1), ua(si));
    }

    /* wallFlux for scalar */
    template<class ET, class EG1, class EG2> Array<DIM, ET>
    static inline wallFlux(const Array<DIM, ET, EG1>& ua, const Array<DIM, ET, EG2>& ub){
        int ifirst = ua.domain().first();
        int ilast = ua.domain().last();
        Array<DIM, ET> result(ua.domain());
        for (int i = ifirst; i <= ilast; i++)
            if (ua(i) < ub(i)) result(i) = SYS::minFlux(ua(i), ub(i));
            else result(i) = SYS::maxFlux(ub(i), ua(i));
        return result;
    }

    /* wallFlux for vector */
    template<int SIZE, class ET, class EG1, class EG2> Array<DIM, Vector<SIZE, ET> >
    static inline wallFlux(const Array<DIM, Vector<SIZE, ET>, EG1>& ua, const Array<DIM, Vector<SIZE, ET>, EG2>& ub){
        int ifirst = ua.domain().first();
        int ilast = ua.domain().last();
        Array<DIM, Vector<SIZE, ET> > result(ua.domain());
        for (int i = ifirst; i <= ilast; i++){
            result(i) = 0.5 * (SYS::cellFlux(ua(i)) + SYS::cellFlux(ub(i))) 
                - 0.5 * dot(entropyFix(ua(i), ub(i)), ub(i) - ua(i));
                //- 0.5 * dot(absRoe(ua(i), ub(i)), ub(i) - ua(i));
        }
        return result;
    }
};
#endif
