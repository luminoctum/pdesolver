#ifndef GODUNOV
#define GODUNOV
#include "Pooma/Arrays.h"
/* assume DIM = 1 */
template<int DIM>
struct Godunov{
    /* wallFlux for scalar */
    template<class A, class W, class I> void
    static inline wallFlux(A, const W& windm, const W& windp, I si, Int2Type<0>){
        int ifirst = si.first();
        int ilast = si.last();
        for (int i = ifirst; i <= ilast; i++){
            if (windm(i) + windp(i) > 0)
                A::flux(i) = A::wallm(i) * windm(i);
            else
                A::flux(i) = A::wallp(i) * windp(i);
        }
    }

    /* wallFlux fot decoupled vector */
    template<class A, class W, class I> void
    static inline wallFlux(A, const W& windm, const W& windp, I si, Int2Type<1>){
        wallFlux(A(), windm, windp, si, Int2Type<0>());
    }

    /* wallFlux for coupled vector */
    template<class A, class I> void
    static inline wallFlux(A, I si, Int2Type<2>){
        int ifirst = si.first();
        int ilast = si.last();
        for (int i = ifirst; i <= ilast; i++){
            A::flux(i) = 0.5 * (A::wallm(i) + A::wallp(i))
                - 0.5 * dot(entropyFix(A(), A::wallm(i), A::wallp(i)), A::wallp(i) - A::wallm(i));
        }
    }

    template<class A, int SIZE, class ET> Tensor<SIZE, ET>
    static inline entropyFix(A, const Vector<SIZE, ET>& ua, const Vector<SIZE, ET>& ub){
        Tensor<SIZE, ET> lp1, D1, rp1, lp2, D2, rp2, lp, D, rp;
        ET beta;
        EigenSystem<SIZE>::solve(A::dFluxdu(ua), lp1, D1, rp1);
        EigenSystem<SIZE>::solve(A::dFluxdu(ub), lp2, D2, rp2);
        EigenSystem<SIZE>::solve(A::roeMatrix(ua, ub), lp, D, rp);
        for (int s = 0; s < SIZE; s++)
            if ((D1(s, s) < 0.) and (D2(s, s) > 0.)){
                beta = (D2(s, s) - D(s, s))/(D2(s, s) - D1(s, s));
                D(s, s) = (1 - beta) * D2(s, s) - beta * D1(s, s);
            } else{
                D(s, s) = fabs(D(s, s));
            }
        return dot(dot(rp, D), lp);
    }

    template<class A, int SIZE, class ET> Tensor<SIZE, ET>
    static inline absRoe(A, const Vector<SIZE, ET>& ua, const Vector<SIZE, ET>& ub){
        Tensor<SIZE, ET> lp, D, rp;
        EigenSystem<SIZE>::solve(A::roeMatrix(ua, ub), lp, D, rp);
        return dot(dot(rp, fabs(D)), lp);
    }
};
#endif
