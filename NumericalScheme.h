#ifndef NUMERICALSCHEME
#define NUMERICALSCHEME
#include "Stencil.h"

/* 1d EnoScheme */
template<int O>
class EnoScheme{
#define _maxorder_ 5
protected:
    static double mpoly[_maxorder_][_maxorder_];
    static Array<2> enoc, polyc;
public:
    EnoScheme(){ initialize(); }
    static void initialize(){
        enoc.initialize(Interval<1>(-1, O - 1), Interval<1>(0, O - 1));
        double coeff1, coeff2, coeff3, result;
        for(int r = -1; r < O; r++) for (int j = 0; j < O; j++){
            result = 0.;
            for (int m = j + 1; m < O + 1; m++){
                coeff2 = 0.; coeff3 = 1.;
                for (int l = 0; l < O + 1; l++){
                    if (l == m) continue;
                    coeff1 = 1.;
                    for (int q = 0; q < O + 1; q++){
                        if (q == m || q == l) continue;
                        coeff1 *= r - q + 1;
                    }
                    coeff2 += coeff1;
                    coeff3 *= m - l;
                }
                result += coeff2 / coeff3;
            }
            enoc(r, j) = result;
        }
        polyc.initialize(Interval<1>(0, _maxorder_ - 1), Interval<1>(0, _maxorder_ - 1));
        for (int i = 0; i < _maxorder_; i++) for (int j = 0; j < _maxorder_; j++)
            polyc(i, j) = mpoly[i][j];
    }

    /* scalar reconstruction */
    template<class A, class B> void
    operator()(const A& ua, const B& result, int i, int offset = 0){
        int imin, imax;
        imin = i; imax = i;
        for (int m = 1; m < O; m++){
            if (fabs(sum(
                    ua(Interval<1>(imin - 1, imax)) 
                    * polyc(m, Interval<1>(_maxorder_ - 1 - m, _maxorder_ - 1))
                )) < 
                fabs(sum(
                    ua(Interval<1>(imin, imax + 1)) 
                    * polyc(m, Interval<1>(_maxorder_ - 1 - m, _maxorder_ - 1))
                ))) imin--;
            else imax++;
        }
        result(i + offset) = 0.;
        for (int m = 0; m < O; m++){
            result(i + offset)(0) += enoc(i - imin - 1, m) * ua(imin + m);
            result(i + offset)(1) += enoc(i - imin, m) * ua(imin + m);
        }
    }

    /* component-wise reconstruction */
    template<int S, class E, class G, class B> inline void
    operator()(const Array<1, Vector<S, E>, G>& ua, const B& result, int i, int offset = 0){
        int imin, imax;
        result(i + offset) = 0.;
        for (int s = 0; s < S; s++){
            imin = i; imax = i;
            for (int m = 1; m < O; m++){
                if (fabs(sum(
                        ua(Interval<1>(imin - 1, imax)).comp(s)
                        * polyc(m, Interval<1>(_maxorder_ - 1 - m, _maxorder_ - 1))
                    )) < 
                    fabs(sum(
                        ua(Interval<1>(imin, imax + 1)).comp(s)
                        * polyc(m, Interval<1>(_maxorder_ - 1 - m, _maxorder_ - 1))
                    ))) imin--;
                else imax++;
            }
            for (int m = 0; m < O; m++){
                result(i + offset)(0)(s) += enoc(i - imin - 1, m) * ua(imin + m)(s);
                result(i + offset)(1)(s) += enoc(i - imin, m) * ua(imin + m)(s);
            }
        }
    }

    /* characteristic reconstruction ... TODO */

    template<class A, class B> inline void
    operator()(const A& ua, const B& result, Interval<1> ci){
        for (int i = ci.first() - 1; i <= ci.last() + 1; i++) operator()(ua, result, i);
    }

    /* 2D construction */
    template<class A, class B, class C> inline void
    operator()(const A& ua, const B& resultX, const C& resultY, int i, int j){
        // wall value is one grid wider than interier cell value
        operator()(ua(AllDomain<1>(), j), resultX(AllDomain<1>(), j), i - ua.domain()[0].first(), ua.domain()[0].first() + 1);
        operator()(ua(i, AllDomain<1>()), resultY(i, AllDomain<1>()), j - ua.domain()[1].first(), ua.domain()[1].first() + 1);
    }

    template<class A, class B, class C> inline void
    operator()(const A& ua, const B& resultX, const C& resultY, Interval<2> cij){
        //std::cout << " calculating ... (" << i << ", " << j << ") " << std::endl;
        int i, j;
        #pragma omp for private(i)
        for (j = cij[1].first(); j <= cij[1].last(); j++)
            for (i = cij[0].first() - 1; i <= cij[0].last() + 1; i++)
                operator()(ua(AllDomain<1>(), j), resultX(AllDomain<1>(), j), i - ua.domain()[0].first(), ua.domain()[0].first() + 1);
        #pragma omp for private(j)
        for (i = cij[0].first(); i <= cij[0].last(); i++)
            for (j = cij[1].first() - 1; j <= cij[1].last() + 1; j++)
                operator()(ua(i, AllDomain<1>()), resultY(i, AllDomain<1>()), j - ua.domain()[1].first(), ua.domain()[1].first() + 1);
        //construct(ua, resultX, resultY, i, j);
    }
};
template<int O> double EnoScheme<O>::mpoly[_maxorder_][_maxorder_] = {
    {0., 0., 0., 0., 1.},
    {0., 0., 0., -1., 1.},
    {0., 0., 1., -2., 1.},
    {0., -1., 3., -3., 1.},
    {1., -4., 6., -4., 1.}
};
template<int O> Array<2> EnoScheme<O>::enoc;
template<int O> Array<2> EnoScheme<O>::polyc;

#undef _maxorder_

/* Godunov Flux */
template<class T, int Dim = DimX>
struct Godunov{
    Array<2, Vector<2, T> > *xwind;
    template<class A> inline typename A::Element_t::Element_t
    operator()(const A& a, int i) const {
        if (xwind->read(i - 1)(1) < 0 && xwind->read(i)(0) > 0) return typename A::Element_t::Element_t(0.);
        if (xwind->read(i - 1)(1) + xwind->read(i)(0) > 0) return a.read(i - 1)(1) * xwind->read(i - 1)(1);
        else return a.read(i)(0) * xwind->read(i)(0); 
    }
    template<class A> inline typename A::Element_t::Element_t
    operator()(const A& a, int i, int j) const {
        if (xwind->read(i - 1, j)(1) < 0 && xwind->read(i, j)(0) > 0) return typename A::Element_t::Element_t(0.);
        if (xwind->read(i - 1, j)(1) + xwind->read(i, j)(0) > 0) return a.read(i - 1, j)(1) * xwind->read(i - 1, j)(1);
        else return a.read(i, j)(0) * xwind->read(i, j)(0); 
    }
    inline int lowerExtent(int d) const {return 1 - d;}
    inline int upperExtent(int d) const {return 0;}
};

template<class T>
struct Godunov<T, DimY>{
    Array<2, Vector<2, T> > *ywind;
    template<class A> inline typename A::Element_t::Element_t
    operator()(const A& a, int i, int j) const {
        if (ywind->read(i, j - 1)(1) < 0 && ywind->read(i, j)(0) > 0) return typename A::Element_t::Element_t(0.);
        if (ywind->read(i, j - 1)(1) + ywind->read(i, j)(0) > 0) return a.read(i, j - 1)(1) * ywind->read(i, j - 1)(1);
        else return a.read(i, j)(0) * ywind->read(i, j)(0); 
    }
    inline int lowerExtent(int d) const {return d;}
    inline int upperExtent(int d) const {return 0;}
};

template<class E, int Dim, class T>
struct FunctorResult<Godunov<E, Dim>, T>{
    typedef typename T::Element_t Type_t;
};

/* Euler Forward */
struct Euler{};

#endif
