#ifndef ENOSCHEME
#define ENOSCHEME
#include "Pooma/Arrays.h"

/* 1d EnoScheme */
template<int O>
class EnoScheme{
#define _maxorder_ 5
protected:
    static double mpoly[_maxorder_][_maxorder_];
    static Array<2> enoc, polyc;
public:
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

    EnoScheme(){}
    inline int lowerExtent(int) const {return O - 1;}
    inline int upperExtent(int) const {return O - 1;}
    /* scalar reconstruction */
    template<class A> Vector<2, typename A::Element_t>
    operator()(const A& ua, int i) const {
        int imin, imax;
        Vector<2, typename A::Element_t> result;
        result = 0.;
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
        for (int m = 0; m < O; m++){
            result(0) += enoc(i - imin - 1, m) * ua(imin + m);
            result(1) += enoc(i - imin, m) * ua(imin + m);
        }
        return result;
    }

    /* component-wise reconstruction */
    template<int S, class E, class G> Vector<2, Vector<S, E> >
    inline operator()(const Array<1, Vector<S, E>, G>& ua, int i) const {
        // relies on actual domain starts from zero
        Vector<2, E> element;
        Vector<2, Vector<S, E> > result;
        for (int s = 0; s < S; s++){
            element = operator()(ua.comp(s), i);
            result(0)(s) = element(0);
            result(1)(s) = element(1);
        }
        return result;
    }

    /* characteristic reconstruction ... TODO */
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

template<int O>
class EnoSchemeX : public EnoScheme<O>{
public:
    inline int lowerExtent(int d) const {return O - 1 + (1 - O) * d;}
    inline int upperExtent(int d) const {return O - 1 + (1 - O) * d;}
    // overload parent function
    using EnoScheme<O>::operator();
    template<class A> Vector<2, typename A::Element_t>
    inline operator()(const A& ua, int i, int j) const {
        // relies on actual domain starts from zero
        Interval<1> cx(ua.domain()[0]);
        return operator()(ua(cx, j), i);
    }
};

template<int O>
class EnoSchemeY : public EnoScheme<O>{
public:
    inline int lowerExtent(int d) const {return (O - 1) * d;}
    inline int upperExtent(int d) const {return (O - 1) * d;}
    // overload parent function
    using EnoScheme<O>::operator();
    template<class A> Vector<2, typename A::Element_t>
    inline operator()(const A& ua, int i, int j) const {
        // relies on actual domain starts from zero
        Interval<1> cy(ua.domain()[1]);
        return operator()(ua(i, cy), j);
    }
};
#undef _maxorder_

/* Partial specialize Return of Functor */
template<int O, class T>
struct FunctorResult<EnoScheme<O>, T>{
    typedef Vector<2, T> Type_t;
};

template<int O, class T>
struct FunctorResult<EnoSchemeX<O>, T>{
    typedef Vector<2, T> Type_t;
};

template<int O, class T>
struct FunctorResult<EnoSchemeY<O>, T>{
    typedef Vector<2, T> Type_t;
};

#endif
