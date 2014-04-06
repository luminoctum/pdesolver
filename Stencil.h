#ifndef STENCIL
#define STENCIL
#include "Pooma/Arrays.h"
enum { DimX, DimY, DimZ, D1, D2, O1, O2, Forward, Backward, Trapz, Staggered};

struct FiniteDifferenceBase{
    static double dx, dy;
};
double FiniteDifferenceBase::dx = 1.;
double FiniteDifferenceBase::dy = 1.;

/* Finite Difference */
template<int drv, int order, int dim, int offset = 0> 
struct FiniteDifference {};

template<>
struct FiniteDifference<D1, O1, DimX> : FiniteDifferenceBase{
    template<class A> inline typename A::Element_t
        operator()(const A& a, int i) const{
            return (a.read(i + 1) - a.read(i)) / dx;
        }
    template<class A> inline typename A::Element_t
        operator()(const A& a, int i, int j) const{
            return (a.read(i + 1, j) - a.read(i, j)) / dx;
        }
    inline int lowerExtent(int d) const { return 0;}
    inline int upperExtent(int d) const { return 1 - d;}
};

template<>
struct FiniteDifference<D1, O1, DimY> : FiniteDifferenceBase{
    template<class A> inline typename A::Element_t
        operator()(const A& a, int i, int j) const{
            return (a.read(i, j + 1) - a.read(i, j)) / dy;
        }
    inline int lowerExtent(int) const { return 0;}
    inline int upperExtent(int d) const { return d;}
};

template<>
struct FiniteDifference<D1, O2, DimX> : FiniteDifferenceBase{
    template<class A> inline typename A::Element_t
        operator()(const A& a, int i) const{
            return (a.read(i + 1) - a.read(i - 1)) / (2. * dx);
        }
    template<class A> inline typename A::Element_t
        operator()(const A& a, int i, int j) const{
            return (a.read(i + 1, j) - a.read(i - 1, j)) / (2. * dx);
        }
    inline int lowerExtent(int d) const { return 1 - d;}
    inline int upperExtent(int d) const { return 1 - d;}
};

template<>
struct FiniteDifference<D1, O2, DimY> : FiniteDifferenceBase{
    template<class A> inline typename A::Element_t
        operator()(const A& a, int i, int j) const{
            return (a.read(i, j + 1) - a.read(i, j - 1)) / (2. * dy);
        }
    inline int lowerExtent(int d) const { return d;}
    inline int upperExtent(int d) const { return d;}
};

template<>
struct FiniteDifference<D2, O2, DimX> : FiniteDifferenceBase{
    template<class A> inline typename A::Element_t
        operator()(const A& a, int i, int j) const{
            return (a.read(i + 1, j) - 2. * a.read(i, j) + a.read(i - 1, j)) / (dx * dx);
        }
    inline int lowerExtent(int d) const { return 1 - d;}
    inline int upperExtent(int d) const { return 1 - d;}
};

template<>
struct FiniteDifference<D2, O2, DimY> : FiniteDifferenceBase{
    template<class A> inline typename A::Element_t
        operator()(const A& a, int i, int j) const{
            return (a.read(i, j + 1) - 2. * a.read(i, j) + a.read(i, j -1)) / (dy * dy);
        }
    inline int lowerExtent(int d) const { return d;}
    inline int upperExtent(int d) const { return d;}
};

// Though currently it is not a stencil yet, it is still placed here
template <int dir, int type, int dim> struct Integrate{};

template<>
struct Integrate<Forward, Trapz, DimY> : FiniteDifferenceBase{
    template<class A, class B> void inline
    operator()(const A& a, const B& result, const Interval<2>& cij) const{
        Interval<1> ci = cij[0], cj = cij[1];
        result(ci, cj.first()) = 0.;
        for (int j = cj.first() + 1; j <= cj.last(); j++)
            result(ci, j) = result(ci, j - 1) + 0.5 * dy * (a(ci, j - 1) + a(ci, j));
    }
};

template<>
struct Integrate<Backward, Trapz, DimY> : FiniteDifferenceBase{
    template<class A, class B> void inline
    operator()(const A& a, const B& result, const Interval<2>& cij) const{
        Interval<1> ci = cij[0], cj = cij[1];
        result(ci, cj.last()) = 0.;
        for (int j = cj.last() - 1; j >= cj.first(); j--)
            result(ci, j) = result(ci, j + 1) - 0.5 * dy * (a(ci, j) + a(ci, j + 1));
    }
};

template<>
struct Integrate<Backward, Staggered, DimY> : FiniteDifferenceBase{
    template<class A, class B> void inline
    operator()(const A& a, const B& result, const Interval<2>& cij) const{
        Interval<1> ci = cij[0], cj = cij[1];
        result(ci, cj.last()) = 0.;
        for (int j = cj.last(); j >= cj.first(); j--){
            result(ci, j - 1) = result(ci, j) - dy * a(ci, j);
        }
    }
    template<class A, class T, class G> void inline
    operator()(const A& a, const Array<2, Vector<2, T>, G>& result, const Interval<2>& cij) const{
        Interval<1> ci = cij[0], cj = cij[1];
        result(ci, cj.last() + 1).comp(0) = 0.;
        result(ci, cj.last()).comp(1) = 0.;
        for (int j = cj.last(); j >= cj.first(); j--){
            result(ci, j).comp(0) = result(ci, j).comp(1) - dy * a(ci, j);
            result(ci, j - 1).comp(1) = result(ci, j).comp(0);
        }
    }
};

/* miscellaneous */
struct Average{
    template<class A> inline typename A::Element_t
        operator()(const A& a, int i) const{
            return 0.5 * (a.read(i) + a.read(i + 1));
        }
    inline int lowerExtent(int) const { return 0;}
    inline int upperExtent(int) const { return 1;}
};

struct Laplace{
    template<class A> inline typename A::Element_t
        operator()(const A& a, int i, int j) const{
            return (a.read(i, j + 1) + a.read(i + 1, j) + a.read(i, j - 1) + a.read(i - 1, j) - 4. * a.read(i, j));
        }
    inline int lowerExtent(int) const { return 1;}
    inline int upperExtent(int) const { return 1;}
};

struct Curvature : FiniteDifferenceBase{
    template<class A> inline typename A::Element_t
    operator()(const A& a, int i, int j) const{
        typename A::Element_t ax, ay, axx, axy, ayy;
        ax = (a.read(i + 1, j) - a.read(i - 1, j)) / (2. * dx);
        ay = (a.read(i, j + 1) - a.read(i, j - 1)) / (2. * dy);
        axx = (a.read(i - 1, j) + a.read(i + 1, j) - 2. * a.read(i, j)) / (dx * dx);
        ayy = (a.read(i, j - 1) + a.read(i, j + 1) - 2. * a.read(i, j)) / (dy * dy);
        axy = (a.read(i + 1, j + 1) + a.read(i - 1, j - 1) - a.read(i - 1, j + 1) - a.read(i - 1, j - 1)) / (4. * dx * dy);
        return (ay*ay*axx - 2.*ax*ay*axy + ax*ax*ayy) / pow(ax*ax+ay*ay, 1.5);
    }
    inline int lowerExtent(int) const {return 1;}
    inline int upperExtent(int) const {return 1;}
};

struct ElementSum{
    template<class A> inline typename A::Element_t::Element_t
    operator()(const A& a, int i) const {return sum(a.read(i));}
    template<class A> inline typename A::Element_t::Element_t
    operator()(const A& a, int i, int j) const {return sum(a.read(i, j));}
    inline int lowerExtent(int) const {return 0;}
    inline int upperExtent(int) const {return 0;}
};

struct ElementAverage{
    template<class A> inline typename A::Element_t::Element_t
    operator()(const A& a, int i) const {
        return 0.5 * (a.read(i)(0) + a.read(i)(1));
    }
    template<class A> inline typename A::Element_t::Element_t
    operator()(const A& a, int i, int j) const {
        return 0.5 * (a.read(i, j)(0) + a.read(i, j)(1));
    }
    inline int lowerExtent(int) const {return 0;}
    inline int upperExtent(int) const {return 0;}
};


/* Partial specialize Return of Functor */
template<class T>
struct FunctorResult<ElementSum, T>{
    typedef typename T::Element_t Type_t;
};

template<class T>
struct FunctorResult<ElementAverage, T>{
    typedef typename T::Element_t Type_t;
};

#endif
