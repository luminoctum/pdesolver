#ifndef STENCIL
#define STENCIL
enum { DimX, DimY, DimZ, Forward, Backward, Trapz, Staggered};

struct FiniteDifferenceBase{
    static double dx, dy;
};
double FiniteDifferenceBase::dx = 1.;
double FiniteDifferenceBase::dy = 1.;

/* Finite Difference */
template<int drv, int order, int dim, int offset = 0> 
struct FiniteDifference {};

template<>
struct FiniteDifference<1, 1, DimX> : FiniteDifferenceBase{
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
struct FiniteDifference<1, 1, DimY> : FiniteDifferenceBase{
    template<class A> inline typename A::Element_t
        operator()(const A& a, int i, int j) const{
            return (a.read(i, j + 1) - a.read(i, j)) / dy;
        }
    inline int lowerExtent(int) const { return 0;}
    inline int upperExtent(int d) const { return d;}
};

template<>
struct FiniteDifference<1, 2, DimX> : FiniteDifferenceBase{
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
struct FiniteDifference<1, 2, DimY> : FiniteDifferenceBase{
    template<class A> inline typename A::Element_t
        operator()(const A& a, int i, int j) const{
            return (a.read(i, j + 1) - a.read(i, j - 1)) / (2. * dy);
        }
    inline int lowerExtent(int d) const { return d;}
    inline int upperExtent(int d) const { return d;}
};

template<>
struct FiniteDifference<2, 2, DimX> : FiniteDifferenceBase{
    template<class A> inline typename A::Element_t
        operator()(const A& a, int i, int j) const{
            return (a.read(i + 1, j) - 2. * a.read(i, j) + a.read(i - 1, j)) / (dx * dx);
        }
    inline int lowerExtent(int d) const { return 1 - d;}
    inline int upperExtent(int d) const { return 1 - d;}
};

template<>
struct FiniteDifference<2, 2, DimY> : FiniteDifferenceBase{
    template<class A> inline typename A::Element_t
        operator()(const A& a, int i, int j) const{
            return (a.read(i, j + 1) - 2. * a.read(i, j) + a.read(i, j -1)) / (dy * dy);
        }
    inline int lowerExtent(int d) const { return d;}
    inline int upperExtent(int d) const { return d;}
};

// Though currently it is not a stencil, it is also placed here
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
        result(ci, cj.last() - 1) = - 2. * dy * a(ci, cj.last());
        for (int j = cj.last() - 1; j >= cj.first(); j--){
            result(ci, j - 1) = result(ci, j + 1) - 2. * dy * a(ci, j);
        }
    }
};

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

/* misc */

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


/* Partial specialize Return of Functor */
template<class T>
struct FunctorResult<ElementSum, T>{
    typedef typename T::Element_t Type_t;
};

template<class E, int Dim, class T>
struct FunctorResult<Godunov<E, Dim>, T>{
    typedef typename T::Element_t Type_t;
};

#endif
