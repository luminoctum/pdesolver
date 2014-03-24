#ifndef STENCIL
#define STENCIL
struct ForwardDifference{
    static double dx;
    ForwardDifference(){}
    template<class A> inline typename A::Element_t
        operator()(const A& a, int i) const{
            return (a.read(i + 1) - a.read(i)) / dx;
        }
    inline int lowerExtent(int) const { return 0;}
    inline int upperExtent(int) const { return 1;}
};
double ForwardDifference::dx = 1.;

struct ForwardDifferenceX{
    static double dx;
    ForwardDifferenceX(){}
    template<class A> inline typename A::Element_t
        operator()(const A& a, int i, int j) const{
            return (a.read(i + 1, j) - a.read(i, j)) / dx;
        }
    inline int lowerExtent(int) const { return 0;}
    inline int upperExtent(int d) const { return 1 - d;}
};
double ForwardDifferenceX::dx = 1.;

struct ForwardDifferenceY{
    static double dy;
    ForwardDifferenceY(){}
    template<class A> inline typename A::Element_t
        operator()(const A& a, int i, int j) const{
            return (a.read(i, j + 1) - a.read(i, j)) / dy;
        }
    inline int lowerExtent(int) const { return 0;}
    inline int upperExtent(int d) const { return d;}
};
double ForwardDifferenceY::dy = 1.;

struct CenterDifferenceX{
    static double dx;
    CenterDifferenceX(){}
    template<class A> inline typename A::Element_t
        operator()(const A& a, int i, int j) const{
            return (a.read(i + 1, j) - a.read(i - 1, j)) / (2. * dx);
        }
    inline int lowerExtent(int d) const { return 1 - d;}
    inline int upperExtent(int d) const { return 1 - d;}
};
double CenterDifferenceX::dx = 1.;

struct CenterDifferenceY{
    static double dy;
    CenterDifferenceY(){}
    template<class A> inline typename A::Element_t
        operator()(const A& a, int i, int j) const{
            return (a.read(i, j + 1) - a.read(i, j - 1)) / (2. * dy);
        }
    inline int lowerExtent(int d) const { return d;}
    inline int upperExtent(int d) const { return d;}
};
double CenterDifferenceY::dy = 1.;

struct Difference2X{
    static double dx;
    Difference2X(){}
    template<class A> inline typename A::Element_t
        operator()(const A& a, int i, int j) const{
            return (a.read(i + 1, j) - 2. * a.read(i, j) + a.read(i - 1, j)) / (dx * dx);
        }
    inline int lowerExtent(int d) const { return 1 - d;}
    inline int upperExtent(int d) const { return 1 - d;}
};
double Difference2X::dx = 1.;

struct Difference2Y{
    static double dy;
    Difference2Y(){}
    template<class A> inline typename A::Element_t
        operator()(const A& a, int i, int j) const{
            return (a.read(i, j + 1) - 2. * a.read(i, j) + a.read(i, j -1)) / (dy * dy);
        }
    inline int lowerExtent(int d) const { return d;}
    inline int upperExtent(int d) const { return d;}
};
double Difference2Y::dy = 1.;

// Though it is not a stencil, assume 2d
// assume the initial value is zero
struct ForwardIntegralY{
    static double dy;
    ForwardIntegralY(){}
    template<class ET, class EG> 
    inline Array<2, ET> operator()(const Array<2, ET, EG>& a, const Interval<2>& cij) const{
        Interval<1> ci = cij[0];
        int jfirst = cij[1].first();
        int jlast = cij[1].last();
        Array<2> result(cij);
        result(ci, jfirst) = 0.;
        for (int j = jfirst + 1; j <= jlast; j++)
            result(ci, j) = result(ci, j - 1) + 0.5 * dy * (a(ci, j - 1) + a(ci, j));
        return result;
    }
};
double ForwardIntegralY::dy = 1.;

struct BackwardIntegralY{
    static double dy;
    BackwardIntegralY(){}
    template<class ET, class EG> 
    inline Array<2, ET> operator()(const Array<2, ET, EG>& a, const Interval<2>& cij) const{
        Interval<1> ci = cij[0];
        int jfirst = cij[1].first();
        int jlast = cij[1].last();
        Array<2> result(cij);
        result(ci, jlast) = 0.;
        for (int j = jlast - 1; j >= jfirst; j--)
            result(ci, j) = result(ci, j + 1) - 0.5 * dy * (a(ci, j) + a(ci, j + 1));
        return result;
    }
};
double BackwardIntegralY::dy = 1.;

struct Average{
    Average(){}
    template<class A> inline typename A::Element_t
        operator()(const A& a, int i) const{
            return 0.5 * (a.read(i) + a.read(i + 1));
        }
    inline int lowerExtent(int) const { return 0;}
    inline int upperExtent(int) const { return 1;}
};

struct Curvature{
    static double dx, dy;
    Curvature(){};
    template<class A> inline typename A::Element_t
    operator()(const A&a, int i, int j) const{
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
double Curvature::dx = 1.;
double Curvature::dy = 1.;
#endif
