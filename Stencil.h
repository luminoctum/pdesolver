#ifndef STENCIL
#define STENCIL
struct Difference{
    static double dx;
    Difference(){}
    template<class A> inline typename A::Element_t
        operator()(const A& a, int i) const{
            return (a.read(i + 1) - a.read(i)) / dx;
        }
    inline int lowerExtent(int) const { return 0;}
    inline int upperExtent(int) const { return 1;}
};
double Difference::dx = 1.;

struct DifferenceX{
    static double dx;
    DifferenceX(){}
    template<class A> inline typename A::Element_t
        operator()(const A& a, int i, int j) const{
            return (a.read(i + 1, j) - a.read(i, j)) / dx;
        }
    inline int lowerExtent(int) const { return 0;}
    inline int upperExtent(int d) const { return 1 - d;}
};
double DifferenceX::dx = 1.;

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

struct DifferenceY{
    static double dy;
    DifferenceY(){}
    template<class A> inline typename A::Element_t
        operator()(const A& a, int i, int j) const{
            return (a.read(i, j + 1) - a.read(i, j)) / dy;
        }
    inline int lowerExtent(int) const { return 0;}
    inline int upperExtent(int d) const { return d;}
};
double DifferenceY::dy = 1.;

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
