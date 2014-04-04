#ifndef BOUNDARY
#define BOUNDARY
#include "Pooma/Arrays.h"
struct Periodic{
    template<class A> void fix(A& a, const Interval<1>& ci){
        Interval<1> cx = a.domain();
        Interval<1> lower(ci.first(), cx.last() - ci.length());
        Interval<1> upper(ci.length() + cx.first(), ci.last());
        Interval<1> lowerfix(cx.first(), ci.first() - 1);
        Interval<1> upperfix(ci.last() + 1, cx.last());
        a(lowerfix) = a(upper);
        a(upperfix) = a(lower);
    }
    template<class A> void fix(A a, const Interval<2>& cij){
        Interval<2> cxy = a.domain();
        Interval<1> cj = cij[1];
        Interval<1> cx = cxy[0];
        int xfirst = cxy[0].first(), xlast = cxy[0].last();
        int yfirst = cxy[1].first(), ylast = cxy[1].last();
        int ifirst = cij[0].first(), ilast = cij[0].last();
        int jfirst = cij[1].first(), jlast = cij[1].last();
        for (int i = xfirst; i < ifirst; i++) a(i, cj) = a(i + ilast - ifirst + 1, cj);
        for (int i = xlast; i > ilast; i--) a(i, cj) = a(i - ilast + ifirst - 1, cj);
        for (int j = yfirst; j < jfirst; j++) a(cx, j) = a(cx, j + jlast - jfirst + 1);
        for (int j = ylast; j > jlast; j--) a(cx, j) = a(cx, j - jlast + jfirst - 1);
    }
    friend std::ostream& operator<< (std::ostream &os, const Periodic &other){
        os << "Periodic"; return os;
    }
};

struct LinearExtrap{
    template<class A> void fix(A& a, const Interval<1>& ci){
        Interval<1> cx = a.domain();
        int ifirst = ci.first(), ilast = ci.last();
        int xfirst = cx.first(), xlast = cx.last();
        for (int i = ifirst - 1; i >= xfirst; i--) a(i) = a(i + 1) - a(ifirst + 1) + a(ifirst);
        for (int i = ilast + 1; i <= xlast; i++) a(i) = a(i - 1) + a(ilast) - a(ilast - 1);
    }
    template<class A> void fix(A& a, const Interval<2>& cij){
        Interval<2> cxy = a.domain();
        int xfirst = cxy[0].first(), xlast = cxy[0].last();
        int yfirst = cxy[1].first(), ylast = cxy[1].last();
        int ifirst = cij[0].first(), ilast = cij[0].last();
        int jfirst = cij[1].first(), jlast = cij[1].last();
        for (int i = ifirst - 1; i >= xfirst; i--) 
            a(i, cij[1]) = a(i + 1, cij[1]) - a(ifirst + 1, cij[1]) + a(ifirst, cij[1]);
        for (int i = ilast + 1; i <= xlast; i++)
            a(i, cij[1]) = a(i - 1, cij[1]) + a(ilast, cij[1]) - a(ilast - 1, cij[1]);
        for (int j = jfirst - 1; j >= yfirst; j--) 
            a(cxy[0], j) = a(cxy[0], j + 1) - a(cxy[0], jfirst + 1) + a(cxy[0], jfirst);
        for (int j = jlast + 1; j <= ylast; j++)
            a(cxy[0], j) = a(cxy[0], j - 1) + a(cxy[0], jlast) - a(cxy[0], jlast - 1);
    }
    friend std::ostream& operator<< (std::ostream &os, const LinearExtrap &other){
        os << "LinearExtrap"; return os;
    }
};

struct ConstExtrap{
    template<class A> void fix(A& a, const Interval<1>& ci){
        Interval<1> cx = a.domain();
        int ifirst = ci.first(), ilast = ci.last();
        int xfirst = cx.first(), xlast = cx.last();
        for (int i = ifirst - 1; i >= xfirst; i--) a(i) = a(ifirst);
        for (int i = ilast + 1; i <= xlast; i++) a(i) = a(ilast);
    }
    template<class A> void fix(A& a, const Interval<2>& cij){
        Interval<2> cxy = a.domain();
        int xfirst = cxy[0].first(), xlast = cxy[0].last();
        int yfirst = cxy[1].first(), ylast = cxy[1].last();
        int ifirst = cij[0].first(), ilast = cij[0].last();
        int jfirst = cij[1].first(), jlast = cij[1].last();
        Interval<1> cj = cij[1];
        Interval<1> cx = cxy[0];
        for (int i = xfirst; i < ifirst; i++) a(i, cj) = a(ifirst, cj);
        for (int i = xlast; i > ilast; i--) a(i, cj) = a(ilast, cj);
        for (int j = yfirst; j < jfirst; j++) a(cx, j) = a(cx, jfirst);
        for (int j = ylast; j > jlast; j--) a(cx, j) = a(cx, jlast);
    }
    friend std::ostream& operator<< (std::ostream &os, const ConstExtrap &other){
        os << "ConstExtrap"; return os;
    }
};

template<class T = double>
struct FixedConst{
    typedef T Element_t;
    bool unset;
    Element_t lower, upper;
    Array<1, Element_t> left, right, bottom, top;
    FixedConst(){
        unset = true;
    }
    template<class A> void fix(A& a, const Interval<1>& ci) {
        Interval<1> cx = a.domain();
        int ifirst = ci.first(), ilast = ci.last();
        int xfirst = cx.first(), xlast = cx.last();
        if (unset){ 
            lower = a(ifirst); 
            upper = a(ilast); 
            unset = false;
        }
        for (int i = ifirst; i >= xfirst; i--) a(i) = lower;
        for (int i = ilast; i <= xlast; i++) a(i) = upper;
    }
    template<class A> void fix(A& a, const Interval<2>& cij){
        Interval<2> cxy = a.domain();
        int xfirst = cxy[0].first(), xlast = cxy[0].last();
        int yfirst = cxy[1].first(), ylast = cxy[1].last();
        int ifirst = cij[0].first(), ilast = cij[0].last();
        int jfirst = cij[1].first(), jlast = cij[1].last();
        Interval<1> cj = cij[1];
        Interval<1> cx = cxy[0];
        if (unset){ 
            left.initialize(cj);
            right.initialize(cj);
            bottom.initialize(cx);
            top.initialize(cx);
            left = a(ifirst, cj); 
            right = a(ilast, cj); 
            for (int i = xfirst; i < ifirst; i++) a(i, cj) =  left;
            for (int i = xlast; i > ilast; i--) a(i, cj) =  right;
            bottom = a(cx, jfirst);
            top = a(cx, jlast);
            for (int j = yfirst; j < jfirst; j++) a(cx, j) = bottom;
            for (int j = ylast; j > jlast; j--) a(cx, j) = top;
            unset = false;
        }else{
            for (int i = xfirst; i <= ifirst; i++) a(i, cj) =  left;
            for (int i = xlast; i >= ilast; i--) a(i, cj) =  right;
            for (int j = yfirst; j <= jfirst; j++) a(cx, j) = bottom;
            for (int j = ylast; j >= jlast; j--) a(cx, j) = top;
        }
    }
    friend std::ostream& operator<< (std::ostream &os, const FixedConst<T> &other){
        os << "FixedConst"; return os;
    }
};

struct Reflective{
    template<class A> void fix(A& a, const Interval<2>& cij){
        Interval<2> cxy = a.domain();
        int xfirst = cxy[0].first(), xlast = cxy[0].last();
        int yfirst = cxy[1].first(), ylast = cxy[1].last();
        int ifirst = cij[0].first(), ilast = cij[0].last();
        int jfirst = cij[1].first(), jlast = cij[1].last();
        Interval<1> cj = cij[1];
        Interval<1> cx = cxy[0];
        for (int i = xfirst; i < ifirst; i++) a(i, cj) = - a(2 * ifirst - i - 1, cj);
        for (int i = xlast; i > ilast; i--) a(i, cj) = - a(2 * ilast - i + 1, cj);
        for (int j = yfirst; j < jfirst; j++) a(cx, j) = - a(cx, 2 * jfirst - j - 1);
        for (int j = ylast; j > jlast; j--) a(cx, j) = - a(cx, 2 * jlast - j + 1);
    }
    friend std::ostream& operator<< (std::ostream &os, const Reflective &other){
        os << "Reflective"; return os;
    }
};

struct Mirror{
    template<class A> void fix(A& a, const Interval<2>& cij){
        Interval<2> cxy = a.domain();
        int xfirst = cxy[0].first(), xlast = cxy[0].last();
        int yfirst = cxy[1].first(), ylast = cxy[1].last();
        int ifirst = cij[0].first(), ilast = cij[0].last();
        int jfirst = cij[1].first(), jlast = cij[1].last();
        Interval<1> cj = cij[1];
        Interval<1> cx = cxy[0];
        for (int i = xfirst; i < ifirst; i++) a(i, cj) = a(2 * ifirst - i - 1, cj);
        for (int i = xlast; i > ilast; i--) a(i, cj) = a(2 * ilast - i + 1, cj);
        for (int j = yfirst; j < jfirst; j++) a(cx, j) = a(cx, 2 * jfirst - j - 1);
        for (int j = ylast; j > jlast; j--) a(cx, j) = a(cx, 2 * jlast - j + 1);
    }
    friend std::ostream& operator<< (std::ostream &os, const Mirror &other){
        os << "Mirror"; return os;
    }
};

struct Dependent{
    template<class A> void fix(A& a, const Interval<1>& ci){}
    template<class A> void fix(A& a, const Interval<2>& cij){}
    friend std::ostream& operator<< (std::ostream &os, const Dependent &other){
        os << "Dependent"; return os;
    }
};

template<class L = Dependent, class B = Dependent, class R = L, class T = B>
class Boundary{
public:
    L left; B bottom; R right; T top;
    template<class A> void fix(A& a, const Interval<2>& cij){
        // corner might got messed up
        Interval<2> cxy = a.domain();
        int ifirst = cij[0].first(), ilast = cij[0].last();
        int jfirst = cij[1].first(), jlast = cij[1].last();
        int xfirst = cxy[0].first(), xlast = cxy[0].last();
        int yfirst = cxy[1].first(), ylast = cxy[1].last();
        left.fix(a, Interval<2>(Interval<1>(ifirst, xlast), Interval<1>(yfirst, ylast)));
        right.fix(a, Interval<2>(Interval<1>(xfirst, ilast), Interval<1>(yfirst, ylast)));
        bottom.fix(a, Interval<2>(Interval<1>(xfirst, xlast), Interval<1>(jfirst, ylast)));
        top.fix(a, Interval<2>(Interval<1>(xfirst, xlast), Interval<1>(yfirst, jlast)));
    }
    friend std::ostream& operator<< (std::ostream &os, const Boundary<L,B,R,T> &other){
        os  << "-Left- " << other.left << ", "
            << "-Bottom- " << other.bottom << ", "
            << "-Right- " << other.right << ", "
            << "-Top- " << other.top; 
        return os;
    }
};
#endif
