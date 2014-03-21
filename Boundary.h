#ifndef BOUNDARY
#define BOUNDARY
#include "Pooma/Arrays.h"
struct Periodic{
    template<class A> static void fix(A& a, const Interval<1>& ci){
        Interval<1> cx = a.domain();
        Interval<1> lower(ci.first(), cx.last() - ci.length());
        Interval<1> upper(ci.length() + cx.first(), ci.last());
        Interval<1> lowerfix(cx.first(), ci.first() - 1);
        Interval<1> upperfix(ci.last() + 1, cx.last());
        a(lowerfix) = a(upper);
        a(upperfix) = a(lower);
    }
};
struct LinearExtrap{
    template<class A> static void fix(A& a, const Interval<1>& ci){
        Interval<1> cx = a.domain();
        int ifirst = ci.first(), ilast = ci.last();
        int xfirst = cx.first(), xlast = cx.last();
        for (int i = ifirst - 1; i >= xfirst; i--) a(i) = a(i + 1) - a(ifirst + 1) + a(ifirst);
        for (int i = ilast + 1; i <= xlast; i++) a(i) = a(i - 1) + a(ilast) - a(ilast - 1);
    }
    template<class A> static void fix(A& a, const Interval<1>& cix, const Interval<1>& ciy){
        Interval<1> cxx = a.domain()[0];
        Interval<1> cxy = a.domain()[1];
        int ifirst = cix.first(), ilast = cix.last();
        int xfirst = cxx.first(), xlast = cxx.last();
        int jfirst = ciy.first(), jlast = ciy.last();
        int yfirst = cxy.first(), ylast = cxy.last();
        for (int i = ifirst - 1; i >= xfirst; i--) 
            a(i, ciy) = a(i + 1, ciy) - a(ifirst + 1, ciy) + a(ifirst, ciy);
        for (int i = ilast + 1; i <= xlast; i++)
            a(i, ciy) = a(i - 1, ciy) + a(ilast, ciy) - a(ilast - 1, ciy);
        for (int j = jfirst - 1; j >= yfirst; j--) 
            a(cxx, j) = a(cxx, j + 1) - a(cxx, jfirst + 1) + a(cxx, jfirst);
        for (int j = jlast + 1; j <= ylast; j++)
            a(cxx, j) = a(cxx, j - 1) + a(cxx, jlast) - a(cxx, jlast - 1);
    }
    template<class A> static void fix(A& a, const Interval<2>& cij){
        Interval<1> cix = cij[0];
        Interval<1> ciy = cij[1];
        fix(a, cix, ciy);
    }
};
struct ConstExtrap{
    template<class A> static void fix(A& a, const Interval<1>& ci){
        Interval<1> cx = a.domain();
        int ifirst = ci.first(), ilast = ci.last();
        int xfirst = cx.first(), xlast = cx.last();
        for (int i = ifirst - 1; i >= xfirst; i--) a(i) = a(ifirst);
        for (int i = ilast + 1; i <= xlast; i++) a(i) = a(ilast);
    }
    template<class A> static void fix(A& a, const Interval<2>& cij){
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
};
template<class ET, int ID = 0>
struct FixedConst{
    enum {id = ID};
    static bool unset;
    static ET lower, upper;
    static ET left, right, bottom, top;
    template<class A> static void fix(A& a, const Interval<1>& ci) {
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
    template<class A> static void fix(A& a, const Interval<2>& cij){
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
};
template<class ET, int ID>
bool FixedConst<ET, ID>::unset = true;
template<class ET, int ID>
ET FixedConst<ET, ID>::lower;
template<class ET, int ID>
ET FixedConst<ET, ID>::upper;
template<class ET, int ID>
ET FixedConst<ET, ID>::left;
template<class ET, int ID>
ET FixedConst<ET, ID>::right;
template<class ET, int ID>
ET FixedConst<ET, ID>::bottom;
template<class ET, int ID>
ET FixedConst<ET, ID>::top;

struct Reflective{
    template<class A> static void fix(A&a, const Interval<2>& cij){
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
};
struct Mirror{
    template<class A> static void fix(A&a, const Interval<2>& cij){
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
};

struct NONE{
    template<class A> static void fix(A&a, const Interval<1>& ci){}
    template<class A> static void fix(A&a, const Interval<2>& cij){}
};

template<class L = NONE, class B = NONE, class R = L, class T = B>
class Boundary{
public:
    template<class A> static void fix(A& a, const Interval<2>& cij){
        // corner might got messed up
        Interval<2> cxy = a.domain();
        int ifirst = cij[0].first(), ilast = cij[0].last();
        int jfirst = cij[1].first(), jlast = cij[1].last();
        int xfirst = cxy[0].first(), xlast = cxy[0].last();
        int yfirst = cxy[1].first(), ylast = cxy[1].last();
        L::fix(a, Interval<2>(Interval<1>(ifirst, xlast), Interval<1>(yfirst, ylast)));
        R::fix(a, Interval<2>(Interval<1>(xfirst, ilast), Interval<1>(yfirst, ylast)));
        B::fix(a, Interval<2>(Interval<1>(xfirst, xlast), Interval<1>(jfirst, ylast)));
        T::fix(a, Interval<2>(Interval<1>(xfirst, xlast), Interval<1>(yfirst, jlast)));
    }
};
#endif
