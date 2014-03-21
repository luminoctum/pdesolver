#ifndef PDESOLVER
#define PDESOLVER
#include "Stencil.h"
#include "Boundary.h"

template<template<int> class T, int N> struct Loop{
    static void printName(){
        std::cout << T<N>::name << std::endl;
        Loop<T, N - 1>::printName();
    }
    template<class I>
    static void initialize(I cx, I si){
        T<N>(cx, si);
        Loop<T, N - 1>::initialize(cx, si);
    }
    template<class I>
    static void fixBoundary(I ci){
        T<N>::bd::fix(T<N>::cell, ci);
        Loop<T, N - 1>::fixBoundary(ci);
    }
};

// note: if you missed the "STRUCT" you will have infinite loops when compiling
template<template<int> class T> struct Loop<T, 0>{
    static void printName(){
        std::cout << T<0>::name << std::endl;
    }
    template<class I>
    static void initialize(I cx, I si){
        T<0>(cx, si);
    }
    template<class I>
    static void fixBoundary(I ci){
        T<0>::bd::fix(T<0>::cell, ci);
    }
};

template<int DIM, class ET, template<int,class> class OPT>
class PDESolver : 
    public OPT<DIM, ET>::Equation, 
    public OPT<DIM, ET>::Spatial
{
    typedef ET Element_t;
    typedef Array<DIM, ET> Array_t;
    typedef Interval<DIM> Interval_t;
    typedef typename OPT<DIM, ET>::Spatial Spatial;
    typedef typename OPT<DIM, ET>::Equation Equation;
private:
    double cfl, dt, dx;
    int nx, gl;
    Interval_t ci, cx, si;
    Stencil<Difference> diff;
    FluxBuffer<DIM, ET> buffer;
public:
    PDESolver(){}
    PDESolver(int _nx, int _gl, double _cfl) : 
        nx(_nx), gl(_gl), cfl(_cfl), buffer(_nx, _gl)
    {
        Interval_t _ci(0, _nx - 1);
        Interval_t _cx(- gl, _nx + _gl - 1);
        Interval_t _si(0, _nx);
        ci = _ci; si = _si; cx = _cx;
        initialize();
    }
    void initialize(){
        Loop<Variable, 1>::printName();
        Loop<Variable, 1>::initialize(cx, si);
        for (int i = ci.first(); i <= ci.last(); i++){
            Variable<0>::cell(i) = (i + 1) * (i + 1);
            Variable<1>::cell(i) = i + 2;
        }
        Loop<Variable, 1>::fixBoundary(ci);
    };
    void solve(double tend){
    }
    void observe(){ 
        std::cout << Variable<0>::cell << std::endl;
        std::cout << Variable<1>::cell << std::endl;
    }
};

#endif
