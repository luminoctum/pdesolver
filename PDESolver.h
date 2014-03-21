#ifndef PDESOLVER
#define PDESOLVER
#include "Stencil.h"
#include "Boundary.h"


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
    Interval_t ci, cx, si, sx;
    Stencil<Difference> diff;
    //FluxBuffer<DIM, ET> buffer;
    template<template<int> class T, int N> struct Loop{
        static void printName(){
            std::cout << T<N>::name << std::endl;
            Loop<T, N - 1>::printName();
        }
        template<class I>
        static void initialize(I cx, I sx){
            T<N>(cx, sx);
            Loop<T, N - 1>::initialize(cx, sx);
        }
        template<class I>
        static void fixBoundary(I ci){
            T<N>::bd::fix(T<N>::cell, ci);
            Loop<T, N - 1>::fixBoundary(ci);
        }
        template<class I>
        //note: static member function doesnot have a this parameter
        static void wallConstruct(I ci){
            construct(T<N>(), ci, typename T<N>::ctag());
            Loop<T, N - 1>::wallConstruct(ci);
        }
    };
    // note: if you missed the "STRUCT" you will have infinite loops when compiling
    template<template<int> class T> struct Loop<T, 0>{
        static void printName(){
            std::cout << T<0>::name << std::endl;
        }
        template<class I>
        static void initialize(I cx, I sx){
            T<0>(cx, sx);
        }
        template<class I>
        static void fixBoundary(I ci){
            T<0>::bd::fix(T<0>::cell, ci);
        }
        template<class I>
        static void wallConstruct(I ci){
            construct(T<0>(), ci, typename T<0>::ctag());
        }
    };
public:
    PDESolver(){}
    PDESolver(int _nx, int _gl, double _cfl) : 
        nx(_nx), gl(_gl), cfl(_cfl)
    {
        Interval_t _ci(0, _nx - 1);
        Interval_t _cx(- gl, _nx + _gl - 1);
        Interval_t _si(0, _nx);
        Interval_t _sx(-1, _nx);
        ci = _ci; si = _si; cx = _cx, sx = _sx;
        initialize();
    }
    void initialize(){
        Loop<Variable, 1>::printName();
        Loop<Variable, 1>::initialize(cx, sx);
        for (int i = ci.first(); i <= ci.last(); i++){
            Variable<0>::cell(i) = (i + 1) * (i + 1);
            Variable<1>::cell(i) = i + 2;
        }
        Loop<Variable, 1>::fixBoundary(ci);
        Loop<Variable, 1>::wallConstruct(ci);
    };
    void solve(double tend){
    }
    void observe(){ 
        std::cout << Variable<0>::cell  << std::endl;
        std::cout << Variable<0>::wallm << std::endl;
        std::cout << Variable<0>::wallp << std::endl;
        std::cout << Variable<1>::cell  << std::endl;
        std::cout << Variable<1>::wallm << std::endl;
        std::cout << Variable<1>::wallp << std::endl;
    }
};

#endif
