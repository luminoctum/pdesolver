#ifndef PDESOLVER
#define PDESOLVER
#include "Stencil.h"
#include "Boundary.h"


template<int DIM, class ET, template<int,class> class OPT>
class PDESolver : 
    public OPT<DIM, ET>::Equation, 
    public OPT<DIM, ET>::Spatial,
    public Godunov<DIM>
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
        template<class I>
        static void observe(I ci, I si){
            std::cout << T<N>::name << std::endl;
            std::cout << T<N>::cell(ci);
            std::cout << T<N>::wallm(si);
            std::cout << T<N>::wallp(si);
            std::cout << T<N>::flux(si) << std::endl;
            Loop<T, N - 1>::observe(ci, si);
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
        static void cell2wall(I ci){
            construct(T<N>(), ci, typename T<N>::ctag());
            Loop<T, N - 1>::cell2wall(ci);
        }
        template<class W, class I>
        static void wall2flux(const W& windm, const W& windp, I si){
            wallFlux(T<N>(), windm, windp, si, typename T<N>::ctag());
            Loop<T, N - 1>::wall2flux(windm, windp, si);
        }
    };
    // note: if you missed the "STRUCT" you will have infinite loops when compiling
    template<template<int> class T> struct Loop<T, 0>{
        template<class I>
        static void observe(I ci, I si){
            std::cout << T<0>::name << std::endl;
            std::cout << T<0>::cell(ci);
            std::cout << T<0>::wallm(si);
            std::cout << T<0>::wallp(si);
            std::cout << T<0>::flux(si) << std::endl;
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
        static void cell2wall(I ci){
            construct(T<0>(), ci, typename T<0>::ctag());
        }
        template<class W, class I>
        static void wall2flux(const W& windm, const W& windp, I si){
            wallFlux(T<0>(), windm, windp, si, typename T<0>::ctag());
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
        Interval_t _sx(-1, _nx + 1);
        ci = _ci; si = _si; cx = _cx, sx = _sx;
        initialize();
    }
    void initialize(){
        Loop<Variable, 2>::initialize(cx, sx);
        for (int i = ci.first(); i <= ci.last(); i++){
            Variable<0>::cell(i) = (i + 1) * (i + 1) - 10.;
            Variable<1>::cell(i) = (i + 1) * (i + 2);
            Variable<2>::cell(i)(0) = (i + 1) * (i + 1);
            Variable<2>::cell(i)(1) = (i + 1) * (i + 2);
        }
        Loop<Variable, 2>::fixBoundary(ci);
        Loop<Variable, 2>::cell2wall(ci);
        Loop<Variable, 2>::wall2flux(Variable<0>::wallm, Variable<0>::wallp, si);
    };
    void solve(double tend){
    }
    void observe(){ 
        Loop<Variable, 2>::observe(ci, si);
    }
};

#endif
