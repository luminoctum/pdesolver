#ifndef PDESOLVER
#define PDESOLVER
#include "NcFileIO.h"
#include "Stencil.h"
#include "Boundary.h"


template<int DIM, class ET, template<int,class> class OPT>
class PDESolver : 
    public OPT<DIM, ET>::Equation, 
    public OPT<DIM, ET>::Spatial,
    public Godunov<DIM>,
    public NcFileIO
{
    typedef ET Element_t;
    typedef Array<DIM, ET> Array_t;
    typedef Interval<DIM> Interval_t;
    typedef typename OPT<DIM, ET>::Spatial Spatial;
    typedef typename OPT<DIM, ET>::Equation Equation;
private:
    template<template<int> class T, int N> struct Loop{
        // output vector into ncfile
        template<class D, class I, int SZ, class VET>
        static void observe(long current, const D& nc, I cij, Vector<SZ, VET>){
            for (int i = 0; i < T<N>::Element_t::d1; i++){
                buffer_out = T<N>::cell(cij).comp(i);
                nc.get_var(T<N>::name[i].c_str())->put_rec(&buffer_out(0, 0), current);
            }
            Loop<T, N - 1>::observe(current, nc, cij, typename T<N - 1>::Element_t());
        }
        // output scalar into ncfile
        template<class D, class I, class SET>
        static void observe(long current, const D& nc, I cij, SET){
            buffer_out = T<N>::cell(cij);
            nc.get_var(T<N>::name.c_str())->put_rec(&buffer_out(0, 0), current);
            Loop<T, N - 1>::observe(current, nc, cij, typename T<N - 1>::Element_t());
        }
        template<class I>
        static void observe(I cij){
            std::cout << T<N>::name << std::endl;
            std::cout << T<N>::cell(cij) << std::endl;
            std::cout << T<N>::wallXm(cij) << std::endl;
            std::cout << T<N>::wallXp(cij) << std::endl;
            std::cout << T<N>::wallYm(cij) << std::endl;
            std::cout << T<N>::wallYp(cij) << std::endl;
            Loop<T, N - 1>::observe(cij);
        }
        template<class I>
        static void initialize(I cxy, I sXy, I sxY){
            T<N>(cxy, sXy, sxY);
            Loop<T, N - 1>::initialize(cxy, sXy, sxY);
        }
        template<class I>
        static void fixBoundary(I cij){
            T<N>::bd::fix(T<N>::cell, cij);
            Loop<T, N - 1>::fixBoundary(cij);
        }
        template<class I>
        //note: static member function doesnot have a this parameter
        static void cell2wall(I cij){
            construct(T<N>(), cij, typename T<N>::ctag());
            Loop<T, N - 1>::cell2wall(cij);
        }
        template<class W, class I>
        static void wall2fluxX(const W& windm, const W& windp, I sIj){
            //FIXME: wallFluxX
            wallFlux(T<N>(), windm, windp, sIj, typename T<N>::ctag());
            Loop<T, N - 1>::wall2fluxX(windm, windp, sIj);
        }
        template<class W, class I>
        static void wall2fluxY(const W& windm, const W& windp, I siJ){
            //FIXME: wallFluxY
            wallFlux(T<N>(), windm, windp, siJ, typename T<N>::ctag());
            Loop<T, N - 1>::wall2fluxY(windm, windp, siJ);
        }
    };
    // note: if you missed the "STRUCT" you will have infinite loops when compiling
    template<template<int> class T> struct Loop<T, 0>{
        // output vector into ncfile
        template<class D, class I, int SZ, class VET>
        static void observe(long current, const D& nc, I cij, Vector<SZ, VET>){
            for (int i = 0; i < T<0>::Element_t::d1; i++){
                buffer_out = T<0>::cell(cij).comp(i);
                nc.get_var(T<0>::name[i].c_str())->put_rec(&buffer_out(0, 0), current);
            }
        }
        // output scalar into ncfile
        template<class D, class I, class SET>
        static void observe(long current, const D& nc, I cij, SET){
            buffer_out = T<0>::cell(cij);
            nc.get_var(T<0>::name.c_str())->put_rec(&buffer_out(0, 0), current);
        }
        template<class I>
        static void observe(I cij){
            std::cout << T<0>::name << std::endl;
            std::cout << T<0>::cell(cij) << std::endl;
            std::cout << T<0>::wallXm(cij) << std::endl;
            std::cout << T<0>::wallXp(cij) << std::endl;
            std::cout << T<0>::wallYm(cij) << std::endl;
            std::cout << T<0>::wallYp(cij) << std::endl;
        }
        template<class I>
        static void initialize(I cxy, I sXy, I sxY){
            T<0>(cxy, sXy, sxY);
        }
        template<class I>
        static void fixBoundary(I cij){
            T<0>::bd::fix(T<0>::cell, cij);
        }
        template<class I>
        static void cell2wall(I cij){
            construct(T<0>(), cij, typename T<0>::ctag());
        }
        template<class W, class I>
        static void wall2fluxX(const W& windm, const W& windp, I sIj){
            wallFlux(T<0>(), windm, windp, sIj, typename T<0>::ctag());
        }
        template<class W, class I>
        static void wall2fluxY(const W& windm, const W& windp, I siJ){
            wallFlux(T<0>(), windm, windp, siJ, typename T<0>::ctag());
        }
    };

    double model_time;
    Interval<1> ci, cj, cx, cy, si, sj, sx, sy;
    Interval_t cij, cxy, sIj, siJ, sXy, sxY, ijpeek;
    Stencil<ForwardDifferenceX> diffx;
    Stencil<ForwardDifferenceY> diffy;
    static Array_t buffer_out;

public:
    PDESolver(){}
    PDESolver(std::string _str) : NcFileIO(_str)
    {
        ci = Interval<1>(0, nx - 1);
        cj = Interval<1>(0, ny - 1);
        cx = Interval<1>(- gl, nx + gl - 1);
        cy = Interval<1>(- gl, ny + gl - 1);
        si = Interval<1>(0, nx);
        sj = Interval<1>(0, ny);
        sx = Interval<1>(- 1, nx + 1);
        sy = Interval<1>(- 1, ny + 1);
        cij = Interval_t(ci, cj);
        sIj = Interval_t(si, cj);
        siJ = Interval_t(ci, sj);
        cxy = Interval_t(cx, cy);
        sXy = Interval_t(sx, cy);
        sxY = Interval_t(cx, sy);
        initialize();
        buffer_out.initialize(cij);
        ForwardDifferenceX::dx = dx;
        ForwardDifferenceY::dy = dy;
        std::cout << dx << std::endl;
        std::cout << dy << std::endl;
        this->finty.dy = dy;
        this->binty.dy = dy;
    }
    void initialize(){
        Loop<Variable, 8>::initialize(cxy, sXy, sxY);
        this->uwind.cell(cij) = ncvar["uwind"];
        this->vwind.cell(cij) = ncvar["vwind"];
        this->theta.cell(cij) = ncvar["tc"];
        this->mixr.cell(cij).comp(0) = ncvar["xH2O"];
        this->mixr.cell(cij).comp(1) = ncvar["xNH3"];
        this->phi.cell(cij) = ncvar["phi"];
        this->wwind.cell(cij) = ncvar["wwind"];
        this->temp.cell(cij).comp(0) = ncvar["temp"];
        this->temp.cell(cij).comp(1) = ncvar["tempv"];
        this->moist.cell(cij).comp(0) = ncvar["svpH2O"];
        this->moist.cell(cij).comp(1) = ncvar["svpNH3"];
        this->moist.cell(cij).comp(2) = ncvar["hH2O"];
        this->moist.cell(cij).comp(3) = ncvar["hNH3"];
        this->moist.cell(cij).comp(4) = ncvar["qH2O"];
        this->moist.cell(cij).comp(5) = ncvar["qNH3"];
        this->pdry.cell(cij) = ncvar["pdry"];

        ijpeek = Interval_t(Interval<1>(5, 10), Interval<1>(5, 10));
        Loop<Variable, 8>::fixBoundary(cij);
        //std::cout << this->mixr.cell.comp(0) << std::endl;
        Loop<Variable, 3>::cell2wall(cij);
        //Loop<Variable, 2>::wall2flux(Variable<0>::wallm, Variable<0>::wallp, sij);
        //Loop<Variable, 3>::observe(ijpeek);
    };
    void solve(double tend){
        updateDiagnostics(ncvar, cij);
    }
    void observe(){ 
        this->current++;
        NcFile dataFile(this->filename.c_str(), NcFile::Write);
        Loop<Variable, 8>::observe(this->current, dataFile, cij, typename Variable<8>::Element_t());
        dataFile.get_var("time")->put_rec(&model_time, this->current);
    }
};

template<int DIM, class ET, template<int,class> class OPT>
Array<DIM, ET> PDESolver<DIM, ET, OPT>::buffer_out;

#endif
