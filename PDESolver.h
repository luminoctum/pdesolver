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
        template<class A, class D, class I, int SZ, class VET>
        static void observe(A& ncvar, long current, const D& nc, I cij, Vector<SZ, VET>, Int2Type<0>){
            for (int i = 0; i < T<N>::Element_t::d1; i++){
                buffer_out = T<N>::cell(cij).comp(i);
                nc.get_var(T<N>::name[i].c_str())->put_rec(&buffer_out(0, 0), current);
            }
            Loop<T, N - 1>::observe(ncvar, current, nc, cij, typename T<N - 1>::Element_t(), typename T<N - 1>::mtag());
        }
        // output scalar into ncfile
        template<class A, class D, class I, class SET>
        static void observe(A& ncvar, long current, const D& nc, I cij, SET, Int2Type<0>){
            buffer_out = T<N>::cell(cij);
            nc.get_var(T<N>::name.c_str())->put_rec(&buffer_out(0, 0), current);
            Loop<T, N - 1>::observe(ncvar, current, nc, cij, typename T<N - 1>::Element_t(), typename T<N - 1>::mtag());
        }
        template<class A, class D, class I, class SET>
        static void observe(A& ncvar, long current, const D& nc, I cij, SET, Int2Type<1>){
            buffer_out = T<N>::cell(cij) / ncvar["mass"];
            nc.get_var(T<N>::name.c_str())->put_rec(&buffer_out(0, 0), current);
            Loop<T, N - 1>::observe(ncvar, current, nc, cij, typename T<N - 1>::Element_t(), typename T<N - 1>::mtag());
        }
        template<class I>
        static void observe(I cij){
            std::cout << T<N>::name << std::endl;
            std::cout << T<N>::cell(cij) << std::endl;
            std::cout << T<N>::wallX(Interval<2>(Interval<1>(cij[0].first(), cij[0].last() + 1),cij[1])) << std::endl;
            std::cout << T<N>::wallY(Interval<2>(cij[0],Interval<1>(cij[1].first(), cij[1].last() + 1))) << std::endl;
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
        static void wall2flux(const W& windX, const W& windY, I sIj, I siJ){
            wallFlux(T<N>(), windX, windY, sIj, siJ, typename T<N>::ctag());
            Loop<T, N - 1>::wall2flux(windX, windY, sIj, siJ);
        }
        template<class A, class I>
        static void advect(A& ncvar, I cij, Int2Type<1>){
            T<N>::dcelldt(cij) += fdx(T<N>::fluxX, cij) + fdy(T<N>::fluxY, cij);
            Loop<T, N - 1>::advect(ncvar, cij, typename T<N - 1>::mtag());
        }
        template<class A, class I>
        static void advect(A& ncvar, I cij, Int2Type<0>){
            T<N>::dcelldt(cij) += (fdx(T<N>::fluxX, cij) + fdy(T<N>::fluxY, cij)) / ncvar["mass"];
            //Loop<T, N - 1>::advect(ncvar, cij, typename T<N - 1>::mtag());
        }
        template<class I>
        static void updateTendency(I cij, double dt){
            T<N>::cell(cij) += T<N>::dcelldt(cij) * dt;
            T<N>::dcelldt = 0;
            Loop<T, N - 1>::updateTendency(cij, dt);
        }
    };
    // note: if you missed the "STRUCT" you will have infinite loops when compiling
    template<template<int> class T> struct Loop<T, 0>{
        // output vector into ncfile
        template<class D, class I, int SZ, class VET>
        static void observe(long current, const D& nc, I cij, Vector<SZ, VET>, Int2Type<0>){
            for (int i = 0; i < T<0>::Element_t::d1; i++){
                buffer_out = T<0>::cell(cij).comp(i);
                nc.get_var(T<0>::name[i].c_str())->put_rec(&buffer_out(0, 0), current);
            }
        }
        // output scalar into ncfile
        /*
        template<class A, class D, class I, class SET>
        static void observe(A& long current, const D& nc, I cij, SET, Int2Type<0>){
            buffer_out = T<0>::cell(cij);
            nc.get_var(T<0>::name.c_str())->put_rec(&buffer_out(0, 0), current);
        }*/
        template<class A, class D, class I, class SET>
        static void observe(A& ncvar, long current, const D& nc, I cij, SET, Int2Type<1>){
            buffer_out = T<0>::cell(cij) / ncvar["mass"];
            nc.get_var(T<0>::name.c_str())->put_rec(&buffer_out(0, 0), current);
        }
        template<class I>
        static void observe(I cij){
            std::cout << T<0>::name << std::endl;
            std::cout << T<0>::cell(cij) << std::endl;
            std::cout << T<0>::wallX(Interval<2>(Interval<1>(cij[0].first(), cij[0].last() + 1),cij[1])) << std::endl;
            std::cout << T<0>::wallY(Interval<2>(cij[0],Interval<1>(cij[1].first(), cij[1].last() + 1))) << std::endl;
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
        static void wall2flux(const W& windX, const W& windY, I sIj, I siJ){
            wallFlux(T<0>(), windX, windY, sIj, siJ, typename T<0>::ctag());
        }
        template<class A, class I>
        static void advect(A&, I cij, Int2Type<1>){
            T<0>::dcelldt(cij) += fdx(T<0>::fluxX, cij) + fdy(T<0>::fluxY, cij);
        }
        /*
        template<class A, class I>
        static void advect(A& ncvar, I cij, Int2Type<0>){
            T<0>::dcelldt(cij) += (fdx(T<0>::fluxX, cij) + fdy(T<0>::fluxY, cij)) / ncvar["mass"];
        }*/
        template<class I>
        static void updateTendency(I cij, double dt){
            T<0>::cell(cij) += T<0>::dcelldt(cij) * dt;
            T<0>::dcelldt = 0;
        }
    };

    double model_time;
    Interval<1> ci, cj, cx, cy, si, sj, sx, sy;
    Interval_t cij, cxy, sIj, siJ, sXy, sxY, ijpeek;
    static Stencil<ForwardDifferenceX> fdx;
    static Stencil<ForwardDifferenceY> fdy;
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
        CenterDifferenceX::dx = dx;
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
        //std::cout << this->mixr.cell.comp(0) << std::endl;
        updateDiagnostics(ncvar, cij);
        Loop<Variable, 8>::fixBoundary(cij);
    };
    void solve(double tend){
        Loop<Variable, 4>::cell2wall(cij);
        Loop<Variable, 3>::wall2flux(Variable<0>::wallX, Variable<4>::wallY, sIj, siJ);
        Loop<Variable, 2>::advect(ncvar, cij, typename Variable<2>::mtag());
        bodyForce(ncvar, cij);
        Loop<Variable, 3>::updateTendency(cij, dt);
        Loop<Variable, 8>::fixBoundary(cij);
        updateDiagnostics(ncvar, cij);
        //Loop<Variable, 3>::observe(ijpeek);
        observe();
    }
    void observe(){ 
        this->current++;
        NcFile dataFile(this->filename.c_str(), NcFile::Write);
        Loop<Variable, 8>::observe(ncvar, this->current, dataFile, cij, typename Variable<8>::Element_t(), typename Variable<8>::mtag());
        dataFile.get_var("time")->put_rec(&model_time, this->current);
    }
};

template<int DIM, class ET, template<int,class> class OPT>
Array<DIM, ET> PDESolver<DIM, ET, OPT>::buffer_out;
template<int DIM, class ET, template<int,class> class OPT>
Stencil<ForwardDifferenceX> PDESolver<DIM, ET, OPT>::fdx;
template<int DIM, class ET, template<int,class> class OPT>
Stencil<ForwardDifferenceY> PDESolver<DIM, ET, OPT>::fdy;

#endif
