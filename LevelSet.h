#ifndef LEVELSET
#define LEVELSET
#include "netcdf.hh"
#include "utils.h"
#include "Stencil.h"
#include "Boundary.h"

class LevelSet{
private:
    double cfl, dt, dx, time;
    int nx, ny, gl;
    long current;
    double eps;
    Array<2> psi, psi0, psi_s, buffer;
    Array<2> psi_xm, psi_xp, psi_ym, psi_yp;
    Array<1> xa, ya;
    Interval<1> ci, cx, si;
    Interval<2> cij, cxy;
    Stencil<DifferenceX> diffx;
    Stencil<DifferenceY> diffy;
    Stencil<Difference2X> diff2x;
    Stencil<Difference2Y> diff2y;
    Stencil<Curvature> curv;
    LinearExtrap psibound, xbound, ybound;
    struct Positive{
        inline double operator()(double x) const { return _max_(x, 0);}
    };
    struct Negative{
        inline double operator()(double x) const { return _min_(x, 0);}
    };
    template<class EG1, class EG2>
    Array<2> minMod(const Array<2, double, EG1>& ua, const Array<2, double, EG2>& ub) const{
        Array<2> result(ua.domain());
        //Pooma::blockAndEvaluate();
        for (int i = ua.first(0); i <= ua.last(0); i++)
            for (int j = ua.first(1); j <= ua.last(1); j++){
                result(i, j) = _minmod_(ua.read(i, j), ub.read(i, j));
            }
        return result;
    }
    struct ReInitTag{
        inline double operator()(double x) const { 
            if ((x > 1.1) || (x < 0.9)) return 1;
            else return 0;
        }
    };
    UserFunction<Positive> pos_;
    UserFunction<Negative> neg_;
    UserFunction<ReInitTag> reinit_;
public:
    LevelSet(){}
    LevelSet(int _nx, int _gl, double _cfl) :
        nx(_nx), ny(_nx), gl(_gl), cfl(_cfl)
    {
        ci = Interval<1>(0, nx - 1);
        si = Interval<1>(0, nx);
        cx = Interval<1>(- gl, nx + gl - 1);
        cij = Interval<2>(ci, ci);
        cxy = Interval<2>(cx, cx);
        psi.initialize(cxy);
        psi0.initialize(cxy);
        psi_s.initialize(cxy);
        buffer.initialize(cij);
        psi_xm.initialize(cij);
        psi_xp.initialize(cij);
        psi_ym.initialize(cij);
        psi_yp.initialize(cij);
        xa.initialize(cx);
        ya.initialize(cx);
        initialize();
        psibound.fix(psi, cij);
        xbound.fix(xa, ci);
        ybound.fix(ya, ci);
        dx = xa(1) - xa(0);
        //dt = cfl * dx;
        dt = 0.005;
        current = 0;
        time = dt;
        eps = 0.025;
        //eps = 0.;
        DifferenceX::dx = dx;
        DifferenceY::dy = dx;
        Difference2X::dx = dx;
        Difference2Y::dy = dx;
        std::cout << dt << std::endl;
        std::cout << dx << std::endl;
    }
    void initialize(){
        //Stencils are not allowed to use here
        for (int i = 0; i < nx; i++) xa(i) = -2. + 4. * i / (nx - 1.);
        for (int j = 0; j < ny; j++) ya(j) = -2. + 4. * j / (nx - 1.);
        for (int i = 0; i < nx; i++) for (int j = 0; j < ny; j++){
            psi(i, j) = - pow(xa(i) * xa(i) + ya(j) * ya(j), 2) 
                + 2. * pow(xa(i) * xa(i) + ya(j) * ya(j), 2.5)
                - 4. * pow(xa(i) * ya(j), 2);
        }
    }
    void reInitialize(){
        psi0 = psi / sqrt(psi * psi + dx * dx);
        std::cout << "reInitializing" << std::endl;
        for (double t = 0; t < 0.4; t += dt){
            psi_xm = diffx(psi, cij + Loc<2>(-1, 0));
            psi_xp = diffx(psi, cij);
            psi_ym = diffy(psi, cij + Loc<2>(0, -1));
            psi_yp = diffy(psi, cij);
            psi(cij) = psi(cij) 
                - dt * pos_(psi0) * (sqrt(
                    pow(pos_(psi_xm), 2) + pow(neg_(psi_xp), 2) + pow(pos_(psi_ym), 2) + pow(neg_(psi_yp), 2)
                    ) - 1.)
                - dt * neg_(psi0) * (sqrt(
                    pow(pos_(psi_xp), 2) + pow(neg_(psi_xm), 2) + pow(pos_(psi_yp), 2) + pow(neg_(psi_ym), 2)
                    ) - 1.);
        }
    }
    int checkReInit(){
        //reInitialize();
        buffer = reinit_(sqrt(pow(diffx(psi, cij), 2) + pow(diffy(psi, cij), 2)));
        Interval<1> center(0.75 * ci.first() + 0.25 * ci.last(), 0.25 * ci.first() + 0.75 * ci.last());
        if (sum(buffer(center, center)) > 200){
            std::cout << sum(buffer(center, center)) << std::endl;
            reInitialize();
        }
    }
    void solve(double tend){
        for (; time < tend + 1.E-6; time += dt){
            std::cout << time << std::endl;
            psi_xm = diffx(psi, cij + Loc<2>(-1, 0)) + 0.5 * dx * minMod(diff2x(psi, cij + Loc<2>(-1, 0)), diff2x(psi, cij));
            psi_xp = diffx(psi, cij) - 0.5 * dx * minMod(diff2x(psi, cij), diff2x(psi, cij + Loc<2>(1, 0)));
            psi_ym = diffy(psi, cij + Loc<2>(0, -1)) + 0.5 * dx * minMod(diff2y(psi, cij + Loc<2>(0, -1)), diff2y(psi, cij));
            psi_yp = diffy(psi, cij) - 0.5 * dx * minMod(diff2y(psi, cij), diff2y(psi, cij + Loc<2>(0, 1)));
            psi_s(cij) = psi(cij) - dt * (1. - eps * curv(psi, cij)) * sqrt(
                    pow(pos_(psi_xm), 2) + pow(neg_(psi_xp), 2) + pow(pos_(psi_ym), 2) + pow(neg_(psi_yp), 2)
                    );
            psibound.fix(psi_s, cij);

            psi_xm = diffx(psi_s, cij + Loc<2>(-1, 0)) + 0.5 * dx * minMod(diff2x(psi_s, cij + Loc<2>(-1, 0)), diff2x(psi_s, cij));
            psi_xp = diffx(psi_s, cij) - 0.5 * dx * minMod(diff2x(psi_s, cij), diff2x(psi_s, cij + Loc<2>(1, 0)));
            psi_ym = diffy(psi_s, cij + Loc<2>(0, -1)) + 0.5 * dx * minMod(diff2y(psi_s, cij + Loc<2>(0, -1)), diff2y(psi_s, cij));
            psi_yp = diffy(psi_s, cij) - 0.5 * dx * minMod(diff2y(psi_s, cij), diff2y(psi_s, cij + Loc<2>(0, 1)));
            psi(cij) = 0.5 * psi(cij) + 0.5 * (psi_s(cij) - dt * (1. - eps * curv(psi_s, cij)) * sqrt(
                    pow(pos_(psi_xm), 2) + pow(neg_(psi_xp), 2) + pow(pos_(psi_ym), 2) + pow(neg_(psi_yp), 2)
                    ));
            psibound.fix(psi, cij);
            //reInitialize();
            checkReInit();
            observe();
        }
    }
    void observe(){ 
        current++;
        buffer = psi(cij);
        NcFile dataFile("levelset.nc", NcFile::Write);
        dataFile.get_var("psi")->put_rec(&buffer(0, 0), current);
        dataFile.get_var("time")->put_rec(&time, current);
    }
};
#endif
