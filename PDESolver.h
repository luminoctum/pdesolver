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
    void initialize();
    double get_cfl() const { return cfl; }
    void solve(double tend){
    }
    void observe(){ 
        std::cout << Variable<0>::cell << std::endl;
    }
};

#endif
