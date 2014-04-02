#ifndef MICROPHYSICS
#define MICROPHYSICS
#include "Pooma/Arrays.h"
#include <vector>
struct Condensate{
    double va, vb, mu, Rgas, Lv, mmr;
    std::string name;
    Condensate(){}
    Condensate(std::string _name) : name(_name){}
    struct SvpFromT{
        double va, vb;
        void initialize(double _va, double _vb){ va = _va; vb = _vb; };
        double inline operator()(double temp) const { return pow(10., va + vb / temp); }
    };
    UserFunction<SvpFromT> svp_from_t;
};

struct Water : public Condensate{
    Water(double _mmr = 1.E9) : Condensate("H2O"){
        mu = 18.0E-3; Rgas = 461.67;
        va = 11.079; vb = -2261.1;
        Lv = 2.38E6;
        mmr = _mmr;
        svp_from_t.function().initialize(va, vb);
    }
};

struct Ammonia : public Condensate{
    Ammonia(double _mmr = 1.E9) : Condensate("NH3"){
        mu = 17.E-3; Rgas = 488.82;
        va = 10.201; vb = -1248.;
        Lv = 1.4E6;
        mmr = _mmr;
        svp_from_t.function().initialize(va, vb);
    }
};

template<int N>
struct CondensateList{
    Vector<N, double> va, vb, mu, Rgas, Lv, mmr;
    std::vector<Condensate*> speci;
    std::string name[N];
    CondensateList(){}
    CondensateList(Condensate _a, Condensate _b){
        va(0) = _a.va; va(1) = _b.va;
        vb(0) = _a.vb; vb(1) = _b.vb;
        mu(0) = _a.mu; mu(1) = _b.mu;
        Rgas(0) = _a.Rgas; Rgas(1) = _b.Rgas;
        Lv(0) = _a.Lv; Lv(1) = _b.Lv;
        mmr(0) = _a.mmr; mmr(1) = _b.mmr;
        speci = {&_a, &_b};
    }
};

#endif
