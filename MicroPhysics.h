#ifndef MICROPHYSICS
#define MICROPHYSICS
#include "Pooma/Arrays.h"
#include <vector>
struct Condensate{
    double va, vb, mu, Rgas, Lv, mmr;
    std::string name;
    Condensate(){}
    Condensate(std::string _name) : name(_name){}
    // should use a more genervaized user function
    template<class A>
    Array<2> svp_from_t(const A& temp){
        Array<2> result(temp.domain());
        result = pow(10., va + vb / temp);
        return result;
    }
    double inline smr_from_tp(double temp, double ptol){
        double svp = pow(10., va + vb / temp);
        double smr = svp / (ptol - svp);
        if ((smr < 0) || (smr > mmr)) smr = mmr;
        return smr;
    }

};

struct Water : public Condensate{
    Water(double _mmr = 1.E9) : Condensate("H2O")
    {
        mu = 18.0E-3; Rgas = 461.67;
        va = 11.079; vb = -2261.1;
        Lv = 2.38E6;
        mmr = _mmr;
    }
};

struct Ammonia : public Condensate{
    Ammonia(double _mmr = 1.E9) : Condensate("NH3")
    {
        mu = 17.E-3; Rgas = 488.82;
        va = 10.201; vb = -1248.;
        Lv = 1.4E6;
        mmr = _mmr;
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
    template<class A>
    Array<2, Vector<N, double> > svp_from_t(const A& temp){
        Array<2, Vector<N, double> > result(temp.domain());
        for (int i = 0; i < N; i++) result.comp(i) = speci[i]->svp_from_t(temp);
        return result;
    }
};

#endif
