#ifndef MICROPHYSICS
#define MICROPHYSICS

struct Species{
    double al, bl, mu, Rgas, Lv, mmr;
    std::string name;
    Species(){}
    Species(std::string _name) : name(_name){}
    // should use a more generalized user function
    template<class A>
    Array<2> svp_from_t(const A& temp){
        Array<2> result(temp.domain());
        result = pow(10., al + bl / temp);
        return result;
    }
    double inline smr_from_tp(double temp, double ptol){
        double svp = pow(10., al + bl / temp);
        double smr = svp / (ptol - svp);
        if ((smr < 0) || (smr > mmr)) smr = mmr;
        return smr;
    }

};

struct Water : public Species{
    Water(double _mmr = 1.E9) : 
        Species("H2O")
    {
        mu = 18.0E-3;
        Rgas = 461.67;
        al = 11.079;
        bl = -2261.1;
        Lv = 2.38E6;
        mmr = _mmr;
    }
};

struct Ammonia : public Species{
    Ammonia(double _mmr = 1.E9) : 
        Species("NH3")
    {
        mu = 17.E-3;
        Rgas = 488.82;
        al = 10.201;
        bl = -1248.;
        Lv = 1.4E6;
        mmr = _mmr;
    }
};

#endif
