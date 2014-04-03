#ifndef FUNCTIONS
#define FUNCTIONS
#include "Pooma/Arrays.h" 

/* eigensystem
 * sovel the eigen value problem A * rp = dp * rp
 * and lp * A = lp * dp
 * where rp is the right eigenvector, lp is the left eigenvector
 * and lp * rp = I
 */
template<int DIM>
struct EigenSystem{
    template<class M>
    void solve(M m, M& lp, M& D, M& rp){}
};

// noets:
// The C++ standard says operator=, operator[], operator(), and operator-> must be non-static member functions
// default template paramter is not allowed in partial specialization
template<> 
struct EigenSystem<2>{
    template<class ET>
    static void solve(Tensor<2, ET> m, Tensor<2, ET>& lp, Tensor<2, ET>& D, Tensor<2, ET>& rp){
        ET a = m(0, 0), b = m(0, 1), c = m(1, 0), d = m(1, 1), e;
        e = sqrt(a * a + 4 * b * c - 2 * a * d + d * d);
        D(0, 0) = 0.5 * (a + d - e);
        D(0, 1) = 0.;
        D(1, 0) = 0.;
        D(1, 1) = 0.5 * (a + d + e);
        lp(0, 0) = - c / e;
        lp(0, 1) = (e + a - d) / (2. * e);
        lp(1, 0) = c / e;
        lp(1, 1) = (e - a + d) / (2. * e);
        rp(0, 0) = - (- a + d + e) / (2. * c);
        rp(0, 1) = - (- a + d - e) / (2. * c);
        rp(1, 0) = 1.;
        rp(1, 1) = 1.;
    }
};

template<class A>
void printarray(const A& a, Interval<2> cij){
    for (int i = cij[0].first(); i < cij[0].last(); i++){
        for (int j = cij[1].first(); j < cij[1].last(); j++)
            printf("%6.1f", a(i, j));
        printf("\n");
    }
    printf("\n");
}

#endif
