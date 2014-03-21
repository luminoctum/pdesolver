#include <iostream>
#include "math.h"
#include "utils.h"
#include "functions.h"
#include "Pooma/Arrays.h"
#include "stdio.h"
#include "BurgerEquation.h"
#include "ShallowWater.h"
#include "MultiVar.h"
#include "Godunov.h"
#include "LaxFriedrich.h"
#include "EnoWeno.h"
#include "PDESolver.h"
using namespace std;

template<int DIM, class ET>
struct Method{
    //typedef ShallowWater<DIM, ET> Equation;
    typedef MultiVar<DIM, ET> Equation;
    typedef WENOScheme<3, DIM, Equation> Spatial;
};

int main(int argc, char **argv){
    Pooma::initialize(argc, argv);

    PDESolver<1, double, Method> pde(100, 3, 0.1);
    pde.solve(2);
    pde.observe();

    Pooma::finalize();
}
