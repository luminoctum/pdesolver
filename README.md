pdesolver
=========

A PDE Solver for conservation law written in C++

Current functionality:
----------------------

| Scheme | 1d1v (Implemented?/Tested?) | 1dnv | 2d1v  | 2dnv  | 3d1v | 3dnv |
| ------ | ---- | ---- | ----- | ---- | ---- | --- |
| LaxFriedrich | Yes/No | No/NO | No/No | No/No | No/No | No/NO |
| ENO/WENO | Yes/Yes | No/NO | No/No | No/No | No/No | No/NO |

Current work:
-------------

Adding support for 1dnv

1. Each variable is a distinct partial specialization of a tempalted class ... DONE
2. Each variable has its own (cell, wall) value, boundary condtion and property tags ... DONE
3. Based on loop unrooling techniques to loop over variables ...DONE
4. modify ENO/WENO reconstruct process to store the wall flux ... ON GOING
    - kind of finished the scalar version and decoupled vector version. runs ok but needs to check the result ... DONE
        - scalar version tested ... OK
        - decoupled vector version tested ... OK
    - finished modifying the Godunov Scheme. Did not run it, compiles ok ... ON GOING
        - Upwinding Godunov scheme tested ... OK
        - Entropy fix ... TODO
    - In order to use loop unrolling for variables, all the functions has to be static
    - this means that you want to generate a look up table for ENO/WENO coefficients, using mathematcia ... TODO
5. 2D multivariable model ... ON GOING
    - NCfile IO ... DONE
    - 2D Reconstruction using ENO Scheme ... DONE
    - update diaganostic variable ... ON GOING
