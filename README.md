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

1. Each variable is a distinct partial specialization of a tempalted class ... done
2. Each variable has its own (cell, wall) value, boundary condtion and property tags ... done
3. Based on loop unrooling techniques to loop over variables ...done
4. modify ENO/WENO reconstruct process to store the wall flux ... on going
    - kind of finished the scalar version and decoupled vector version. runs ok but needs to check the result
    - In order to use loop unrolling for variables, all the functions has to be static
    - this means that you want to generate a look up table for ENO/WENO coefficients, using mathematcia ... TODO
