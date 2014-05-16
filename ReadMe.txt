More detailed description can be found at
http://readerunner.wordpress.com/2014/04/30/multigrid-methods-with-repa/

The main test example is

poissonTest

which requires user input of a grid stack size (probably less than 12 to avoid very large arrays)
and a filename for a bmp output.

This uses the module

MultiGrid_Poisson

which exports a full multigrid (fullMG) function to solve poissons equation.
MultiGrid_Poisson imports two finite difference calculations (approxOp and residualOp) from the module

PosissonOps

Solving other linear PDEs is possible with appropriate alternatives to the PoissonOps

Only square 2D grids are supported with equal spacing in each dimension.



