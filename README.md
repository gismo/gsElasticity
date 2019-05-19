# gsElasticity
The gsElasticity module provides tools for structural mechanics analysis of 2D and 3D solids. This currently includes solvers for the following systems of equations:
- linear elasticity ~~(with eigenmode analysis)~~
- nonlinear elasticity (with St. Venant-Kirchhoff or neo-Hookean material models)
- incremental loading and adaptive algoithms for problems with large deformations
- mixed displacement-pressure formulation for nearly incompressible elasticity
- thermal expansion
- isogeometric parametrization algorithms based on mesh deformation

## Installing
To use the module, configure G+Smo with `GISMO_ELASTICITY=ON`. The module source code will be automatically downloaded from GitHub.com if an internet connection is available.

## Using
For detailed examples of using the module we refer to demo applications in `/extensions/gsElasticity/examples`. 
