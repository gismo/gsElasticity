# gsElasticity
The gsElasticity module provides tools for structural mechanics analysis of 2D and 3D solids. This currently includes solvers for the following systems of equations:
- linear elasticity 
- nonlinear elasticity (with St. Venant-Kirchhoff or neo-Hookean material models)
- time-dependent problems (with ~~explicit~~ and implicit time integration) 
- mixed displacement-pressure formulation for (nearly) incompressible elasticity
- thermal expansion
Among additional features:
- incremental loading 
- isogeometric parametrization algorithms based on mesh deformation
- ~~eigenmode analysis~~

## Installing
To use the module, configure G+Smo with `GISMO_ELASTICITY=ON`. The module source code will be automatically downloaded from GitHub.com if an internet connection is available.

## Using
For detailed examples of using the module we refer to demo applications in `/extensions/gsElasticity/examples`. 
