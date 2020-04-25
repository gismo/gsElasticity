
# gsElasticity
gsElasticity is a **submodule** of G+Smo which has started as a collection of nonlinear elasticity solvers for 2D and 3D solids. Since then, its focus has shifted towards mesh deformation, and now gsElasticity includes isogeometric solvers for nonlinear elasticity and incompressible Navier-Stokes equations, various PDE-based mesh deformations algorithms and a partitioned fluid-structure interaction solver. 

------------------  Nonlinear elastic deformation of a 3D object  ----------------------------------------------------

<img src="../media/images/terrific.png" width="500">   

------------------  Fluid-structure interaction in 2D  ---------------------------------------------------------------

<img src="../media/images/FSI.png" width="450">

------------------  Mesh deformation for isogeometric domain parametrization -----------------------------------------

<img src="../media/images/meshDeform.png" width="550">

## Solvers
gsElasticity currently includes the following solvers for 2D and 3D multi-patch tensor-product NURBS geometries:
* Elasticity solvers 
  * linear elasticity
  * nonlinear elasticity with St.Venant-Kirchhoff and neo-Hookean material laws
  * implicit time integration with Newmark method
  * pure displacement and mixed displacement-pressure formulations
  * thermal expansion
* Incompressible Navier-Stokes solver
  * Stokes equation
  * stationary INSE
  * explicit and implicit time integration with the one-step theta-scheme
  * suitable for arbitrary Lagrangian-Eulerian (ALE) mappings
  * ~~SUPG stabilization~~
  * ~~turbulence model~~
* Fluid-structure interaction solver
  * partitioned approach
  * strong coupling
  * Aitken relaxation for convergence speed-up
* Bi-harmonic equation solver in mixed formulation
* Poisson's equation solver

## Installing
To use the module, configure G+Smo with `GISMO_ELASTICITY=ON`. The module source code will be automatically downloaded from GitHub.com if an internet connection is available.

mkdir, cd, cmake gsElasticty only, cmake with openmp and pardiso, release

update gsElasticity

part of gismo

## Using
Predefined assemblers style. Into nonlinear or corresponding time integration. 
For detailed examples of using the module we refer to demo applications in `/extensions/gsElasticity/examples`. 
Comment to each example. Maybe with pictures and references

Link to the RG profile

### Flapping beam

Link to the FSI video

### Plate with hole

### Cooks



