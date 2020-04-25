
# gsElasticity
gsElasticity is a **submodule of G+Smo** which has started as a collection of nonlinear elasticity solvers for 2D and 3D solids. Since then, its focus has shifted towards mesh deformation, and now gsElasticity contains isogeometric solvers for nonlinear elasticity and incompressible Navier-Stokes equations, various PDE-based mesh deformations algorithms and a partitioned fluid-structure interaction solver. Additionally, gsElasticity includes numerous application examples and the corresponding NURBS geometries.

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

## Installation
Since gsElasticity is a submodule of G+Smo, you should download G+Smo first:
```
git clone https://github.com/gismo/gismo.git
```
Then, configure G+Smo with `GISMO_ELASTICITY=ON`:
```
cd gismo
mkdir build
cd build
cmake .. -DGISMO_ELASTICITY=ON
```
This will trigger a download of gsElasticity from GitHub. Once gsElasticity is downloaded, you can compile G+Smo with gsElasticity:
```
make
```
Once complete, you can find the compiled library in `/path/to/gismo/build/lib` and examples in `/path/to/gismo/build/bin`.

#### Advanced installation
gsElasticity is an independent .git repository. By default, the procedure described above downloads the version of gsElasticity __that is linked to G+Smo__. Usually, it is the latest stable version. However, if you want to get other versions/branches of gsElastisity, or even contribute to it, you should access gsElasticity via its .git repository located in `/path/to/gismo/extenstions/gsElasticity`. For example, to get the latest version of gsElasticity, do
```
cd /path/to/gismo/extenstion/gsElasticity
git pull
```

Configuration with OpenMP
```
cd /path/to/gismo/build
cmake .. -DGISMO_ELASTICITY=ON -DGISMO_WITH_OPENMP=ON
```
Optimal configuration
```
cd /path/to/gismo/build
source /path/to/intel/compilers_and_libraries/linux/bin/compilervars.sh intel64
cmake .. -DGISMO_ELASTICITY=ON -DGISMO_WITH_OPENMP=ON -DEIGEN_USE_MKL_ALL=ON -DGISMO_WITH_PARDISO=ON -DPARDISO_USE_MKL=ON -DINTEL_ROOT=/path/to/intel -DCMAKE_C_COMPILER=icc -DCMAKE_CXX_COMPILER=icpc
```

Alternatively,
```
cd /path/to/gismo/build
cmake .. -DGISMO_ELASTICITY=ON -DGISMO_WITH_OPENMP=ON -DGISMO_WITH_PARDISO=ON -DPardiso_DIR=/path/to/pardiso
```



For faster compilation
```
make -j8
```

## Using
Predefined assemblers style. Into nonlinear or corresponding time integration. 
For detailed examples of using the module we refer to demo applications in `/extensions/gsElasticity/examples`. 
Comment to each example. Maybe with pictures and references

Link to the RG profile

### Flapping beam

Link to the FSI video

### Plate with hole

### Cooks



