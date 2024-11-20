/** @file composite_elasticity_example.cpp

  @brief Expression assembler to solve the a linear elasticity problem

  Based on the poisson_example

  The test case is based on the proposed manufactured solution found here:
  https://web.mit.edu/16.20/homepage/4_ElasticityBVP/ElasticityBVP_files/module_4_with_solutions.pdf
  page: 82

  The solution must be degree 2 minimum to produce accurate results, use
  degree elevation and refinement to modify accuracy

  This file is part of the G+Smo library.

  This Source Code Form is subject to the terms of the Mozilla Public
  License, v. 2.0. If a copy of the MPL was not distributed with this
  file, You can obtain one at http://mozilla.org/MPL/2.0/.

  To run script:
  ./bin/composite_elasticity_example -f composite_example_mp.xml

*/

#include <gismo.h>
#include <gsElasticity/gsElasticityAssembler.h>
#include <gsElasticity/gsMaterialBase.h>
#include <gsElasticity/gsLinearMaterial.h>
#include <gsElasticity/gsCompositeMaterial.h>
#include <gsElasticity/gsCompositeMatrix.cpp>


using namespace gismo;

int main(int argc, char *argv[]) {
  //! [Parse command line]
  bool plot = false;
  index_t n_h_refinements = 0;
  index_t n_deg_elevations = 0;
  bool only_last = false;
  bool compute_error{false};
  index_t sample_rate{9};

  std::string file_name("composite_elasticity_example.xml");

  gsCmdLine cmd("Tutorial on solving a Linear Elasticity problem.");
  cmd.addInt("e", "degreeElevation",
             "Number of degree elevation steps to perform before solving (0: "
             "equalize degree in all directions)",
             n_deg_elevations);
  cmd.addInt("r", "uniformRefine", "Number of Uniform h-refinement loops",
             n_h_refinements);
  cmd.addString("f", "file", "Input XML file", file_name);
  cmd.addSwitch("only_last",
                "Solve solely for the only_last level of h-refinement",
                only_last);
  cmd.addSwitch("plot",
                "Create a ParaView visualization file with the solution", plot);
  cmd.addInt("q", "sample-rate", "Samples per spline in paraview export",
             sample_rate);
  cmd.addSwitch("compute-error",
                "Evaluate the error with respect to the analytical solution "
                "(evaluation with default options and default file required)",
                compute_error);

  // Error is computed for an analytical solution with the following
  // Lame constants
  real_t lame_lambda{80000.0}, lame_mu{80000.0};
  cmd.addReal("L", "firstLame", "First Lame constant, material parameter",
              lame_lambda);
  cmd.addReal("M", "secondLame", "Second Lame constant, material parameter",
              lame_mu);

  //! [Parse command line]
  try {
    cmd.getValues(argc, argv);
  } catch (int rv) {
    return rv;
  }

  //! [Read input file]
  gsFileData<> file_data(file_name);
  gsInfo << "Loaded file " << file_data.lastPath() << "\n";

  gsMultiPatch<> mp;
  file_data.getId(0, mp);  // id=0: Multipatch domain

  gsBoundaryConditions<> bc;
  file_data.getId(1, bc);  // id=2: boundary conditions
  bc.setGeoMap(mp);
  gsInfo << "Boundary conditions:\n" << bc << "\n";

  // gsOptionList assembler_options = gsAssembler<>::defaultOptions();
  // if (file_data.hasId(4)) {
  //   file_data.getId(4, assembler_options);  // id=4: assembler options
  // }
  //! [Read input file]

  //! [Refinement]
  gsMultiBasis<> dbasis(mp, true);
  const int solution_field_dimension = mp.geoDim();

  // Elevate and p-refine the basis
  dbasis.degreeElevate(n_deg_elevations);

  // h-refine each basis
  for (int r = 0; r < n_h_refinements; ++r) {
    dbasis.uniformRefine();
  }
  //! [Refinement]

  //! [Problem setup]
  gsExprAssembler<> A(1, 1);
  // if (file_data.hasId(4)) {
  //   A.setOptions(assembler_options);
  // }

  typedef gsExprAssembler<>::geometryMap GeometryMap;
  typedef gsExprAssembler<>::variable Variable;
  typedef gsExprAssembler<>::space Space;
  typedef gsExprAssembler<>::solution Solution;

  // Elements used for numerical integration
  A.setIntegrationElements(dbasis);
  gsExprEvaluator<> ev(A);

  // Set the geometry map
  GeometryMap G = A.getMap(mp);

  // Set the discretization space
  Space w = A.getSpace(dbasis, solution_field_dimension);

  // // Set the source term ?? do i need this?
  // auto source_function =
  //     A.getCoeff(source_function_expression, G);

  // Solution vector and solution variable
  gsMatrix<> solution_vector;
  Solution u_solution_expression = A.getSolution(w, solution_vector);
  //! [Problem setup]

  // ![User output]
  gsInfo << "Problem Overview:\t"
         << "\nGeometric Dim:\t" << solution_field_dimension << "\nPatches:\t"
         << mp.nPatches() << "\nMinimum degree:\t"
         << dbasis.minCwiseDegree() << "\nMaximum degree:\t"
         << dbasis.maxCwiseDegree() << "\n\n";
#ifdef _OPENMP
  gsInfo << "Available threads: " << omp_get_max_threads() << "\n";
#endif
  gsInfo << "Active assembly options:\n"
         << A.options() << std::endl;
  // ![User output]

  //! [Assembly]
  gsSparseSolver<>::CGDiagonal solver;

  real_t l2_err{}, h1_err{}, setup_time{}, ma_time{}, slv_time{}, err_time{};
  gsStopwatch timer;
  // dbasis.uniformRefine();

  // Apply dirichlet boundary conditions
  // w.setup(bc, dirichlet::l2Projection, 0);

  // Initialize the system
  A.initSystem();
  
  setup_time += timer.stop();
  // User output
  gsInfo << "Number of degrees of freedom:\t" << A.numDofs()
         << std::endl;
  gsInfo << "\nAssembling linear system ... " << std::flush;
  timer.restart();

  // Compute the system matrix and right-hand side
  auto phys_jacobian = ijac(w, G);

  // Material properties (0 deg ply)
  real_t E11 = 25000; // in MPa (N/mm^2)
  real_t E22 = E11/25;
  real_t E33 = E11/25;
  real_t G23 = 500;
  real_t G12 = G23/2.5;
  real_t G13 = G23/2.5;
  real_t nu12 = 0.25;
  real_t nu13 = 0.25;
  real_t nu23 = 0.25;

  // Use gsCompositeMaterial class!?!?!?!?1
  // 0 degree ply (como se a que gsCompositeMatrix accede????? getMaterialMatrix?)
  // Function in gsCompositeMatrix.cpp -- returns the stiffness matrix!
  gsMatrix<> Gmat = gsCompositeMatrix(E11,E22,E33,G12,G23,G13,nu12,nu23,nu13);

  // Ideally with sine cosine?
  // Otherwise, entering the values?
  // Rotation matrix for 90 deg ply (orthotropic)
  // gsMatrix<> T(6, 6);
  // T << 0, 1, 0, 0, 0, 0,
  //      1, 0, 0, 0, 0, 0,
  //      0,0,1,0,0,0,
  //      0,0,0,0,-1,0,
  //      0,0,0,1,0,0,
  //      0,0,0,0,0,-1;

  // gsMatrix<> T_inv(6, 6); // same matrix evaluated angle -90deg
  // T_inv << 0, 1, 0, 0, 0, 0,
  //         1, 0, 0, 0, 0, 0,
  //         0,0,1,0,0,0,
  //         0,0,0,0,1,0,
  //         0,0,0,-1,0,0,
  //         0,0,0,0,0,-1;

  // gsMatrix<> Gmat_90= T_inv*Gmat_0*T;

  gsConstantFunction<> fun_rot_0_z(0.0, 0.0, 0.0, mp.parDim()); // 3 is the dimension?
  gsConstantFunction<> fun_rot_90_z(0.0, 0.0, M_PI/2.0, mp.parDim()); // in radians?? or d

  Gmat.resize(Gmat.rows()*Gmat.cols(),1); //  needed for gs constant function!!!
  
  gsConstantFunction<> G_0(Gmat,mp.parDim()); //?
  // gsConstantFunction<> G_90(Gmat,mp.parDim()); //?

  // gsLinearMaterial<real_t> materialMat(youngsModulus,poissonsRatio);
  // Check the rotation matrix!!!!
  gsCompositeMaterial<real_t> materialMat_0(G_0,fun_rot_0_z); //??
  gsCompositeMaterial<real_t> materialMat_90(G_0,fun_rot_90_z); //??

  gsMaterialContainer<real_t> container(3); // preallocated

  container.set(0,materialMat_90);
  container.set(1,materialMat_0);
  container.set(2,materialMat_90);

  //gsMaterialContainer<real_t> container(materialMat_90,3); // preallocated
  //container.set(1,materialMat_0); //override it!

  // container.add(materialMat_90); // material 0 --> patch 0 (piece)
  // container.add(materialMat_0);
  // container.add(materialMat_90);

  // materialMat_0

  // creating assembler
  // auto g = 0;

  gsConstantFunction<> g(0,mp.parDim()); // it is a gsFunction
  gsElasticityAssembler<real_t> A(mp,dbasis,bc,g,container); // not with pointer
  
  A.assemble();

  //! [Solver]
  ma_time += timer.stop();
  gsInfo << "Done" << std::endl;
  gsInfo << "Solving Linear System ...... " << std::flush;
  timer.restart();

  // Solve linear system
  solver.compute(A.matrix());
  solution_vector = solver.solve(A.rhs());
  //! [Solver]

  slv_time += timer.stop();
  gsInfo << "Done" << std::endl;


  //! [Summarize User Output]
  gsInfo << "\n\nTotal time: " << setup_time + ma_time + slv_time + err_time
         << "\n";
  gsInfo << "     Setup: " << setup_time << "\n";
  gsInfo << "  Assembly: " << ma_time << "\n";
  gsInfo << "   Solving: " << slv_time << "\n";
  gsInfo << "     Norms: " << err_time << "\n";


  //! [Export visualization in ParaView]
  if (plot) {
    gsInfo << "Plotting in Paraview ... ";

    gsParaviewCollection collection("ParaviewOutput_composite/solution", &ev);
    collection.options().setSwitch("plotElements", true);
    collection.options().setInt("plotElements.resolution", sample_rate);
    collection.newTimeStep(&mp);
    collection.addField(u_solution_expression, "numerical solution");
    collection.saveTimeStep();
    collection.save();

    gsFileManager::open("ParaviewOutput_composite/solution.pvd");
    gsInfo << "Done" << std::endl;
  } else {
    gsInfo << "No output created, re-run with --plot to get a ParaView "
              "file containing the solution.\n";
  }
  //! [Export visualization in ParaView]

  return EXIT_SUCCESS;

}  // end main
