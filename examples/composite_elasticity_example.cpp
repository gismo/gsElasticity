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
  ./bin/composite_elasticity_example -r 4

*/

#include <gismo.h>
#include <gsElasticity/gsElasticityAssembler.h>
#include <gsElasticity/gsWriteParaviewMultiPhysics.h>
#include <gsElasticity/gsGeoUtils.h>

#include <gsElasticity/gsMaterialBase.h>
#include <gsElasticity/gsLinearMaterial.h>
// #include <gsElasticity/gsCompositeMaterial.h>
// #include <gsElasticity/gsCompositeMatrix.cpp>
#include <gsElasticity/gsCompositeMat.h>
#include <gsElasticity/gsMaterialEval.h>
#include <gsElasticity/gsVisitorElUtils.h>





using namespace gismo;

int main(int argc, char *argv[]) {
  //! [Parse command line]
  bool plot = false;
  index_t n_h_refinements = 0;
  index_t n_deg_elevations = 0;
  bool only_last = false;
  bool compute_error{false};
  index_t sample_rate{9};
  index_t numPlotPoints = 10000;


  std::string file_name("composite_example_mp.xml");

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
  cmd.addInt("p","points","Number of points to plot to Paraview",numPlotPoints);


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

  gsMultiPatch<> mp, mp_def;
  file_data.getId(0, mp);  // id=0: Multipatch domain

  gsDebugVar(mp);

  gsBoundaryConditions<> bc;
  //file_data.getId(1, bc);  // id=2: boundary conditions
  
  // Dirichlet BCs
  // Dirichlet x
  bc.addCondition(0,boundary::east,condition_type::dirichlet,0,0);
  bc.addCondition(0,boundary::west,condition_type::dirichlet,0,0);
  bc.addCondition(1,boundary::east,condition_type::dirichlet,0,0);
  bc.addCondition(1,boundary::west,condition_type::dirichlet,0,0);
  bc.addCondition(2,boundary::east,condition_type::dirichlet,0,0);
  bc.addCondition(2,boundary::west,condition_type::dirichlet,0,0);
  // Dirichlet y
  bc.addCondition(0,boundary::north,condition_type::dirichlet,0,1);
  bc.addCondition(0,boundary::south,condition_type::dirichlet,0,1);
  bc.addCondition(1,boundary::north,condition_type::dirichlet,0,1);
  bc.addCondition(1,boundary::south,condition_type::dirichlet,0,1);
  bc.addCondition(2,boundary::north,condition_type::dirichlet,0,1);
  bc.addCondition(2,boundary::south,condition_type::dirichlet,0,1);
  // Dirichlet z
  bc.addCondition(0,boundary::north,condition_type::dirichlet,0,2);
  bc.addCondition(0,boundary::south,condition_type::dirichlet,0,2);
  bc.addCondition(0,boundary::east,condition_type::dirichlet,0,2);
  bc.addCondition(0,boundary::west,condition_type::dirichlet,0,2);
  bc.addCondition(1,boundary::north,condition_type::dirichlet,0,2);
  bc.addCondition(1,boundary::south,condition_type::dirichlet,0,2);
  bc.addCondition(1,boundary::east,condition_type::dirichlet,0,2);
  bc.addCondition(1,boundary::west,condition_type::dirichlet,0,2);
  bc.addCondition(2,boundary::north,condition_type::dirichlet,0,2);
  bc.addCondition(2,boundary::south,condition_type::dirichlet,0,2);
  bc.addCondition(2,boundary::east,condition_type::dirichlet,0,2);
  bc.addCondition(2,boundary::west,condition_type::dirichlet,0,2);

  // Apply Neumann BC
  gsFunctionExpr<>func("0","0","1",mp.parDim()); // tiene sentido?????
  bc.addCondition(2,boundary::back,condition_type::neumann,&func);

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

  // for (int r = 0; r < n_h_refinements; ++r) {
  //   dbasis.uniformRefineComponent(0);
  //   dbasis.uniformRefineComponent(1);
  // }

  // for (int r = 0; r < 2; ++r) {
  //   dbasis.uniformRefineComponent(2);
  // }

  // h-refine each basis
  for (int r = 0; r < n_h_refinements; ++r) {
    dbasis.uniformRefine();
  }
  //! [Refinement]

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


  // Donde las doy como input???????????????? Importante! :(
  gsConstantFunction<> rot_0(0.0, 0.0, 0.0, mp.parDim()); // 3 is the dimension?
  gsConstantFunction<> rot_90(0.0, 0.0, M_PI/2.0, mp.parDim()); // in radians?? or d

  // ========== Material class ==========
  gsMaterialBase<real_t> * ply_mat_0;
  gsMaterialBase<real_t> * ply_mat_90;

  if      (mp.geoDim()==2)
  {
    ply_mat_0  = new gsCompositeMat<2,real_t>(E11,E22,G12,nu12,rot_0,mp);
    ply_mat_90 = new gsCompositeMat<2,real_t>(E11,E22,G12,nu12,rot_90,mp);
  }
  else if (mp.geoDim()==3)
  {
    ply_mat_0  = new gsCompositeMat<3,real_t>(E11,E22,E33,G12,G23,G13,nu12,nu23,nu13,rot_0,mp);
    ply_mat_90 = new gsCompositeMat<3,real_t>(E11,E22,E33,G12,G23,G13,nu12,nu23,nu13,rot_90,mp);
  }     
  else
      GISMO_ERROR("Dimension not supported");

  // if      (mp.geoDim()==2)
  // {
  //   mp_def = mp;
  //   ply_mat_0  = new gsLinearMaterial<2,real_t>(E11,0.25,mp,mp_def);
  //   ply_mat_90 = new gsLinearMaterial<2,real_t>(E22,0.25,mp,mp_def);
  // }
  // else if (mp.geoDim()==3)
  // {
  //   mp_def = mp;
  //   ply_mat_0  = new gsLinearMaterial<3,real_t>(E11,0.25,mp,mp_def);
  //   ply_mat_90 = new gsLinearMaterial<3,real_t>(E22,0.25,mp,mp_def);
  // }     
  // else
  //     GISMO_ERROR("Dimension not supported");

  // ======================================

  // Option 1
  // gsMaterialContainer<real_t> container(3); // preallocated
  // container.add(&ply_mat_90); // material 0 --> patch 0 (piece)
  // container.add(&ply_mat_0);
  // container.add(&ply_mat_90);

  // Option 2
  //  gsDebugVar(mp.nPatches());
  gsMaterialContainer<real_t> container(ply_mat_90,mp.nPatches()); // 3 containers! preallocated (all patches have the properties of layered_mat_90)
  container.set(1,ply_mat_0); //overload it to assign the material props of the ply 0 degrees
  
  gsConstantFunction<> g(0,0,0,mp.parDim()); // it is a gsFunction (body load?) comprobar??????? no entiendo!!

  // Call the elasticity assembler!
  gsElasticityAssembler<real_t> A(mp,dbasis,bc,g,container); // not with pointer


  // // In the elasticity assemble.,
  // gsInfo<<"Material1\n";
  // gsMaterialEval<real_t, gsMaterialOutput::C, true> ply_mat_0_eval(ply_mat_0.get(),&mp); // angle?????
  // gsInfo<<"Material2\n";
  // gsMaterialEval<real_t, gsMaterialOutput::C, true> ply_mat_90_eval(ply_mat_90.get(),&mp); 

  //=============================================//
  //                    Assembly                 //
  //=============================================//
  gsStopwatch clock;
  clock.restart();
  A.assemble();
  gsInfo << "Assembled a system (matrix and load vector) with "
          << A.numDofs() << " dofs in " << clock.stop() << "s.\n";
  gsInfo << "Solving...\n";
  clock.restart();  

  // Solve linear system
  gsSparseSolver<>::SimplicialLDLT solver(A.matrix());
  gsVector<> solVector = solver.solve(A.rhs()); // displacements!?
  gsInfo << "Solved the system with EigenLDLT solver in " << clock.stop() <<"s.\n";

  //=============================================//
  //                    Output                   //
  //=============================================//
  
  // constructing solution as an IGA function
  gsMultiPatch<> sol;
  A.constructSolution(solVector,A.allFixedDofs(),sol); //??? esto no se si esta bien... que resuelve esto?????? son displacements or deformation
  // constructing stresses
  gsPiecewiseFunction<> stresses;
  A.constructCauchyStresses(sol,stresses,stress_components::all_3D_matrix);
  // Mirar que soluciones puedo tener ??

  // gsWriteParaview(mp,"multip");
  // gsDebugVar(mp.nPatches());

  if (numPlotPoints > 0)
  {
        // constructing an IGA field (mp + solution)
        gsField<> solutionField(A.patches(),sol);
        gsField<> stressField(A.patches(),stresses,true);
        // creating a container to plot all fields to one Paraview file
        std::map<std::string,const gsField<> *> fields;
        fields["Displacement"] = &solutionField; // not deformation!!!
        fields["Stresses"] = &stressField;
        gsWriteParaviewMultiPhysics(fields,"composite_elast",numPlotPoints);
        gsInfo << "Open \"composite_elast.pvd\" in Paraview for visualization.\n";
  }

  delete ply_mat_0;
  delete ply_mat_90;
  return EXIT_SUCCESS;

}  // end main
