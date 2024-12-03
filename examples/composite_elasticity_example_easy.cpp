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
  ./bin/composite_elasticity_example -r 3

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
  index_t numPlotPoints = 5;


  std::string file_name("composite_example_singlepatch.xml");

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
  // real_t lame_lambda{80000.0}, lame_mu{80000.0};
  // cmd.addReal("L", "firstLame", "First Lame constant, material parameter",
  //             lame_lambda);
  // cmd.addReal("M", "secondLame", "Second Lame constant, material parameter",
  //             lame_mu);

  //! [Parse command line]
  try {
    cmd.getValues(argc, argv);
  } catch (int rv) {
    return rv;
  }

  //! [Read input file]
  gsDebugVar(file_name);
  gsFileData<> file_data(file_name);
  gsInfo << "Loaded file " << file_data.lastPath() << "\n";

  gsMultiPatch<> mp, mp_def;

  file_data.getId(0, mp);  // id=0: Multipatch domain

  gsBoundaryConditions<> bc;
  //file_data.getId(1, bc);  // id=2: boundary conditions
  
  // Dirichlet BCs
  // Dirichlet x
  bc.addCondition(0,boundary::east,condition_type::dirichlet,0,0);
  bc.addCondition(0,boundary::west,condition_type::dirichlet,0,0);
  // Dirichlet y
  bc.addCondition(0,boundary::north,condition_type::dirichlet,0,1);
  bc.addCondition(0,boundary::south,condition_type::dirichlet,0,1);
  // Dirichlet z
  bc.addCondition(0,boundary::north,condition_type::dirichlet,0,2);
  bc.addCondition(0,boundary::south,condition_type::dirichlet,0,2);
  bc.addCondition(0,boundary::east,condition_type::dirichlet,0,2);
  bc.addCondition(0,boundary::west,condition_type::dirichlet,0,2);


  // Apply Neumann BC
  gsFunctionExpr<>func("0","0","sin(pi*x/60)*sin(pi*y/60)",mp.parDim()); // tiene sentido?????
  bc.addCondition(0,boundary::back,condition_type::neumann,&func);

  bc.setGeoMap(mp);
  gsInfo << "Boundary conditions:\n" << bc << "\n";

  //! [Refinement]
  gsMultiBasis<> dbasis(mp, true);
  const int solution_field_dimension = mp.geoDim();

  // Elevate and p-refine the basis
  // dbasis.degreeElevate(n_deg_elevations);
  // dbasis.degreeElevate(5,0); // degree 6 in plane
  // dbasis.degreeElevate(5,1); // degree 6 in plane
  // dbasis.degreeElevate(3,2); // degree 4 through the thickness
  // dbasis.degreeElevate(1,0); // degree 6 in plane
  // dbasis.degreeElevate(1,1); // degree 6 in plane
  // dbasis.degreeElevate(1,2); // degree 4 through the thickness

  // h-refine each basis
  dbasis.uniformRefine(20,1,0); // 5 elements in 0 direction
  dbasis.uniformRefine(20,1,1); // 5 elements in 1 direction

  // dbasis.uniformRefine(2,1,0); // 5 elements in 0 direction
  // dbasis.uniformRefine(2,1,1); // 5 elements in 1 direction
  // for (int r = 0; r < n_h_refinements; ++r) {
  //   dbasis.uniformRefine();
  // }

  mp_def = mp;
  gsDebugVar(mp_def);
  
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
  gsMaterialBase<real_t> * ply_mat_homogeneous;
  gsMaterialBase<real_t> * ply_mat_elast;

  if      (mp.geoDim()==2)
  {
    ply_mat_0  = new gsCompositeMat<2,real_t>(E11,E22,G12,nu12,rot_0,mp);
    ply_mat_90 = new gsCompositeMat<2,real_t>(E11,E22,G12,nu12,rot_90,mp);
  }
  else if (mp.geoDim()==3)
  {
    ply_mat_homogeneous  = new gsCompositeMat<3,real_t>(true,rot_0,mp);
    ply_mat_elast  = new gsLinearMaterial<3,real_t>(E11,0.25,mp,mp_def);
    ply_mat_0  = new gsCompositeMat<3,real_t>(E11,E22,E33,G12,G23,G13,nu12,nu23,nu13,rot_0,mp);
    ply_mat_90 = new gsCompositeMat<3,real_t>(E11,E22,E33,G12,G23,G13,nu12,nu23,nu13,rot_90,mp);
  }     
  else
      GISMO_ERROR("Dimension not supported");
  // ======================================

  // Option 2
  gsDebugVar(mp.nPatches());
  //gsMaterialContainer<real_t> container(ply_mat_homogeneous,mp.nPatches());
  gsMaterialContainer<real_t> container(ply_mat_homogeneous,mp.nPatches());

  gsConstantFunction<> g(0,0,0,mp.parDim()); // it is a gsFunction (body load?) comprobar??????? no entiendo!!

  // Call the elasticity assembler!
  gsElasticityAssembler<real_t> A(mp,dbasis,bc,g,container); // not with pointer

  //=============================================//
  //                    Assembly                 //
  //=============================================//
  gsStopwatch clock;
  clock.restart();
  gsInfo << "Starting assembly...\n";
  A.assemble();
  gsInfo << "Assembled a system (matrix and load vector) with "
          << A.numDofs() << " dofs in " << clock.stop() << "s.\n";
  gsInfo << "Solving...\n";
  clock.restart();  

  // // Solve linear system
  gsSparseSolver<>::SimplicialLDLT solver(A.matrix());
  gsVector<> solVector = solver.solve(A.rhs()); // displacements!?
  gsInfo << "Solved the system with EigenLDLT solver in " << clock.stop() <<"s.\n";

  //=============================================//
  //                    Output                   //
  //=============================================//
  
  // For the post-processing
  index_t materialLaw = material_law::hooke;
  // index_t materialLaw = material_law::composite;
  A.options().setInt("MaterialLaw",materialLaw);
  gsDebugVar(A.options());

  // constructing solution as an IGA function
  gsMultiPatch<> sol;
  A.constructSolution(solVector,A.allFixedDofs(),sol); //??? esto no se si esta bien... que resuelve esto?????? son displacements or deformation
  // constructing stresses
  gsPiecewiseFunction<> stresses;
  //A.constructCauchyStresses(sol,stresses,stress_components::all_3D_matrix); // for hooke!
  A.constructCauchyStresses(sol,container,stresses,stress_components::all_3D_matrix);
  // Mirar que soluciones puedo tener ??

  // Sol is the multipatch with the solution!
  gsMaterialEval<real_t,gsMaterialOutput::E, true> strains(container,&sol); // ????
  gsInfo<<"Strain\n";
  gsMaterialEval<real_t,gsMaterialOutput::S, true> stresses_para(container,&sol); // ????

  // short_t n = mp.targetDim();
  // short_t d = mp.domainDim();

  // // gsMatrix<> ab = mp.support();
  // gsDebugVar(n);
  // gsDebugVar(d);

  // gsVector<> a = ab.col(0);
  // gsVector<> b = ab.col(1);
  // // gsVector<unsigned> np = distributePoints<real_t>(mp,10000);

  // gsVector<> a(3), b(3);
  // a.col(0)<<0.,0.,0.;
  // b.col(0)<<1.,1.,1.;

  // gsMatrix<> pts = gsPointGrid(a,b,10000);


  // gsVector<> pt(3);
  // pt.col(0)<<0.,0.,0.;
  // strains.piece(0).eval(pt);
  //gsDebugVar(strains.piece(0).eval(pt)); // degraded C

  // gsWriteParaview(strains,"strains"); // how to make difference in eps_xx eps yy etc!


  //gsMaterialEvalMatrix<real_t,gsMaterialOutput::E, true> strains(container,&sol); // ????

  // gsWriteParaview(mp,"multip");
  // gsDebugVar(mp.nPatches());

  // gsInfo<<"HOLA BEFORE WRITING PARAVIEW\n";

  if (numPlotPoints > 0)
  {
        // constructing an IGA field (mp + solution)
        gsField<> solutionField(A.patches(),sol);
        gsField<> stressField(A.patches(),stresses,true);
        gsField<> stressField2(A.patches(),stresses_para,true);
        gsField<> strainField(A.patches(),strains,true);
        // creating a container to plot all fields to one Paraview file
        std::map<std::string,const gsField<> *> fields;
        fields["Displacement"] = &solutionField; // not deformation!!!
        fields["Stresses"] = &stressField;
        //fields["Strains"] = &strainField;
        gsInfo<<"====================================\n"; 
        gsInfo<<"         Stresses_evaluator         \n"; 
        gsInfo<<"====================================\n"; 
        gsWriteParaview(stressField,"stresses",numPlotPoints,true);
        gsWriteParaview(stressField2,"stresses2",numPlotPoints,true);
        gsInfo<<"==============================================\n"; 
        gsInfo<<"         Stresses ElasticityFunctions         \n";
        gsInfo<<"==============================================\n"; 
        gsDebugVar(numPlotPoints);
        //gsWriteParaviewMultiPhysics(fields,"composite_elast",numPlotPoints,true);
        gsInfo << "Open \"composite_elast.pvd\" in Paraview for visualization.\n";
  }

  delete ply_mat_0;
  delete ply_mat_90;
  delete ply_mat_homogeneous;
  delete ply_mat_elast;
  return EXIT_SUCCESS;

}  // end main
