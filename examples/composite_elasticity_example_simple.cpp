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
#include <gsElasticity/gsCompositeMat.h>
#include <gsElasticity/gsHomogenizedCompositeMat.h>

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
  index_t numPlotPoints = 5000; //!!


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

  // p-refinement
  dbasis.degreeElevate(1); // quadratic (for integragtion?????)
  gsInfo<<"Maximum degree of the basis: "<<dbasis.basis(0).maxDegree()<<"\n";

  // h-refine each basis (knot insertion)
  dbasis.uniformRefine(20,1,0); // 21 elements in 0 direction
  dbasis.uniformRefine(20,1,1); // 21 elements in 1 direction

  mp_def = mp;

  //=============================================//
  //                    Material                 //
  //=============================================//  
  real_t E11 = 25000; // in MPa (N/mm^2)
  real_t E22 = E11/25;
  real_t E33 = E11/25;
  real_t G23 = 500;
  real_t G12 = G23/2.5;
  real_t G13 = G23/2.5;
  real_t nu12 = 0.25;
  real_t nu13 = 0.25;
  real_t nu23 = 0.25;

  gsConstantFunction<> rot_0(0.0, 0.0, 0.0, mp.parDim()); // 3 is the dimension?
  gsConstantFunction<> rot_90(0.0, 0.0, M_PI/2.0, mp.parDim()); // in radians?? or d
  gsFunctionExpr<> ply_angle("0","0","if ((z>1) and (z<2)) {0.0} else {1.57079632679};",mp.geoDim()); // for ply 90/0/90

  gsMaterialBase<real_t> * material_homogenized;
  gsMaterialBase<real_t> * material_laminate;

  material_homogenized  = new gsHomogenizedCompositeMat<3,real_t>(true,rot_0,mp); // ideally computation of the homogenized material matrix will be inside the class!
  material_laminate  = new gsCompositeMat<3,real_t>(E11,E22,E33,G12,G23,G13,nu12,nu23,nu13,ply_angle,mp);

  gsMaterialContainer<real_t> container_homogenized(material_homogenized,mp.nPatches());
  gsMaterialContainer<real_t> container_laminate(material_laminate,mp.nPatches());

  // Check if the matrix C changes with angle (atencion a matrix m_C_hom and m_C...)
  gsMaterialEval<real_t, gsMaterialOutput::C, true> C_mat(container_laminate, &mp_def);

  gsMatrix<> pts(3,3);
  pts.col(0)<< 0,0,0.1;
  pts.col(1)<< 0,0,1.1;
  pts.col(2)<< 0,0,2.1;

  gsMatrix<> C1, C2, C3;
  C_mat.piece(0).eval_into(pts.col(0),C1);
  C_mat.piece(0).eval_into(pts.col(1),C2);
  C_mat.piece(0).eval_into(pts.col(2),C3);

  //=============================================//
  //                    Assembly                 //
  //=============================================//
  
  gsConstantFunction<> g(0,0,0,mp.parDim()); // body load

  gsElasticityAssembler<real_t> A(mp,dbasis,bc,g,container_homogenized); // not with pointer

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

  gsDebugVar(solVector.cols());  
  gsDebugVar(solVector.rows());

  gsDebugVar(mp.patch(0).coefs().cols());
  gsDebugVar(mp.patch(0).coefs().rows());

  // ===== Update coefficients multipatch? =====
  //mp_def.patch(0).coefs() = mp.patch(0).coefs() + solVector; // deformation *2 in the x direction (col(0))

  // Construct solution as an IGA function
  gsMultiPatch<> sol;
  A.constructSolution(solVector,A.allFixedDofs(),sol); //??? esto no se si esta bien... que resuelve esto?????? son displacements or deformation
  gsField<> solutionField(A.patches(),sol);
  // gsWriteParaview(sol,"displacement",numPlotPoints,true);
  std::map<std::string,const gsField<> *> fields;
  fields["Displacement"] = &solutionField; // not deformation!!!
  gsWriteParaviewMultiPhysics(fields,"composite_disp",numPlotPoints,true);

  // gsMaterialEval<real_t,gsMaterialOutput::S, true> stresses(container_homogenized,&sol); // ????
  // gsField<> stressField(A.patches(),stresses,true);
  // gsWriteParaview(stressField,"stresses",numPlotPoints,true); // add mesh!!!!!
  // gsWriteParaview(mp, stresses.piece(0).coord(0),"stressesxx"); // stress sigma_11 with the homogenized matrix!

  //=============================================//
  //     Post-processing (stress recovery)       //
  //=============================================//

  gsMaterialEval<real_t,gsMaterialOutput::S, true> stresses(container_laminate,&sol); // ????

  // gsMatrix<> s1, s2, s3;
  // stresses.piece(0).eval_into(pts.col(0),s1);
  // stresses.piece(0).eval_into(pts.col(1),s2);
  // stresses.piece(0).eval_into(pts.col(2),s3);

  std::vector <index_t> index = {0,1,3}; // indices of in-plane stresses 

  gsMatrix<> coefs(dbasis.basis(0).size(),index.size()), tmp;

  for(index_t i=0; i<index.size(); i++)
  {
    // QI the stress values index[i] (position in Voigt stress)
    gsQuasiInterpolate<real_t>::localIntpl(dbasis.basis(0),stresses.piece(0).coord(index[i]),tmp);  
    coefs.col(i) = tmp.col(0); //save the result in coefs
    // stress_plane[i]=dbasis.basis(0).makeGeometry(give(coefs));
  }  

  gsGeometry<>::uPtr stress_plane =  dbasis.basis(0).makeGeometry(give(coefs)); // create spline function with coefficients ()
  gsWriteParaview(mp, *stress_plane,"stress_qi", numPlotPoints);
  // Check the values of ths in plane stresses! are they correct? 

  // ================== Integration ==================
  // Quadrature rule: selection of points???? Per ply????????????

  //Gauss Legendre
  gsGaussRule<real_t> rule(6); // number of quadrature points!!!!!!!
  gsMatrix<> quPoints;
  gsVector<> quWeights;
  rule.mapTo(0.0,3.0,quPoints,quWeights); 
  // gsDebugVar(quPoints);
  // gsDebugVar(quWeights);

  gsMatrix<> dS_in_plane;
  stress_plane->deriv_into(quPoints,dS_in_plane); // in plane stresses evaluated at quPoints
  // Body load!?!?!?!?

  gsMatrix<> dS11_1 = dS_in_plane.row(0); // not sure .row()
  gsMatrix<> dS22_2 = dS_in_plane.row(1);
  gsMatrix<> dS12_2 = dS_in_plane.row(2); // no se 
  gsMatrix<> dS12_1 = dS_in_plane.row(2); // no se 

  gsMatrix<> S13 = - (dS11_1 + dS12_2) * quWeights; // missing body force and sigma at 0?????????
  gsMatrix<> S23 = - (dS12_1 + dS22_2) * quWeights; // missing body force and sigma at 0?????????

 // Create spline function with coefficients of out of plane stresses (sigma 23 and sigma 13)
  gsGeometry<>::uPtr stress_13 =  dbasis.basis(0).makeGeometry(give(S13)); 
  gsGeometry<>::uPtr stress_23 =  dbasis.basis(0).makeGeometry(give(S23)); 
  // gsWriteParaview(mp, *stress_13,"stress_13", numPlotPoints);
  // gsWriteParaview(mp, *stress_23,"stress_23", numPlotPoints);

  gsMatrix<> dS13_1, dS23_2;
  stress_13->deriv_into(quPoints,dS13_1); // in plane stresses evaluated at quPoints
  stress_23->deriv_into(quPoints,dS23_2); // in plane stresses evaluated at quPoints

  gsMatrix<> S33 = - (dS13_1 + dS23_2) * quWeights;

  gsGeometry<>::uPtr stress_33 =  dbasis.basis(0).makeGeometry(give(S33)); 

  // Then --> plot stresses (gsWriteParaview)


  //Integration in the quasi interpolation
  //   gsMatrix<T> bev, fev, pts, tmp;
  //   gsVector<T> weights;
  //   gsVector<index_t> nNodes = gsQuadrature::numNodes(bb,(T)1.0,1);
  //   gsQuadRule<T>  qRule     = gsQuadrature::get<T>(gsQuadrature::GaussLegendre,nNodes);
  //   qRule.mapTo(ab.col(0),ab.col(1), pts, weights);//map points and weights on element (for quadrature)
    
  //   bb .eval_into(pts, bev);//evaluate basis on quadrature points (numbasis*pts)
  //   fun.eval_into(pts, fev);//evaluate function on quadrature points (1*pts)
  //   gsMatrix<T> M = bev * weights.asDiagonal() * bev.transpose(); // Mass matrix
  //   gsMatrix<T> RHS = bev * weights.asDiagonal() * fev.transpose();
  //   tmp = M.fullPivLu().solve(RHS); //solve on element


  //   // find the i-th BS:
  //   gsMatrix<index_t> act = bb.active(pts.col(0)); // the same basis functions will be active for the element in consideration!
  //   index_t c = std::lower_bound(act.data(), act.data()+act.size(), i) - act.data();
  //   GISMO_ASSERT(c<act.size(), "Problem with basis function index");
  //   return tmp.row(c);

  // // Gauss Legendre
  // gsOptionList legendreOpts;
  // legendreOpts.addInt   ("quRule","Quadrature rule used (1) Gauss-Legendre; (2) Gauss-Lobatto; (3) Patch-Rule",gsQuadrature::GaussLegendre);
  // legendreOpts.addReal("quA", "Number of quadrature points: quA*deg + quB", 1.0  );
  // legendreOpts.addInt ("quB", "Number of quadrature points: quA*deg + quB", 1    );
  // legendreOpts.addSwitch("overInt","Apply over-integration or not?",false);
  // gsQuadRule<real_t>::uPtr legendre = gsQuadrature::getPtr(dbasis.basis(0), legendreOpts);


  // gsMatrix<> pts(3,rule.nPoints());
  // pts.row(0).setConstant(pt[0]);
  // pts.row(1).setConstant(pt[1]);
  // pts.row(2) = quPoints.row(0);

  // gsMatrix<> dS = SvK_S.deriv(pts);

  // gsMatrix<> dS_int = dS * quWeights;


  // gsExprAssembler<> A2(1,1);
  // A2.setIntegrationElements(dbasis);
  // gsExprEvaluator<> ev(A2);
  // auto G = ev.getMap(mp); // is this correct?
  // auto c_sinus_L2_local_fine = ev.getVariable(*fine_L2_local); // pointer to the L2 local


  // For the approximation of the stresses!
  // gsMatrix<> coefs;
  // gsQuasiInterpolate<real_t>::localIntpl(mp.basis(0),SvK_S,coefs);
  // gsGeometry<>::uPtr geo = mp.basis(0).makeGeometry(coefs);
  // gsWriteParaview(mp_def.patch(0),SvK_S,"S");
  // gsWriteParaview(mp_def.patch(0),*geo,"S_int");

  // gsVector<> pt(2);
  // pt.col(0)<<0.125,0.375;

  // Gauss Legendre
  // gsGaussRule<real_t> rule(4);
  // gsMatrix<> quPoints;
  // gsVector<> quWeights;
  // rule.mapTo(0.0,1.0,quPoints,quWeights);

  // gsMatrix<> pts(3,rule.nPoints());
  // pts.row(0).setConstant(pt[0]);
  // pts.row(1).setConstant(pt[1]);
  // pts.row(2) = quPoints.row(0);

  // gsMatrix<> dS = SvK_S.deriv(pts);

  // gsMatrix<> dS_int = dS * quWeights;



  delete material_homogenized;
  delete material_laminate;
  return EXIT_SUCCESS;

}  // end main
