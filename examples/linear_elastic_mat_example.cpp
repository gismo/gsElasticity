/** @file linear_elasticity_example.cpp

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
*/

#include <gismo.h>
#include <gsElasticity/gsMaterialEval.h>
#include <gsElasticity/gsLinearMaterial.h>
#include <gsElasticity/gsLinearDegradedMaterial.h>
#include <gsElasticity/gsVisitorElUtils.h>

namespace gismo 
{
    namespace expr 
    {
        /**
        Transforms a matrix expression into a vector expression by computing the vector
        [ a b c+d]^T
        for each matrix block
        [ a d ]
        [ c b ]
        */
        template<class E>
        class voigt_expr  : public _expr<voigt_expr<E> >
        {
        public:
            typedef typename E::Scalar Scalar;
            enum {ScalarValued = 0, Space = E::Space, ColBlocks= 0}; // to do: ColBlocks
        private:
            typename E::Nested_t _u;
            mutable gsMatrix<Scalar> tmp, result;
            mutable short_t d;

        protected:
            inline short_t _voigt(short_t dim, short_t I, short_t J) const
            {
                if (dim == 2)
                    switch(I)
                    {
                    case 0: return J == 0 ? 0 : 0;
                    case 1: return J == 0 ? 1 : 1;
                    case 2: return J == 0 ? 0 : 1;
                    }
                else if (dim == 3)
                    switch (I)
                    {
                    case 0: return J == 0 ? 0 : 0;
                    case 1: return J == 0 ? 1 : 1;
                    case 2: return J == 0 ? 2 : 2;
                    case 3: return J == 0 ? 0 : 1;
                    case 4: return J == 0 ? 1 : 2;
                    case 5: return J == 0 ? 0 : 2;
                    }
                GISMO_ERROR("voigt notation indices error");
            }

            inline void _voigtStrain(typename gsMatrix<Scalar>::ColXpr & Evec, const gsMatrix<Scalar> & E_mat) const
            {
                short_t dim = E_mat.cols();
                short_t dimTensor = dim*(dim+1)/2;
                Evec.resize(dimTensor);
                for (short i = 0; i < dimTensor; ++i)
                    if (_voigt(dim,i,0) != _voigt(dim,i,1))
                        Evec(i) = E_mat(_voigt(dim,i,1),_voigt(dim,i,0)) + E_mat(_voigt(dim,i,0),_voigt(dim,i,1)); // off-diagonal terms
                    else
                        Evec(i) = E_mat(_voigt(dim,i,0),_voigt(dim,i,1)); // diagonal terms
            }


        public:

            voigt_expr(_expr<E> const& u) : _u(u)
            {
                d = _u.rows(); 
                //GISMO_ASSERT( _u.rows()*_u.cols() == _n*_m, "Wrong dimension"); //
            }

            const gsMatrix<Scalar> & eval(const index_t k) const
            {
                tmp = _u.eval(k);
                const index_t numActives = _u.cardinality(); // basis functions 
                result.resize(d*(d+1)/2,numActives);

                for (index_t i = 0; i<numActives; ++i)
                {
                    typename gsMatrix<Scalar>::ColXpr col = result.col(i);
                    _voigtStrain(col,tmp.block(0,i*_u.cols(),_u.rows(),_u.cols()));
                }
                
                // Why is this done????
                if ( 1==Space )
                    result.transposeInPlace(); 
                else if (2!=Space) // if not colSpan and not rowSpan
                    result.transposeInPlace();


                return result;
            }

            // Per basis function 
            index_t rows() const { return 1; }
            index_t cols() const { return d*(d+1)/2; } // dim*(dim+1)/2

            void parse(gsExprHelper<Scalar> & evList) const
            { _u.parse(evList); }

            const gsFeSpace<Scalar> & rowVar() const { return _u.rowVar(); }
            const gsFeSpace<Scalar> & colVar() const { return _u.colVar(); }
            index_t cardinality_impl() const { return _u.cardinality_impl(); }

            void print(std::ostream &os) const { os << "flat("; _u.print(os); os<<")"; }
        };

        /// Make a matrix 2x2 expression "flat"
        template <typename E> EIGEN_STRONG_INLINE
        voigt_expr<E> const voigt(E const & u)
        { return voigt_expr<E>(u); }

    }
}

using namespace gismo;

int main(int argc, char *argv[]) {
  //! [Parse command line]
  bool plot = false;
  index_t n_h_refinements = 0;
  index_t n_deg_elevations = 0;
  bool only_last = false;
  bool compute_error{false};
  index_t sample_rate{9};

  std::string file_name("pde/linear_elasticity_example_singlepatch.xml");

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

  gsMultiPatch<> mp, mp_def;
  file_data.getId(0, mp);  // id=0: Multipatch domain

  gsFunctionExpr<> source_function_expression;
  file_data.getId(1, source_function_expression);  // id=1: source function
  gsInfo << "Source function " << source_function_expression << "\n";

  gsBoundaryConditions<> bc;
  file_data.getId(2, bc);  // id=2: boundary conditions
  bc.setGeoMap(mp);
  gsInfo << "Boundary conditions:\n" << bc << "\n";

  gsFunctionExpr<> reference_solution_expr;
  file_data.getId(3, reference_solution_expr);  // id=3: reference solution

  gsOptionList assembler_options = gsAssembler<>::defaultOptions();
  if (file_data.hasId(4)) {
    file_data.getId(4, assembler_options);  // id=4: assembler options
  }
  //! [Read input file]

  //! [Refinement]
  gsMultiBasis<> dbasis(mp, true);
  const int dim = mp.geoDim();

  // Elevate and p-refine the basis
  dbasis.degreeElevate(n_deg_elevations);

  // h-refine each basis
  for (int r = 0; r < n_h_refinements; ++r) {
    dbasis.uniformRefine();
  }

  mp_def = mp;
  //! [Refinement]

  //! [Problem setup]
  gsExprAssembler<> A(1, 1);
  if (file_data.hasId(4)) {
    A.setOptions(assembler_options);
  }

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
  Space w =
      A.getSpace(dbasis, dim);

  // Set the source term
  auto source_function =
      A.getCoeff(source_function_expression, G);

  // Recover manufactured solution
  auto reference_solution =
      ev.getVariable(reference_solution_expr, G);

  // Solution vector and solution variable
  gsMatrix<> solution_vector;
  Solution u =
      A.getSolution(w, solution_vector);
  //! [Problem setup]

  // ![User output]
  gsInfo << "Problem Overview:\t"
         << "\nGeometric Dim:\t" << dim << "\nPatches:\t"
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
  gsSparseSolver<>::CGDiagonal linear_solver;

  real_t l2_err{}, h1_err{}, setup_time{}, ma_time{}, slv_time{}, err_time{};
  gsStopwatch timer;
  dbasis.uniformRefine();

  // Apply dirichlet boundary conditions
  w.setup(bc, dirichlet::l2Projection, 0);

  // Initialize the system
  A.initSystem();
  setup_time += timer.stop();

  // User output
  gsInfo << "Number of degrees of freedom:\t" << A.numDofs()
         << std::endl;
  gsInfo << "\nAssembling linear system ... " << std::flush;
  timer.restart();

  // Compute the system matrix and right-hand side
  auto phys_jacobian_space = ijac(w, G);

  real_t thickness = 1.0;
  real_t E_modulus = 210e9;
  real_t PoissonRatio = 0.0;
  real_t Ratio = 1;
  gsFunctionExpr<> t(std::to_string(thickness), 3);
  gsFunctionExpr<> E(std::to_string(E_modulus),3);
  gsFunctionExpr<> nu(std::to_string(PoissonRatio),3);
  gsConstantFunction<> ratio(Ratio,3);

  //gsFunctionExpr<> fun("if ((x> 0) and (x < 0.5) and (y > 0.5-0.005e-3) and (y < 0.5+0.005e-3)) { 1 } else { 0 };",3); //3 dimension
  
  gsFunctionExpr<> fun("if ((x> -0.05) and (x < 0.501) and (y > 0.5-0.05) and (y < 0.5+0.05)) { 1 } else { 0. };",mp.geoDim()); //3 dimension
  


  //gsFunctionExpr<> fun2("if (( 0 < x < 7.5e-3) and (0.5-0.005e-3 < y < 0.5+0.005e-3)) { 1 } else { 0 };",2);
  gsMatrix<> tmp1;
  gsL2Projection<real_t>::projectFunction(dbasis,fun,mp,tmp1); 
  gsGeometry<>::uPtr geom_ptr = dbasis.basis(0).makeGeometry(give(tmp1));   
  gsWriteParaview(mp,*geom_ptr,"fun",1e5);
       

  //gsConstantFunction<> fun(0,3); // to set damage to zero!

  // gsGeometry<>::uPtr geom_ptr = dbasis.basis(0).makeGeometry(give(coefs(dbasis.size(),2)));        
  // gsWriteParaview(*geom_ptr,mp,"fun");
  // gsMatrix<> tmp;
  // gsL2Projection<real_t>::projectFunction(dbasis,fun,mp,tmp); 

  gsLinearDegradedMaterial<3,real_t> SvK(E_modulus, PoissonRatio, mp, mp_def, fun);

  gsWriteParaview(mp,fun,"damage",1e5);

  //gsLinearMaterial<3,real_t> SvK(E_modulus, PoissonRatio, mp, mp_def); 
  gsMaterialEval<real_t, gsMaterialOutput::E, true> SvK_E(&SvK, &mp_def);
  gsMaterialEval<real_t, gsMaterialOutput::S, true> SvK_S(&SvK, &mp_def);
  gsMaterialEval<real_t, gsMaterialOutput::C, true> SvK_C(&SvK, &mp_def);
  //gsMaterialEval<real_t, gsMaterialOutput::Psi, true> SvK_Psi(&SvK, &mp_def);

  gsLinearMaterial<3,real_t> SvK_undeg(E_modulus, PoissonRatio, mp, mp_def); 
  gsMaterialEval<real_t, gsMaterialOutput::C, true> SvK_undeg_C(&SvK_undeg, &mp_def);


  //  =========== Bilinear form I ===========
  real_t lambda = E_modulus * PoissonRatio / ( ( 1. + PoissonRatio ) * ( 1. - 2. * PoissonRatio ) );
  real_t mu     = E_modulus / ( 2. * ( 1. + PoissonRatio ) );
  real_t kappa  = E_modulus / ( 3. * ( 1. - 2.*PoissonRatio ) );
  
  auto bilin_lambda = lambda * idiv(w, G) *
                      idiv(w, G).tr() * meas(G);
  auto bilin_mu_ =
      mu *
      ((phys_jacobian_space.cwisetr() + phys_jacobian_space) % phys_jacobian_space.tr()) *
      meas(G);
  auto bilin_combined2 = (bilin_lambda + bilin_mu_);

  // =========== Bilinear form II ===========
  auto delta_EE = 0.5*(phys_jacobian_space.cwisetr() + phys_jacobian_space);
  auto delta_EE_voigt = voigt(delta_EE);
  auto C = A.getCoeff(SvK_C);
  auto C_elast = A.getCoeff(SvK_undeg_C);

  //auto energy_pos = A.getCoeff(SvK_Psi);

  // A.clearMatrix();
  auto bilin_combined = delta_EE_voigt*reshape(C,6,6)*delta_EE_voigt.tr()*meas(G);
  
  // A.assemble(bilin_combined);
  
  // gsMatrix<> matrix_gsMat = A.matrix().toDense();
  // A.clearMatrix();

  gsVector<> pt(3);
  //pt.col(0)<<0.125,0.375,0.5;
  pt.col(0)<<0.,0.5,1.;

  gsDebugVar(ev.eval(reshape(C,6,6),pt,0)); // degraded C
  gsDebugVar(ev.eval(reshape(C_elast,6,6),pt,0));

  //gsDebugVar(ev.eval(energy_pos,pt,0));

  //hola

  // gsDebugVar(ev.eval(bilin_combined,pt,0)); //bilinear de gsmaterial
  // gsDebugVar(ev.eval(bilin_combined2,pt,0)); //bilinear de linear elasticity

  // A.assemble(bilin_combined2);
  // gsMatrix<> matrix2_elast = A.matrix().toDense();

  //gsDebugVar((matrix_gsMat-matrix2_elast).maxCoeff());

  // Linear Form
  auto linear_form = w * source_function * meas(G);

  A.assemble(bilin_combined,  // matrix
                                linear_form      // rhs vector
  );

  // Compute the Neumann terms (will not be used in default file)
  auto neumann_boundary_function = A.getBdrFunction(G);

  // Neumann conditions
  A.assembleBdr(
      bc.get("Neumann"),
      w * neumann_boundary_function * nv(G).norm());
  //! [Assembly]

  //! [Solver]
  ma_time += timer.stop();
  gsInfo << "Done" << std::endl;
  gsInfo << "Solving Linear System ...... " << std::flush;
  timer.restart();

  // Solve linear system
  linear_solver.compute(A.matrix());
  solution_vector = linear_solver.solve(A.rhs());
  //! [Solver]

  slv_time += timer.stop();
  gsInfo << "Done" << std::endl;

  //! [Evaluate Error]
  if (compute_error) {
    gsInfo << "Evaluating solution ........ " << std::flush;
    timer.restart();
    l2_err = math::sqrt(
        ev.integral((reference_solution - u).sqNorm() *
                    meas(G)));

    h1_err = l2_err +
             math::sqrt(ev.integral((igrad(reference_solution) -
                                     igrad(u, G))
                                        .sqNorm() *
                                    meas(G)));
    err_time += timer.stop();
    gsInfo << "Done" << std::endl;
  }
  //! [Evaluate Error]

  //! [Summarize User Output]
  gsInfo << "\n\nTotal time: " << setup_time + ma_time + slv_time + err_time
         << "\n";
  gsInfo << "     Setup: " << setup_time << "\n";
  gsInfo << "  Assembly: " << ma_time << "\n";
  gsInfo << "   Solving: " << slv_time << "\n";
  gsInfo << "     Norms: " << err_time << "\n";

  if (compute_error) {
    gsInfo << "\nL2 error: " << std::scientific << std::setprecision(3)
           << l2_err << "\nH1 error: " << std::scientific << h1_err << "\n\n";
  }
  //! [Summarize User Output]

  //! [Export visualization in ParaView]
  if (plot) {
    gsInfo << "Plotting in Paraview ... ";

    gsParaviewCollection collection("ParaviewOutput/solution", &ev);
    collection.options().setSwitch("plotElements", true);
    collection.options().setInt("plotElements.resolution", sample_rate);
    collection.newTimeStep(&mp);
    collection.addField(u, "numerical solution");
    if (compute_error) {
      collection.addField(reference_solution - u, "error");
    }
    collection.saveTimeStep();
    collection.save();

    gsFileManager::open("ParaviewOutput/solution.pvd");
    gsInfo << "Done" << std::endl;
  } else {
    gsInfo << "No output created, re-run with --plot to get a ParaView "
              "file containing the solution.\n";
  }
  //! [Export visualization in ParaView]

  return EXIT_SUCCESS;

}  // end main
