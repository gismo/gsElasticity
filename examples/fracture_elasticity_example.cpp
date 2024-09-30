/** @file fracture_elasticity_example.cpp

    @brief Tutorial on how to use expression assembler to solve the Cahn-Hilliard equation

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): 
    
    ./bin/fracture_elasticity_example -f linear_elasticity_example_singlepatch_2d.xml -r 8  --plot
*/

//! [Include namespace]
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
//! [Include namespace]

int main(int argc, char *argv[])
{
    //! [Parse command line]
    bool plot = false;
    index_t n_h_refinements = 0;
    index_t n_deg_elevations = 0;
    bool only_last = false;
    bool compute_error{false};
    index_t sample_rate{9}; 

    std::string file_name("linear_elasticity_example_singlepatch.xml");

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

    try { cmd.getValues(argc,argv); } catch (int rv) { return rv; }
    
    // ========== [Read input file] ==========
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
    // =====================================

    // ========== [Refinement] ==========
    // gsMultiBasis<> dbasis(mp, true);
    // const int dim = mp.geoDim();

    // // Elevate and p-refine the basis
    // dbasis.degreeElevate(n_deg_elevations);

    // // h-refine each basis
    // for (int r = 0; r < n_h_refinements; ++r) {
    //     dbasis.uniformRefine();
    //     mp.uniformRefine();
    // }

    // p-refine
    if (n_deg_elevations!=0)
        mp.degreeElevate(n_deg_elevations);

    // h-refine
    for (int r =0; r < n_h_refinements; ++r)
        mp.uniformRefine();
    
    mp_def = mp;
    // gsWriteParaview(mp,"mp",1000,true);

    gsMultiBasis<> dbasis(mp);
    // =====================================

    // ========== [Problem setup] ==========
    //gsExprAssembler<> A(2,2); // 2 trial 2 test spaces
    gsExprAssembler<> A_u(1,1); // 1 trial 1 test space
    gsExprAssembler<> A_d(1,1); // 1 trial 1 test space

    typedef gsExprAssembler<>::geometryMap geometryMap;
    typedef gsExprAssembler<>::variable    variable;
    typedef gsExprAssembler<>::space       space;
    typedef gsExprAssembler<>::solution    solution;

    // Elements used for numerical integration
    A_u.setIntegrationElements(dbasis); //?? makes sense? no?
    A_d.setIntegrationElements(dbasis); //?? different basis?
   
   // gsExprEvaluator<> ev(A);
    gsExprEvaluator<> ev_u(A_u);
    gsExprEvaluator<> ev_d(A_d);

    // Set the geometry map
    geometryMap G = A_u.getMap(mp); //????? no estoy seguro??? necesito mp or mp_u?

    // Set the discretization space (!!!!!!!!!!different for displacements and phase field!!!!!!!!!)
    space w = A_u.getSpace(dbasis,mp.parDim()); // displacements (3rd argument is the id of the space)
    space b = A_d.getSpace(dbasis,1); // phase field d  (es parDim-1?????????? no estoy seguro)

    // Solution vectors and solution variables
    gsMatrix<> Unew, Uold;
    gsMatrix<> Dnew, Dold;
    gsMatrix<> deltaU, deltaD; // ????? no estoy seguro

    solution uold = A_u.getSolution(w, Uold); 
    solution unew = A_u.getSolution(w, Unew); 

    solution dold = A_d.getSolution(b, Dold); 
    solution dnew = A_d.getSolution(b, Dnew); 

    gsMultiPatch<> mp_unew, mp_dnew;

    auto dnew_u = A_u.getCoeff(mp_dnew); // gets mp as an evaluable object
    auto unew_d = A_d.getCoeff(mp_unew);

    // Setup the spaces (compute Dirichlet BCs)
    w.setup(bc, dirichlet::l2Projection, 0); // computes the values on CPs of dirichlet BC   
    b.setup(); // no damage boundary conditions 

    Dold.setZero(A_d.numDofs(),1);
    Uold.setZero(A_u.numDofs(),1);
    
    // Extract variable to the corresponding multipatch
    unew.extract(mp_unew); 
    //dnew.extract(mp_dnew); // I do this in the initial condition!

    // =====================================

    // Define linear solver (install SuperLUMT-devel)
#ifdef GISMO_WITH_SUPERLU
    gsSparseSolver<>::SuperLU solver;
#   else
    gsSparseSolver<>::LU solver;
#endif
    
    // Generalized-alpha method parameters
    // real_t rho_inf = 0.5;
    // real_t alpha_m = 0.5*(3-rho_inf) / (1+rho_inf);
    // real_t alpha_f = 1 / (1+rho_inf);
    // real_t gamma   = 0.5 + alpha_m - alpha_f;

    // Legend:
    // C_old   = C_n
    // C       = C_n+1,i-1
    // C_alpha = C_{n+alpha_f,i-1}
    // dC

    gsInfo<<"Starting.."<<"\n";
    
    gsInfo<<"Initial condition.."<<"\n";

    // Compute jacobian
    auto phys_jacobian_space = ijac(w, G);

    // AT2 parameters                                        
    real_t c_w = 2.0;
    // Material properties
    real_t G_c = 2.7e3; // in N/m
    real_t ell = 0.015; // in m
    real_t E_modulus = 2.1e11; // Young's modulus
    real_t PoissonRatio = 0.3; // Poisson's ratio
    // Formulas 3D
    real_t mu = E_modulus/(2.0*(1+PoissonRatio)); // shear modulus 
    real_t kappa = E_modulus/(3.0*(1-2.0*PoissonRatio)); // bulk modulus

    // ================ Initialize displacements and phase field ================
    //Uold.setZero(A.numDofs(),mp.parDim());


    // 1) Knot coordinates x and y ?? if not in parametric space?
    // 2) Check if values of x and y are within the bounds [0, 7.5e-3] and [0.5-l0/3, 0.5+l0/3], respectively
    // 3) Assign value of 1 
    // 4) L2 project to the control points (how do I L2 project without a source function)
    // gsL2Projection<real_t>::projectFunction(dbasis,source,mp,tmp); 
    gsFunctionExpr<> fun("if ((x> -0.05) and (x < 0.501) and (y > 0.5-0.01) and (y < 0.5+0.01)) { 1. } else { 0. };",mp.geoDim()); //3 dimension
    gsMatrix<> anchors = dbasis.basis(0).anchors();
    gsMatrix<> vals = fun.eval(anchors);
    vals.transposeInPlace();
    //gsDebugVar(vals);
    gsGeometry<>::uPtr geom_ptr = dbasis.basis(0).makeGeometry(give(vals));   
    gsWriteParaview(mp,*geom_ptr,"fun2",1e5);

    mp_dnew.addPatch(*geom_ptr); // add a patch to mp_dnew
    //gsDebugVar(mp_dnew.patch(0).coefs());
    dold.insert(mp_dnew.patch(0),0); // CHECK! update the coefficients in the solution!
    // evaluation to check that dnew equal mp_dnew (solution updated to initial condition)
    //gsDebugVar(Dold);
    // New solution (check!)
    Dnew = Dold; 
    Unew = Uold;
    // ============================================================================

    // ========== Material class ==========
    gsMaterialBase<real_t>::uPtr SvK;
    gsMaterialBase<real_t>::uPtr SvK_undeg;
    if      (mp.geoDim()==2)
    {
        SvK = memory::make_unique (new gsLinearDegradedMaterial<2,real_t>(E_modulus, PoissonRatio, mp, mp_def, mp_dnew));
        SvK_undeg = memory::make_unique (new gsLinearMaterial<2,real_t>(E_modulus, PoissonRatio, mp, mp_def));
    }
    else if (mp.geoDim()==3)
    {
        SvK = memory::make_unique (new gsLinearDegradedMaterial<3,real_t>(E_modulus, PoissonRatio, mp, mp_def, mp_dnew));
        SvK_undeg = memory::make_unique (new gsLinearMaterial<3,real_t>(E_modulus, PoissonRatio, mp, mp_def));
    }
    else
        GISMO_ERROR("Dimension not supported");
    // ======================================

    // Loop variables
    real_t res_u_norm0, res_u_norm, res_d_norm0, res_d_norm, res_stag_norm0, res_stag_norm;
    // real_t res_norm0 = 1, res_norm = 10;
    real_t tol_NR = 1e-6;
    real_t tol_stag = 1e-4;
    index_t maxStag = 10; // Staggered iterations
    index_t maxSteps = 10; // Time steps
    index_t maxIt = 50; // NR iterations
    real_t penalty_irrev = 2.7e12; // Penalty formulation (T. Gerasimov et al., 2019)

    // gsParaviewCollection collection("ParaviewOutput/solution", &ev);
    gsParaviewCollection collection("ParaviewOutput/solution", &ev_u); // Which evaluator?
    collection.options().setSwitch("plotElements", true);
    collection.options().setInt("plotElements.resolution", 4);
    collection.options().setInt("numPoints", 10000); 

    real_t time = 0;
    bool converged = false;

    // A.initSystem(); // Initialize the system (outside the loops)
    A_u.initSystem(); 
    A_d.initSystem();

    gsSparseMatrix<> K_u, K_d, K_d_old;
    gsMatrix<> Q_u, Q_d;

    // Voigt dimension
    short_t sz = (mp.parDim()+1)*mp.parDim()/2; 

    // umax =  0.04 mm
    // delta_u = 5e-3 mm

    // These work by reference (updating mp_def already updates C/Psi...)
    gsMaterialEval<real_t, gsMaterialOutput::C, true> SvK_C(SvK.get(), &mp_def); // change name of mp?
    gsMaterialEval<real_t, gsMaterialOutput::Psi, true> SvK_Psi(SvK.get(), &mp_def);

    //for (index_t step = 0; step!=maxSteps; step++)
    for (index_t step = 0; step!=1; step++)
    {
        // gsInfo<<"Time step "<<step<<"/"<<maxSteps<<", iteration "<<dt_it<<": dt = "<<dt<<", [t_start,t_end] = ["<<time<<" , "<<time+dt<<"]"<<"\n";
        gsInfo<<"Time step "<<step<<"/"<<maxSteps<<"\n";
        res_u_norm0 = 1;
        res_d_norm0 = 1;
        res_stag_norm0 = 1;
        res_u_norm = 10;
        res_d_norm = 10;
        res_stag_norm = 10;

        // !!!!!!!!!!!!!!!!!!!!!!!!! STEP THE APPLIED DISPLACEMENT !!!!!!!!!!!!!!!!!!!!!!!!!
        //Unew = Uold + (step/maxSteps) * Uold; // can I divide by an index?
        //Dnew += step/maxSteps;

        // If (Dirichlet)
        //     Dnew += step/maxSteps*(applieddisplacement);
        //     Dnew += deltaU;
        // If (neumann)
                // increase value of traction;
                // Assemble boundary inside the loop?
                // auto f = f_increment * step; // in the boundary conditions
                // bc.get("Neumann"), w * f * nv(G)) ????

        // Neumann along one edge???? Neumann on a single degree of freedom?

        for (index_t stag = 0; stag!=1; stag++)
        {
            // if (stag > 0)

            // Before the NR disp --> displacement step???????
            for (index_t it_eq = 0; it_eq!=maxIt; it_eq++) 
            {
                
                // Unew = Uold;

                // Pending: Dirichlet boundary conditions!!!!
                A_u.clearMatrix(); //?
                A_u.clearRhs();

                auto C = A_u.getCoeff(SvK_C);

                auto delta_EE = 0.5*(phys_jacobian_space.cwisetr() + phys_jacobian_space);
                auto delta_EE_voigt = voigt(delta_EE);

                auto bilinear_form = delta_EE_voigt*reshape(C,sz,sz)*delta_EE_voigt.tr();

                gsVector<> pt(mp.parDim());
                if (mp.parDim() == 2) 
                    pt.col(0)<<0.,0.5;
                else if (mp.parDim() == 3)
                    pt.col(0)<<0.,0.5,1.;
                else 
                    GISMO_ERROR("Dimension not supported");

                //gsDebugVar(ev_u.eval(reshape(C,sz,sz),pt,0)); // degraded C
                // gsInfo<<"hola0\n";
                // //gsDebugVar(ev_u.eval(reshape(C,sz,sz),pt,0)); // degraded C
                // gsInfo<<"hola\n";

                // Compute stiffness matrix
                A_u.assemble(bilinear_form * meas(G)); // stiffness matrix
                K_u = A_u.matrix();
                solver.compute(K_u);
                // gsDebugVar(K_u.toDense());

                // Boundary term (Body forces? add missing terms)
                auto f = 0; //change this
                A_u.assembleBdr(bc.get("Neumann"), w * f * nv(G)); // missing the f function!
                
                // Compute RHS
                Q_u = A_u.rhs();

                // Update displacements
                deltaU = solver.solve(-(K_u * Unew - Q_u)); // update vector
                Unew += deltaU; // update the solution vector Unew
                unew.extract(mp_unew);
                mp_def.patch(0).coefs() = mp.patch(0).coefs() + mp_unew.patch(0).coefs(); // update deformed configuration
                // mp_def.patch(0).coefs() += mp_deltaunew.patch(0).coefs();

                A_u.clearMatrix(); 
                A_u.clearRhs();
                A_u.assemble(bilinear_form * meas(G)); // stiffness matrix
                K_u = A_u.matrix();
                Q_u = A_u.rhs();

                // Compute residual norm
                real_t res_u =  (K_u * Unew - Q_u).norm(); ///?????????????

                // ========== Staggered convergence check ==========
                if (stag == 0 && it_eq == 0) 
                {
                    res_stag_norm0 = res_u;
                    gsInfo<<"\tStaggered  iter   "<<stag<<": res_stag = "<<res_stag_norm/res_stag_norm0<<"\n";
                }
                else if (stag>0)
                {
                    res_stag_norm = res_u;
                    gsInfo<<"\tStaggered  iter   "<<stag<<": res_stag = "<<res_stag_norm/res_stag_norm0<<"\n";
                }         
                if (stag>0 && res_u < tol_stag)
                {
                    gsInfo<<"\tStaggered     converged in "<<stag+1<<" iterations\n";
                    converged = true;
                    break;
                }
                else if (stag==maxStag-1)
                {
                    gsInfo<<"\tStaggered          did not converge!\n";
                    converged = false;
                    break;
                }
                // =================================================

                // ========= Elasticity convergence check =========
                if (it_eq == 0) res_u_norm0 = res_u;
                else         res_u_norm = res_u;
                gsInfo<<"\t\tNR EQ iter   "<<it_eq<<": res_u = "<<res_u<<"\n";
                gsInfo<<"\t\tNorm     "<<deltaU.norm()/Unew.norm()<<"\n";


                //if (it_eq>=0 && res_u_norm/res_u_norm0 < tol_NR)
                if (it_eq>0 && res_u_norm/res_u_norm0 < tol_NR)
                {
                    gsInfo<<"\t\tEquilibrium converged in "<<it_eq<<" iterations\n";
                    converged = true; // not using it atm?
                    break;
                }
                else if (it_eq==maxIt-1)
                {
                    gsInfo<<"\t\tEquilibrium did not converge!\n";
                    converged = false;
                    break;
                }
                // =================================================


            } // NR equilibrium equation

            //  Elastic problem -- Plot solution
            if (plot) 
            {
            gsInfo << "Plotting in Paraview ... ";

            gsParaviewCollection collection("ParaviewOutput_trial/disp", &ev_u);
            collection.options().setSwitch("plotElements", true);
            collection.options().setInt("plotElements.resolution", sample_rate);
            collection.options().setInt("numPoints", 10000); 
            collection.newTimeStep(&mp);
            collection.addField(unew, "numerical solution");
            // if (compute_error) {
            // collection.addField(reference_solution - u, "error");
            // }
            collection.saveTimeStep();
            collection.save();

            gsFileManager::open("ParaviewOutput_trial/disp.pvd");
            gsInfo << "Done" << std::endl;
            }

            for (index_t it_pf = 0; it_pf!=3; it_pf++)
            {                
                A_d.clearMatrix();
                A_d.clearRhs();

                real_t penalty = 0.0;
                // if alpha - alpha_old is negative, change penalty parameter to take into account integral in the formulation!
                // if (dnew < dold) // solution from the previous iteration is greater than the new solution
                //     penalty = penalty_irrev;
                // Check ppartval()!!!! muy importante

                auto strain_energy_pos = A_d.getCoeff(SvK_Psi);
                // ev_d.eval(strain_energy_pos,pt,0);

                // Stiffness matrix
                auto K_d_expr = (b*b.tr()) * (2.0 * strain_energy_pos.val() + G_c*2.0/(c_w * ell) + penalty) +   
                                (G_c * 2.0 * ell/c_w) * igrad(b,G)* igrad(b,G).tr();

                A_d.assemble(K_d_expr * meas(G)); // assemble stiffness matrix
                K_d = A_d.matrix();                                 
                solver.compute(K_d);

                // gsDebugVar(K_d.toDense().maxCoeff());
                // gsDebugVar(K_d.toDense().minCoeff());

                K_d_old = K_d;

                // gsDebugVar(K_d.toDense()); //?????? mirar esto con  menos grados d libertad

                // RHS
                auto f_d_expr = b * (2.0 * strain_energy_pos.val());// + penalty * dold); // revisar si es RHS!
                A_d.assemble(f_d_expr * meas(G)); // assemble RHS
                Q_d = A_d.rhs();

                // Compute solution
                deltaD = solver.solve(-(K_d * Dnew - Q_d)); // update damage vector
                Dnew += deltaD;
                dnew.extract(mp_dnew); // extract the solution in mp_dnew (input of gsLinearDegradedMaterial)
                //mp_def.patch(0).coefs() = mp_dnew.patch(0).coefs() + mp_dnew.patch(0).coefs(); // update deformed configuration
                //mp.patch(0).coefs() += mp_dnew.patch(0).coefs(); // update deformed config

                // Compute the residual norm (with the new solution)
                A_d.clearMatrix(); 
                A_d.clearRhs();
                A_d.assemble(K_d_expr * meas(G)); // stiffness matrix
                K_d = A_d.matrix();
                Q_d = A_d.rhs();

                // gsDebugVar(K_d.toDense().maxCoeff());
                // gsDebugVar(K_d.toDense().minCoeff());

                gsDebugVar((K_d.toDense()-K_d_old.toDense()).norm());

                // Compute rdesidual norm
                real_t res_d =  (K_d * Dnew - Q_d).norm(); // esto no se si esta bien???
                
                // ================ Convergence check ================
                if (it_pf == 0) res_d_norm0 = res_d;
                else         res_d_norm = res_d;
                // gsInfo<<"\t\tNR PF iter   "<<it_pf<<": res_d/res_d0 = "<<res_d_norm/res_d_norm0<<"\n";
                gsInfo<<"\t\tNR PF iter   "<<it_pf<<": res_d = "<<res_d<<"\n";

                if (it_pf>0 && res_d_norm/res_d_norm0 < tol_NR)
                {
                    gsInfo<<"\t\tPhase-field   converged in "<<it_pf<<" iterations\n";
                    converged = true;
                    break;
                }
                else if (it_pf==maxIt-1)
                {
                    gsInfo<<"\t\tPhase-field   did not converge!\n";
                    converged = false;
                    break;
                }
                // ====================================================

            } // NR phase field equation

            //  Damage problem -- Plot solution
            if (plot) 
            {
                gsInfo << "Plotting in Paraview ... ";

                gsParaviewCollection collection("ParaviewOutput_dmg/dmg", &ev_d);
                collection.options().setSwitch("plotElements", true);
                collection.options().setInt("plotElements.resolution", sample_rate);
                collection.options().setInt("numPoints", 10000); 
                collection.options().setInt("precision", 7); // 1e-7
                // gsVector<> pt(2);
                // pt << 0, 0.5;
                // gsDebugVar(ev_d.eval(dnew,pt,0)); // gives me something slightly bigger than 1....
                collection.newTimeStep(&mp);
                collection.addField(dnew, "damage");
                // if (compute_error) {
                // collection.addField(reference_solution - u, "error");
                // }
                collection.saveTimeStep();
                collection.save();

                gsFileManager::open("ParaviewOutput_dmg/dmg.pvd");
                gsInfo << "Done" << std::endl;
            }

        } // staggered iteration

    } // step iteration

    // if (plot)
    // {
    //     collection.save();
    // }
    // else
    //     gsInfo << "Done. No output created, re-run with --plot to get a ParaView "
    //               "file containing the solution.\n";

    return EXIT_SUCCESS;

} // end main
    
