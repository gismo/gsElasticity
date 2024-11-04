/** @file fracture_elasticity_example.cpp

    @brief Tutorial on how to use expression assembler to solve the Cahn-Hilliard equation

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): 
    
    To run the script in 2D:
    ./bin/fracture_elasticity_example -f linear_elasticity_example_singlepatch_2d.xml -r 8 --plot
    ./bin/fracture_elasticity_example -f linear_elasticity_example_singlepatch_2d.xml -r 7 --plot
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

    gsBoundaryConditions<> bc, bc_d;
    file_data.getId(2, bc);  // id=2: boundary conditions
    // bc.setGeoMap(mp);
    // gsInfo << "Boundary conditions:\n" << bc << "\n";

    //// ===================== MUY IMPORTANTE =====================
    // gsConstantFunction<> applied_displacement(bcOpt.getReal("stepsize"),mp.geoDim());
    // bc.addCondition(bcOpt.getInt("side"),condition_type::dirichlet,&applied_displacement,0,false,bcOpt.getInt("dir"));
    // ============================================================

    gsFunctionExpr<> reference_solution_expr;
    file_data.getId(3, reference_solution_expr);  // id=3: reference solution

    gsOptionList assembler_options = gsAssembler<>::defaultOptions();
    if (file_data.hasId(4)) {
        file_data.getId(4, assembler_options);  // id=4: assembler options
    }

    //// ===================== MUY IMPORTANTE =====================
    gsOptionList applied_disp_options;
    if (file_data.hasId(5)) {
        file_data.getId(5, applied_disp_options);  // id=4: assembler options
    }

    real_t delta_disp = applied_disp_options.getReal("maxDisp")/applied_disp_options.getInt("steps");
    //gsDebugVar(applied_disp_options.getReal("maxDisp"));
    gsConstantFunction<> applied_displacement(delta_disp,mp.geoDim()); // apply the initial value of the boundary condition (it will be updated with the step)!

    ////////////////// no bc.addCondition(applied_disp_options.getInt("side"),condition_type::dirichlet,&applied_displacement,0,false,applied_disp_options.getInt("direction"));
    bc.addCondition(0,applied_disp_options.getInt("side"),condition_type::dirichlet,&applied_displacement,0,false,applied_disp_options.getInt("component"));
    bc.setGeoMap(mp);
    
    bc_d.addCondition(0,boundary::west,condition_type::clamped,0,0,false,0);
    bc_d.addCondition(0,boundary::east,condition_type::clamped,0,0,false,0);
    bc_d.addCondition(0,boundary::south,condition_type::clamped,0,0,false,0);
    bc_d.addCondition(0,boundary::north,condition_type::clamped,0,0,false,0);
    bc_d.setGeoMap(mp);
    gsInfo << "Boundary conditions disp:\n" << bc << "\n";
    gsInfo << "Boundary conditions dmg:\n" << bc_d << "\n";

    // ============================================================

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
    auto unew_d = A_d.getCoeff(mp_unew); // ??????


    // // Add clamped bc to the phase field (is it right?)
    // for ( gsMultiPatch<>::const_biterator
    // bit = mp.bBegin(); bit != mp.bEnd(); ++bit)
    // {
    //     bc_d.addCondition( *bit, condition_type::clamped,0);
    // }

    // Setup the spaces (compute Dirichlet BCs)
    w.setup(bc, dirichlet::l2Projection, 0); // computes the values on CPs of dirichlet BC   
    b.setup(bc_d, dirichlet::l2Projection, 0); // no damage boundary conditions 
    //b.setup(); // computes the values on CPs of dirichlet BC   

    Dold.setZero(A_d.numDofs(),1);
    Uold.setZero(A_u.numDofs(),1);
    
    Dnew = Dold; 
    Unew = Uold;

    
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
    // Material properties for SEN tensile test (Greco et al.)
    real_t G_c = 2.7e-3; // in kN/mm
    real_t ell = 0.01; // in mm
    //real_t ell = 0.0000015; // in mm
    real_t E_modulus = 210; // Young's modulus in kN/mm2
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

    gsFunctionExpr<> fun("if ((x> -0.05) and (x < 0.50) and (y > 0.5-0.002) and (y < 0.5+0.002)) { 1.0 } else { 0. };",mp.geoDim()); //3 dimension
    //gsFunctionExpr<> fun("if ((x> -0.05) and (x < 0.50) and (y = 0.5)) { 1. } else { 0. };",mp.geoDim()); //3 dimension
    gsMatrix<> anchors = dbasis.basis(0).anchors();
    gsMatrix<> vals = fun.eval(anchors);
    vals.transposeInPlace();
    //gsDebugVar(vals);
    gsGeometry<>::uPtr geom_ptr = dbasis.basis(0).makeGeometry(give(vals));   
    gsWriteParaview(mp,*geom_ptr,"fun2",1e7);

    //L2 INTERPOLATION! (i get negative values how does luigi do it?)
    // gsMatrix<> coefs_interp;
    // gsL2Projection<real_t>::projectFunction(dbasis.basis(0),fun,mp,coefs_interp);
    // gsGeometry<>::uPtr geom_l2 = dbasis.basis(0).makeGeometry(give(coefs_interp));   
    // gsWriteParaview(mp,*geom_l2,"fun2_l2",1e5);

    mp_dnew.addPatch(*geom_ptr); // add a patch to mp_dnew
    dnew.insert(mp_dnew.patch(0),0); // CHECK! update the coefficients in the solution!
    dold.insert(mp_dnew.patch(0),0); // CHECK! update the coefficients in the solution!
    //gsDebugVar(dnew.maxCoeff());
    // evaluation to check that dnew equal mp_dnew (solution updated to initial condition)
    //gsDebugVar(Dold);
    // New solution (check!)
    // Dnew = Dold; 
    // Unew = Uold;
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
    // this is to be put in the xml file?
    real_t tol_NR = 1e-7;
    real_t tol_stag = 1e-6;
    index_t maxStag = 10; // Staggered iterations
    // index_t maxSteps = 10; // Time steps
    index_t maxSteps = applied_disp_options.askInt("steps",10);
    index_t maxIt = 50; // NR iterations
    real_t tol_irrev = 1e-3;
    real_t penalty_irrev = G_c/ell*(1/pow(tol_irrev,2)-1); // Penalty formulation (T. Gerasimov et al., 2019) (2.7e5)
    //real_t penalty_irrev = 0.0; // Penalty formulation (T. Gerasimov et al., 2019) (2.7e5)
    //gsDebugVar(penalty_irrev);
    //real_t penalty_irrev = 1e6; // Penalty formulation (T. Gerasimov et al., 2019)

    // // gsParaviewCollection collection("ParaviewOutput/solution", &ev);
    // gsParaviewCollection collection("ParaviewOutput_fracture/solution", &ev_u); // Which evaluator?
    // collection.options().setSwitch("plotElements", true);
    // collection.options().setInt("plotElements.resolution", 4);
    // collection.options().setInt("numPoints", 10000); 

    real_t time = 0;
    bool converged = false;

    // Initialize the system(s)
    A_u.initSystem(); 
    A_d.initSystem();

    gsSparseMatrix<> K_u, K_d, K_d_old, K_d_jaco;
    gsMatrix<> Q_u, Q_d;
    gsMatrix<> coefs_strain;
    // Voigt dimension
    short_t sz = (mp.parDim()+1)*mp.parDim()/2; 

    // These work by reference (updating mp_def already updates C/Psi...)
    gsMaterialEval<real_t, gsMaterialOutput::C, true> SvK_C(SvK.get(), &mp_def); // change name of mp?
    gsMaterialEval<real_t, gsMaterialOutput::Psi, true> SvK_Psi(SvK.get(), &mp_def);
    gsMaterialEval<real_t, gsMaterialOutput::E, true> SvK_E(SvK.get(), &mp_def);

    real_t lambda = 1/maxSteps;

    for (index_t step = 0; step!=1; step++)
    //for (index_t step = 0; step!=1; step++)
    {
        // gsInfo<<"Time step "<<step<<"/"<<maxSteps<<", iteration "<<dt_it<<": dt = "<<dt<<", [t_start,t_end] = ["<<time<<" , "<<time+dt<<"]"<<"\n";
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
        
        // ================ Update the applied displacement ================ 
        real_t appl_displ_step = delta_disp*(step+1);
        applied_displacement.setValue(appl_displ_step,mp.parDim()); // udpate the applied displacement!
        w.setup(bc, dirichlet::l2Projection, 0); //  // compute again after updating the values of the applied displacement 
        
        gsInfo<<"================================================================================\n";
        gsInfo<<" Time step "<<step<<"/"<<maxSteps<<", du: "<<appl_displ_step<<"\n";
        gsInfo<<"================================================================================\n\n";
        // =================================================================

        for (index_t stag = 0; stag!=1; stag++)
        {
            // if (stag > 0)

            // Before the NR disp --> displacement step???????
            for (index_t it_eq = 0; it_eq!=maxIt; it_eq++) 
            {
                A_u.clearMatrix();
                A_u.clearRhs();

                auto C = A_u.getCoeff(SvK_C);
                auto E_voigt = A_u.getCoeff(SvK_E); 

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
                // gsDebugVar(SvK_E.piece(0).eval(pt));
                
                // Compute stiffness matrix
                A_u.assemble(bilinear_form * meas(G)); // stiffness matrix
                K_u = A_u.matrix();
                solver.compute(K_u);
                // gsDebugVar(K_u.toDense());

                // Boundary term (Body forces? add missing terms)
                auto f = 0; //change this
                A_u.assembleBdr(bc.get("Neumann"), w * f * nv(G)); // missing the f function!
                
                auto res_el = delta_EE_voigt * (reshape(C,sz,sz) * E_voigt);
                // auto K_nl = ;
                //A_u.assemble(res_el * meas(G)); // assemble RHS
                //A_u.assembleJacobian(res_el * meas(G),unew); // assemble RHS

                //A_u.assemble(bilinear_form.dot(dnew) * meas(G));
                // Compute RHS
                //Q_u = K_u * Unew;                 // Perform matrix multiplication
                Q_u = A_u.rhs();

                // Update displacements
                deltaU = solver.solve(-(K_u * Unew - Q_u)); 
                real_t res_u =  (K_u * Unew - Q_u).norm(); ///?????????????

                Unew += deltaU; // update the solution vector Unew
                unew.extract(mp_unew);
                mp_def.patch(0).coefs() = mp.patch(0).coefs() + mp_unew.patch(0).coefs(); // update deformed configuration
                // mp_def.patch(0).coefs() += mp_deltaunew.patch(0).coefs();

                // A_u.clearMatrix(); 
                // A_u.clearRhs();
                // A_u.assemble(bilinear_form * meas(G)); // stiffness matrix
                // //A_u.assemble(res_el * meas(G)); // assemble RHS
                // K_u = A_u.matrix();
                // Q_u = A_u.rhs();

                // Compute residual norm
                /////// real_t res_u =  (K_u * Unew - Q_u).norm(); ///?????????????
                
                //real_t res_u =  (Q_u).norm(); ///?????????????

                // ========== Staggered convergence check ==========
                if (stag == 0 && it_eq == 0) 
                {
                    res_stag_norm0 = res_u;
                    gsInfo<<" Staggered  iter   "<<stag<<": res_stag = "<<res_stag_norm/res_stag_norm0<<"\n";
                    gsInfo<<"--------------------------------------------------------------------------------\n\n";
                }
                else if (stag>0)
                {
                    res_stag_norm = res_u;
                    gsInfo<<" Staggered  iter   "<<stag<<": res_stag = "<<res_stag_norm/res_stag_norm0<<"\n";
                    gsInfo<<"--------------------------------------------------------------------------------\n\n";
                }         
                if (stag>0 && res_u < tol_stag)
                {
                    gsInfo<<" Staggered converged in "<<stag+1<<" iterations\n";
                    converged = true;
                    break;
                }
                else if (stag==maxStag-1)
                {
                    gsInfo<<" Staggered did not converge!\n";
                    converged = false;
                    break;
                }
                // =================================================

                // ========= Elasticity convergence check =========
                gsInfo<<"\tNR EQ iter     "<<it_eq<<" : R_u_norm = "<<res_u<<", deltaU.norm() = "<<deltaU.norm()<<"\n";


                if (it_eq>0 && res_u < tol_NR)
                {
                    gsInfo<<"\tEquilibrium converged in "<<it_eq<<" iterations\n\n";
                    converged = true; // not using it atm?
                    break;
                }
                else if (it_eq==maxIt-1)
                {
                    gsInfo<<"\tEquilibrium did not converge!\n\n";
                    converged = false;
                    break;
                }
                // =================================================


            } // NR equilibrium equation

            //  Elastic problem -- Plot solution
            if (plot) 
            {
                gsInfo << "Plotting in Paraview ... ";

                gsParaviewCollection collection1("ParaviewOutput_fracture/disp", &ev_u);
                collection1.options().setSwitch("plotElements", true);
                collection1.options().setInt("plotElements.resolution", sample_rate);
                collection1.options().setInt("numPoints", 10000); 
                collection1.newTimeStep(&mp);
                collection1.addField(unew, "numerical solution");
                // if (compute_error) {
                // collection.addField(reference_solution - u, "error");
                // }
                collection1.saveTimeStep();
                collection1.save();

                gsFileManager::open("ParaviewOutput_fracture/disp.pvd");
                gsInfo << "Done" << std::endl;
            }

            for (index_t it_pf = 0; it_pf!=13; it_pf++)
            {                
                A_d.clearMatrix();
                A_d.clearRhs();

                auto strain_energy_pos = A_d.getCoeff(SvK_Psi); // I want to plot the energy!!
                // gsWriteParaview(strain_energy_pos,G,"ParaviewOutput_fracture/energy"+util::to_string(it_pf));
                auto strain = A_d.getCoeff(SvK_E); // I want to plot the energy!!
                //ev_d.writeParaview(strain,G,"ParaviewOutput_fracture/strain"+util::to_string(it_pf));
                
                // STRAIN HISTORY?
                // auto strain_old = strain_energy_pos;
                // auto H_energy = ternary((strain_energy_pos-strain_old).val(), strain_energy_pos, strain_old); // no se si tiene sentido....
                
                gsVector<> pt(mp.parDim());
                gsVector<> pt1(mp.parDim());

                if (mp.parDim() == 2) 
                    pt.col(0)<<0.,0.5;
                else if (mp.parDim() == 3)
                    pt.col(0)<<0.,0.5,1.;
                else 
                    GISMO_ERROR("Dimension not supported");

                pt1.col(0)<<0.7,0.5;

                // We add a ternary of dnew-dold; if true, we do not apply the penalty!
                // This penalty term is added to the residual with npart(), which takes the negative part of the expression dnew-dold

                auto penalty_constant = ternary((dnew-dold).val()+1e-5, 0.0*dnew.val(), penalty_irrev + 0.0*dnew.val()); // no se si tiene sentido....
                auto penalty_term_residual = ternary((dnew-dold).val()+1e-5, 0.0*dnew.val(), penalty_irrev*(dnew-dold).val()); 

                // Residual (Newton-Raphson)
                auto res_d_expr = b * (2.0*(-1.0+dnew.val()) * strain_energy_pos.val() + (G_c*2.0/(ell*c_w)) * dnew.val() + penalty_term_residual) 
                                   + G_c*ell*2.0/(c_w)*igrad(b,G)*igrad(dnew,G).tr();
                // auto res_d_expr = b * (2.0*(-1.0+dnew.val()) * strain_energy_pos.val() + (G_c*2.0/(ell*c_w)) * dnew.val() ) 
                //                      + G_c*ell*2.0/(c_w)*igrad(b,G)*igrad(dnew,G).tr() ;

                // Tangent stiffness matrix (Newton-Raphson)
                auto K_d_expr = (b*b.tr()) * (2.0*G_c/(ell*c_w) + 2.0*strain_energy_pos.val() + penalty_constant) 
                                + 2.0*G_c*ell/c_w*igrad(b,G)*igrad(b,G).tr();

                A_d.assemble(K_d_expr * meas(G), res_d_expr * meas(G));
                // A_d.assembleJacobian(res_d_expr * meas(G), dnew); // assemble stiffness matrix  ------ this is not working! it does not give the same as with my derivation!
                // A_d.assemble(res_d_expr * meas(G));
                
                K_d = A_d.matrix();
                Q_d = A_d.rhs();
                solver.compute(K_d); // factorization of K

                // A_d.initSystem();
                // A_d.assembleJacobian(res_d_expr * meas(G), dnew);
                // K_d_jaco = A_d.matrix();
                // gsSparseMatrix<> jaco = A_d.matrix();
                // gsDebugVar((jaco-K_d).norm());

                gsInfo<<"==========================================\n";
                gsDebugVar(ev_d.eval(dnew,pt,0));
                gsDebugVar(ev_d.eval(dnew,pt1,0));

                // Compute solution
                deltaD = solver.solve(-Q_d); // update damage vector
                Dnew += deltaD;

                //gsWriteParaview(dnew,"ParaviewOutput_fracture/damage_evol"); //??? no funciona que cojones?
                // gsBSpline<> solution2(basis_plot,Dnew);
                // gsWriteParaview(solution2,"ParaviewOutput_fracture/solution2"); // interpolated function mySinus

                // ev_d.writeParaview(dnew,G,"hola"); // works
                // ev_d.writeParaview(strain,G,"ParaviewOutput_fracture/strain"+util::to_string(it_pf)); // does not work!
                // gsField<> solField(mp,dnew);
                // gsWriteParaview(solField,"solution_step",1000,true);
                // gsDebugVar(Dnew.norm());
                ev_d.writeParaview(dnew,G,"ParaviewOutput_fracture/damage_"+util::to_string(it_pf)); // does not work!


                gsInfo<<"==========================================\n";


                //Dold = Dnew; // for checking the penalty

                dnew.extract(mp_dnew); // extract the solution in mp_dnew (input of gsLinearDegradedMaterial)                
                
                // gsDebugVar(ev_d.eval(dnew,pt,0));
                // gsDebugVar(ev_d.eval(dnew,pt1,0));
                
                real_t res_d =  Q_d.norm(); // esto no se si esta bien???

                // ================ Convergence check ================
                gsInfo<<"\tNR PF iter     "<<it_pf<<": R_d_norm = "<<res_d<<", deltaD.norm() = "<<deltaD.norm()<<"\n";

                if (it_pf>0 && res_d < tol_NR)
                {
                    gsInfo<<"\tPhase-field converged in "<<it_pf<<" iterations\n\n";
                    converged = true;
                    break;
                }
                else if (it_pf==maxIt-1)
                {
                    gsInfo<<"\tPhase-field   did not converge!\n";
                    converged = false;
                    break;
                }
                // ====================================================



            } // NR phase field equation
            
            // Dold = Dnew; // for checking the penalty

            //  Damage problem -- Plot solution
            if (plot) 
            {
                gsInfo << "Plotting in Paraview ... ";

                gsParaviewCollection collection2("ParaviewOutput_fracture/dmg", &ev_d);
                collection2.options().setSwitch("plotElements", true);
                collection2.options().setInt("plotElements.resolution", sample_rate);
                collection2.options().setInt("numPoints", 10000); 
                collection2.options().setInt("precision", 7); // 1e-7
                collection2.newTimeStep(&mp);
                collection2.addField(dnew, "damage");
                // if (compute_error) {
                // collection.addField(reference_solution - u, "error");
                // }
                collection2.saveTimeStep();
                collection2.save();

                gsFileManager::open("ParaviewOutput_fracture/dmg.pvd");
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
    
