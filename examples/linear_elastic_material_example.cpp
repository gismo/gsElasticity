/** @file gsMaterialMatrix_test.cpp

    @brief Simple example for material matrix evaluations

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): H.M. Verhelst
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

int main (int argc, char** argv)
{
    bool plot = false;
    index_t n_h_refinements = 0;

    //  ====== NO SE SI HACE FALTA ESTO ======
    int material = 0;
    int impl = 1;
    int Compressibility = 0;
    // =======================================

    std::string file_name("pde/linear_elasticity_example_singlepatch.xml");
    gsCmdLine cmd("Tutorial on solving a Linear Elasticity problem with the gsMaterialEval class");
    // cmd.addInt("e", "degreeElevation",
    //             "Number of degree elevation steps to perform before solving (0: "
    //             "equalize degree in all directions)",
    //             n_deg_elevations);
    cmd.addInt("r", "uniformRefine", "Number of Uniform h-refinement loops",
                n_h_refinements);
    // cmd.addString("f", "file", "Input XML file", file_name);
    // cmd.addSwitch("only_last",
    //                 "Solve solely for the only_last level of h-refinement",
    //                 only_last);
    cmd.addSwitch("plot",
                    "Create a ParaView visualization file with the solution", plot);
    // cmd.addInt("q", "sample-rate", "Samples per spline in paraview export",
    //             sample_rate);
    // cmd.addSwitch("compute-error",
    //                 "Evaluate the error with respect to the analytical solution "
    //                 "(evaluation with default options and default file required)",
    //                 compute_error);
    try { cmd.getValues(argc,argv); } catch (int rv) { return rv; }

    gsFileData<> file_data(file_name);
    gsInfo << "Loaded file " << file_data.lastPath() << "\n";

    // Create geometry
    gsMultiPatch<> mp, mp_def;
    mp.addPatch( gsNurbsCreator<>::BSplineCube(1) ); // degree
    mp.addAutoBoundaries();

    // h-refine each basis
    for (int r = 0; r < n_h_refinements; ++r) {
        mp.uniformRefine();
    }

    mp_def = mp;
    //gsDebugVar(mp.patch(0).coefs());
    mp_def.patch(0).coefs().col(0) *= 2; // deformation *2 in the x direction (col(0))

    // ===== Source term =====
    gsFunctionExpr<> f_expr;
    file_data.getId(1, f_expr);  // id=1: source function
    gsInfo << "Source function " << f_expr << "\n";

    // ===== Boundary conditions =====
    gsBoundaryConditions<> bc;
    gsVector<> displ_bc (2); //could be a scalar? check!
    displ_bc<< 0,1; // 0 in x 1(quantity) in y ????
    gsConstantFunction<> displ_bc_fun (displ_bc,2);
    // ERROR: Expecting scalar function
    // bc.addCondition(boundary::side::north,condition_type::dirichlet,&displ_bc_fun,0,false,1); //last componnent is the component x(0), y(1), z(2)
    // bc.addCondition(boundary::side::south,condition_type::dirichlet,0,0,false,-1); //last componnent -1 fixes all displacements
    // bc.addCondition(boundary::side::east,condition_type::dirichlet,0,0,false,1); //last componnent -1 fixes all displacements

    // BC for phase field (fracture)
    // bc.addCondition(boundary::side::north,condition_type::dirichlet,&displ_bc_fun,1); // 1 because of the ID of the space (phase field is 1)
    bc.setGeoMap(mp);
    // ===============================   

    // ===== Material properties =====
    real_t thickness = 1.0;
    real_t E_modulus = 210e9;
    real_t PoissonRatio = 0.3;
    auto lambda = E_modulus * PoissonRatio / ( ( 1. + PoissonRatio ) * ( 1. - 2. * PoissonRatio ) );
    auto mu     = E_modulus / ( 2. * ( 1. + PoissonRatio ) );
    gsFunctionExpr<> t(std::to_string(thickness), 3);
    gsFunctionExpr<> E(std::to_string(E_modulus),3);
    gsFunctionExpr<> nu(std::to_string(PoissonRatio),3);

    // Linear material class
    // New structure: introduce a flag for Voigt notation (true) or tensor notation (false)
    // gsMaterialEval<real_t, gsMaterialOutput::S, true> SvK_S_deg(&SvK_deg, &mp_def);
    gsLinearMaterial<3,real_t> SvK(E_modulus, PoissonRatio, mp, mp_def); 
    gsMaterialEval<real_t, gsMaterialOutput::E, true> SvK_E(&SvK, &mp_def);
    gsMaterialEval<real_t, gsMaterialOutput::S, true> SvK_S(&SvK, &mp_def);
    gsMaterialEval<real_t, gsMaterialOutput::C, true> SvK_C(&SvK, &mp_def);
    
    // Evaluation of the variables at point pt
    gsVector<> pt(3);
    pt.col(0)<<0.125,0.375,0.5;

    gsMatrix<> Sres;

    // // strain
    // gsInfo<<"Strain\n";
    // SvK_E.piece(0).eval_into(pt,Sres);
    // gsDebugVar(Sres);

    // // stress
    // gsInfo<<"Undegraded stress\n";
    // SvK_S.piece(0).eval_into(pt,Sres);
    // gsDebugVar(Sres);

    // //material matrix 
    // gsInfo<<"Undegraded material matrix\n";
    // SvK_C.piece(0).eval_into(pt,Sres);
    // gsDebugVar(Sres);
    
    // ============================================================
    // Expression assembler with new material class
    gsExprAssembler<> A(1,1);
    gsExprEvaluator<> ev(A);

    typedef gsExprAssembler<>::geometryMap geometryMap;
    typedef gsExprAssembler<>::variable    variable;
    typedef gsExprAssembler<>::space       space;
    typedef gsExprAssembler<>::solution    solution;
    
    gsMatrix<> disp;
    //gsMatrix<> disp =  mp_def.patch(0).coefs()- mp.patch(0).coefs(); // displacement vector

    gsMultiBasis<> dbasis(mp,true);
    
    space w = A.getSpace(dbasis,3); //should be dimension-independent

    A.setIntegrationElements(dbasis); // necesito escribir la basis!

    geometryMap G = A.getMap(mp); // al mp o al mp_def? domain dimension 3

    solution u = A.getSolution(w, disp); // 3 rows one column
   
    // Set the source term
    auto f = A.getCoeff(f_expr, G);

    disp = disp.reshape(disp.cols()*disp.rows(),1); // to ensure Assert `_u.coefs().rows()==_u.mapper().freeSize()`

    w.setup(); //setup allocates the degree-of-freedom mapper
    
    A.initSystem();
    
    // Apply dirichlet boundary conditions
    w.setup(bc, dirichlet::l2Projection, 0);

    gsInfo<<"Hola1\n";
    gsMatrix<> I = gsMatrix<>::Identity(3,3); // es 3x3?

    auto phys_jacobian = ijac(u, G);
    auto phys_jacobian_space = ijac(w, G);

    auto delta_EE = 0.5*(phys_jacobian_space.cwisetr() + phys_jacobian_space);
    auto delta_EE_voigt = voigt(delta_EE);
    //auto delta_EE_NL = 0.5*(phys_jacobian_space.cwisetr() + phys_jacobian_space + phys_jacobian*phys_jacobian_space.cwisetr() + phys_jacobian_space*phys_jacobian.cwisetr()); // linearization with geometric non-linearities

    auto C_undgr = A.getCoeff(SvK_C); // undegraded material matrix
    auto bilinear_form = delta_EE_voigt*reshape(C_undgr,6,6)*delta_EE_voigt.tr()*meas(G);

    //auto f = 0; // como lo hago?
    auto linear_form = w * f * meas(G); //
    
    A.clearMatrix();
    A.assemble(bilinear_form);
    //gsDebugVar(A.giveMatrix().toDense());
    A.assemble(linear_form);
    //gsDebugVar(linear_form);
    
    // gsMatrix<> matr_assembly = A.matrix().toDense();
    // gsDebugVar(matr_assembly);

    gsSparseSolver<>::LU solver;
    solver.compute(A.matrix());
    disp = solver.solve(A.rhs()); // would be disp

    // Add Paraview representation
    if (plot) {
        gsInfo << "Plotting in Paraview ... ";

        gsParaviewCollection collection("ParaviewOutput/solution", &ev);
        collection.options().setSwitch("plotElements", true);
        // collection.options().setInt("plotElements.resolution", sample_rate);
        collection.newTimeStep(&mp);
        collection.addField(u, "numerical solution");
        // if (compute_error) {
        // collection.addField(reference_solution - u_solution_expression, "error");
        // }
        collection.saveTimeStep();
        collection.save();

        gsFileManager::open("ParaviewOutput/solution.pvd");
        gsInfo << "Done" << std::endl;
    } else {
        gsInfo << "No output created, re-run with --plot to get a ParaView "
                "file containing the solution.\n";
    } 

    return 0;
}
