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
//#include <gsElasticity/gsMaterialEval_copy.h>
// #include <gsElasticity/gsMaterialBase.h>
#include <gsElasticity/gsLinearMaterial.h>
#include <gsElasticity/gsLinearDegradedMaterial.h>
#include <gsElasticity/gsVisitorElUtils.h>

// #include <gsElasticity/gsNeoHookLogMaterial.h>
// #include <gsElasticity/gsNeoHookQuadMaterial.h>
// #include <gsElasticity/gsMuesliMaterial.h>

// #define WITHEIGEN
// #include <muesli/muesli.h>
// #include <muesli/tensor.h>
// #undef matrix

// #include <gsElasticity/gsWriteParaviewMultiPhysics.h>
// #include <gsElasticity/gsElasticityAssembler.h>
// #include <gsStructuralAnalysis/gsStaticNewton.h>

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
    int material = 0;
    int impl = 1;
    int Compressibility = 0;
    gsCmdLine cmd(".");

    try { cmd.getValues(argc,argv); } catch (int rv) { return rv; }

    // Make geometry
    gsMultiPatch<> mp, mp_def;
    mp.addPatch( gsNurbsCreator<>::BSplineCube(1) ); // degree
    mp.addAutoBoundaries();

    gsDebugVar(mp.basis(0).numElements());
    mp.uniformRefine(); // refines both
    //mp.basis(0).uniformRefine(); // this refines the basis but not the geometry!

    gsDebugVar(mp.basis(0));
    
    mp_def = mp;
    gsDebugVar(mp_def.patch(0).coefs());
    gsDebugVar(mp.patch(0).coefs());
    
    //gsDebugVar(mp.nPatches());

    mp_def.patch(0).coefs().col(0) *= 2; // deformation *2 in the x direction (col(0))
    // mp_def.patch(0).coefs().col(1) += 0.5 * mp_def.patch(0).coefs().col(0);
    // mp_def.patch(0).coefs().col(0) += 0.5 * mp_def.patch(0).coefs().col(1);

    // gsDebugVar(mp.basis(0));
    // I need to construct the solution of the phase field?
     
    real_t thickness = 1.0;
    real_t E_modulus = 210e9;
    real_t PoissonRatio = 0.3;
    real_t PoissonRatio2 = 0.3;
    real_t Ratio = 1;
    gsFunctionExpr<> t(std::to_string(thickness), 3);
    gsFunctionExpr<> E(std::to_string(E_modulus),3);
    gsFunctionExpr<> nu(std::to_string(PoissonRatio),3);
    gsFunctionExpr<> nu2(std::to_string(PoissonRatio2),3);
    gsConstantFunction<> ratio(Ratio,3);

    // ====== Initialize damage ====== !!
    gsMatrix<> tmp = gsMatrix<>::Random(mp_def.patch(0).coefs().rows(),1);
    gsMatrix<> damage = 0.1 * tmp.array() + 0.1; // damage btw 0 and 0.2
    gsGeometry<>::uPtr damage_ = mp_def.basis(0).makeGeometry(give(damage)); // spline (Function)

    // ========================================================================
    // New structure: introduce a flag for Voigt notation (true) or tensor notation (false)
    // gsMaterialEval<real_t, gsMaterialOutput::S, true> SvK_S_deg(&SvK_deg, &mp_def);

    gsLinearMaterial<3,real_t> SvK(E_modulus, PoissonRatio, mp, mp_def); 
    gsMaterialEval<real_t, gsMaterialOutput::E, true> SvK_E(&SvK, &mp_def);
    gsMaterialEval<real_t, gsMaterialOutput::S, true> SvK_S(&SvK, &mp_def);
    gsMaterialEval<real_t, gsMaterialOutput::C, true> SvK_C(&SvK, &mp_def);

    gsLinearDegradedMaterial<3,real_t> SvK_deg(E_modulus, PoissonRatio, mp, mp_def, *damage_);
    gsMaterialEval<real_t, gsMaterialOutput::C, true> SvK_C_deg(&SvK_deg, &mp_def);
    gsMaterialEval<real_t, gsMaterialOutput::S, true> SvK_S_deg(&SvK_deg, &mp_def);
    gsMaterialEval<real_t, gsMaterialOutput::C_pos, true> SvK_C_pos(&SvK_deg, &mp_def);
    // gsNeoHookLogMaterial<real_t> NH(E_modulus, PoissonRatio2, &mp, &mp_def);
    // gsMaterialBase<real_t> base(&mp, &mp_def);

    gsVector<> pt(3);
    pt.col(0)<<0.125,0.375,0.5;

    gsMatrix<> pts(3,6);
    pts.col(0)<<0.125,0.375,0.5;
    pts.col(1)<<0.375,0.125,0.5;
    pts.col(2)<<0.125,0.25,0.5;
    pts.col(3)<<0.25,0.125,0.5;
    pts.col(4)<<0.25,0.25,0.5;
    pts.col(5)<<0.5,0.5,0.5;

    gsMatrix<> Fres, Eres, Sres ;
    
    // // strain
    gsInfo<<"Strain\n";
    SvK_E.piece(0).eval_into(pt,Sres);
    gsDebugVar(Sres);

    // stress
    gsInfo<<"Undegraded stress\n";
    SvK_S.piece(0).eval_into(pt,Sres);
    gsDebugVar(Sres);

    // degraded stress (in Voigt)
    gsInfo<<"Degraded stress\n";
    SvK_S_deg.piece(0).eval_into(pt,Sres);
    gsDebugVar(Sres);

    // //material matrix 
    // gsInfo<<"Undegraded material matrix\n";
    // SvK_C.piece(0).eval_into(pt,Sres);
    // //gsDebugVar(Sres);
    
    // // degraded material matrix
    // gsInfo<<"Degraded material matrix\n";
    // SvK_C_deg.piece(0).eval_into(pt,Sres);
    // //gsDebugVar(Sres);

    // ============================================================
    // Expression assembler check
    gsExprAssembler<> A(1,1);
    gsExprEvaluator<> ev(A);

    typedef gsExprAssembler<>::geometryMap geometryMap;
    typedef gsExprAssembler<>::variable    variable;
    typedef gsExprAssembler<>::space       space;
    typedef gsExprAssembler<>::solution    solution;
    
    gsMatrix<> disp =  mp_def.patch(0).coefs()- mp.patch(0).coefs(); // displacement vector
    gsMatrix<> disp_vec = disp.transpose();

    gsDebugVar(disp);

    gsMultiBasis<> dbasis(mp,true);
    
    space w = A.getSpace(dbasis,3); // 3 dimensional solution field?

    A.setIntegrationElements(dbasis); // necesito escribir la basis!

    geometryMap G = A.getMap(mp); // al mp o al mp_def? domain dimension 3

    solution u = A.getSolution(w, disp); // 3 rows one column

    disp = disp.reshape(disp.rows()*disp.cols(),1); // to ensure Assert `_u.coefs().rows()==_u.mapper().freeSize()`

    w.setup(); //setup allocates the degree-of-freedom mapper
    A.initSystem();

    auto lambda = E_modulus * PoissonRatio / ( ( 1. + PoissonRatio ) * ( 1. - 2. * PoissonRatio ) );
    auto mu     = E_modulus / ( 2. * ( 1. + PoissonRatio ) );
    gsDebugVar(lambda);
    gsDebugVar(mu);

    gsMatrix<> I = gsMatrix<>::Identity(3,3); // es 3x3?

    auto phys_jacobian = ijac(u, G);
    auto phys_jacobian_space = ijac(w, G);

    // Linear strains
    auto EE = 0.5*(phys_jacobian.cwisetr() + phys_jacobian); //?? in matrix??? strain
    // Green lagrange strain
    //auto EE = 0.5*(phys_jacobian.cwisetr() + phys_jacobian + phys_jacobian.cwisetr() * phys_jacobian); //?? in matrix??? strain
    // gsDebugVar(ev.eval(phys_jacobian,pt,0));
    // gsDebugVar(ev.eval(EE,pt,0).cols());
    // gsDebugVar(ev.eval(EE,pt,0).rows());

    auto EE_eval = ev.eval(EE,pt,0); // Strain (small deformations) evaluated at pt (x,y,z)
    gsDebugVar(EE_eval);

    auto SS = lambda * EE.trace().val()*I + 2.0*mu*EE; // stress
    auto SS_eval = ev.eval(SS,pt,0); // Stress evaluated at pt (x,y,z)
    gsDebugVar(SS_eval);

    // ========== Bilinear (linear elasticity script!!!) ==========
    auto bilin_lambda = lambda * idiv(w, G) *
                        idiv(w, G).tr() * meas(G);
    auto bilin_mu_ =
        mu *
        ((phys_jacobian_space.cwisetr() + phys_jacobian_space) % phys_jacobian_space.tr()) *
        meas(G);

    // auto bilin_mu_ =
    //     mu *
    //     ((phys_jacobian_space.cwisetr() + phys_jacobian_space) * 0.5*((phys_jacobian_space.cwisetr() + phys_jacobian_space).tr())) *
    //     meas(G);

    // gsDebugVar((phys_jacobian_space.cwisetr() + phys_jacobian_space).Space);
    // gsDebugVar(((phys_jacobian_space.cwisetr() + phys_jacobian_space) % phys_jacobian_space.tr()).Space);
    // gsDebugVar(((phys_jacobian_space.cwisetr() + phys_jacobian_space) * (0.5*(phys_jacobian_space.cwisetr() + phys_jacobian_space)).tr()).Space);
    // gsDebugVar((phys_jacobian_space.cwisetr()).Space);
    // gsDebugVar((phys_jacobian_space).Space);
    // gsDebugVar((phys_jacobian_space.tr()).Space);

    // gsDebugVar((lambda * idiv(w, G) * idiv(w, G).tr() * meas(G)).Space);
    
    auto bilin_combined = (bilin_lambda + bilin_mu_);
    A.assemble(bilin_combined);
    gsMatrix<> matr_assembly1 = A.matrix().toDense();
    //gsDebugVar(matr_assembly1); 

    // ============================================================

    // ==== Bilinear form !!!!! ====
    auto delta_EE = 0.5*(phys_jacobian_space.cwisetr() + phys_jacobian_space);
    //auto delta_EE_NL = 0.5*(phys_jacobian_space.cwisetr() + phys_jacobian_space + phys_jacobian*phys_jacobian_space.cwisetr() + phys_jacobian_space*phys_jacobian.cwisetr()); // linearization with geometric non-linearities
    
    // gsDebugVar(phys_jacobian_space.cwisetr().Space);
    // gsDebugVar(phys_jacobian_space.Space);
    // gsDebugVar((phys_jacobian*phys_jacobian_space.cwisetr()).Space)
    // gsDebugVar(phys_jacobian_space.cwisetr())

    // delta_EE to Voigt
    auto delta_EE_voigt = voigt(delta_EE);
    auto C_undgr = A.getCoeff(SvK_C);
    auto C_positive = A.getCoeff(SvK_C_pos);
    A.clearMatrix();
    A.assemble(delta_EE_voigt*reshape(C_undgr,6,6)*delta_EE_voigt.tr()*meas(G));
    gsMatrix<> matr_assembly2 = A.matrix().toDense();
    //gsDebugVar(ev.eval(reshape(C_undgr,6,6),pt,0));
    //gsDebugVar(ev.eval(reshape(C_positive,6,6),pt,0));
    // gsDebugVar(matr_assembly2);
    //gsDebugVar(((matr_assembly2-matr_assembly1).array() / matr_assembly1.array()) .abs().maxCoeff()); //relative err
    gsDebugVar((matr_assembly2-matr_assembly1).maxCoeff()); //relative err
    gsDebugVar((matr_assembly2-matr_assembly1).minCoeff()); //relative err

    // =======================

    // Check that voigt is working
    auto EE_voigt = voigt(EE);
    gsDebugVar(ev.eval(EE_voigt,pt,0));

    // NH.eval_vector_into(pt,Sres,0);
    // gsDebugVar(Sres);
    // NH.eval_matrix_into(pt,Sres,0);
    // gsDebugVar(Sres);


    // //////////////////////////////////////////////////////////////////////////////////////////////
    // ////////////////////////////////////////MUESLI////////////////////////////////////////////////
    // //////////////////////////////////////////////////////////////////////////////////////////////

    // muesli::materialProperties properties_svk;
    // properties_svk.insert(std::pair<std::string,real_t>("young",E_modulus));
    // properties_svk.insert(std::pair<std::string,real_t>("poisson",PoissonRatio));
    // muesli::svkMaterial svk("svk", properties_svk);
    // muesli::finiteStrainMP* p_svk = svk.createMaterialPoint();

    // muesli::finiteStrainMaterial * muesliMaterial = new muesli::svkMaterial("svk", properties_svk);

    // gsMuesliMaterial<real_t> muesliSvK(muesliMaterial,&mp);
    // muesliSvK.updateCurrentState(&mp_def);
    // muesliSvK.commitCurrentState();
    // muesliSvK.eval_vector_into(pt,Sres,0);
    // gsDebugVar(Sres);
    // muesliSvK.eval_matrix_into(pt,Sres,0);
    // gsDebugVar(Sres);


    // real_t mu = E_modulus / (2.*(1. + PoissonRatio2));
    // real_t c2 = mu/(Ratio+1);
    // real_t c1 = Ratio*c2;

    // c1 *= 2;
    // c2 *= 2;

    // muesli::materialProperties properties_nh, properties_mr;
    // properties_nh.insert(std::pair<std::string,real_t>("young",E_modulus));
    // properties_nh.insert(std::pair<std::string,real_t>("poisson",PoissonRatio));
    // // properties_nh.insert(std::pair<std::string,real_t>("subtype regularized",true));
    // muesli::neohookeanMaterial nh("nh", properties_nh);
    // muesli::finiteStrainMP* p_nh = nh.createMaterialPoint();

    // properties_mr.insert(std::pair<std::string,real_t>("alpha1",c1));
    // properties_mr.insert(std::pair<std::string,real_t>("alpha2",c2));
    // properties_mr.insert(std::pair<std::string,real_t>("incompressible",true));
    // muesli::mooneyMaterial mr("mr", properties_mr);
    // muesli::finiteStrainMP* p_mr = mr.createMaterialPoint();

    // // std::ofstream stream("data.txt");
    // // mr.print(stream);
    // // // std::cout<<stream;

    // gsMatrix<> F, strain, C;
    // for (index_t k=0; k!= pt.cols(); k++)
    // {
    //     F = Fres.reshapeCol(k,3,3);
    //     strain = Eres.reshapeCol(k,3,3);
    //     gsDebugVar(F);
    //     gsDebugVar(strain);

    //     gsInfo<<"---------------------------------------------------------------------\n";

    //     istensor stress;
    //     itensor4 tangent;
    //     p_svk->updateCurrentState(k,F);
    //     // gsMatrix<real_t,3,3> stress;
    //     p_svk->secondPiolaKirchhoffStress(stress);
    //     gsDebugVar(stress);
    //     p_svk->CauchyStress(stress);
    //     gsDebugVar(stress);
    //     p_svk->convectedTangent(tangent);
    //     gsDebugVar(tangent);

    //     gsDebugVar(p_svk->deformationGradient());


    //     p_svk->commitCurrentState();

    //     gsDebugVar(F.transpose() * F);
    //     gsDebugVar(C);
    //     gsDebugVar(strain);


    //     gsInfo<<"---------------------------------------------------------------------\n";

    //     p_nh->updateCurrentState(k,F);
    //     // gsMatrix<real_t,3,3> stress;
    //     p_nh->secondPiolaKirchhoffStress(stress);
    //     gsDebugVar(stress);
    //     p_nh->CauchyStress(stress);
    //     gsDebugVar(stress);
    //     p_nh->convectedTangent(tangent);
    //     gsDebugVar(tangent);

    //     gsDebugVar(p_nh->deformationGradient());


    //     p_nh->commitCurrentState();

    //     gsInfo<<"---------------------------------------------------------------------\n";

    //     p_mr->updateCurrentState(k,F);
    //     // gsMatrix<real_t,3,3> stress;
    //     p_mr->secondPiolaKirchhoffStress(stress);
    //     gsDebugVar(stress);
    //     p_mr->CauchyStress(stress);
    //     gsDebugVar(stress);

    //     gsDebugVar(p_mr->deformationGradient());


    //     p_mr->commitCurrentState();

    // }

    return 0;
}