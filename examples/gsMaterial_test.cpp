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
// #include <gsElasticity/gsMaterialBase.h>
#include <gsElasticity/gsLinearMaterial.h>
#include <gsElasticity/gsLinearDegradedMaterial.h>

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

    mp_def = mp;
    mp_def.patch(0).coefs().col(0) *= 2; // deformation *2 in the x direction (col(0))


     
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


    gsLinearMaterial<3,real_t> SvK(E_modulus, PoissonRatio2, mp, mp_def); // why twice?
    gsMaterialEval<real_t, gsMaterialOutput::E> SvK_E(&SvK, &mp_def);
    gsMaterialEval<real_t, gsMaterialOutput::S> SvK_S(&SvK, &mp_def);
    gsMaterialEval<real_t, gsMaterialOutput::C> SvK_C(&SvK, &mp_def);

    // gsLinearDegradedMaterial<3,real_t> SvK_deg(E_modulus, PoissonRatio2, mp, mp_def, damage_);
    // gsMaterialEval<real_t, gsMaterialOutput::C> SvK_C(&SvK_deg, &mp_def);

    // gsLinearDegradedMaterial<3,real_t> SvK_deg(E_modulus, PoissonRatio, mp, mp_def, damage);

    // gsMaterialEval<real_t, gsMaterialOutput::C> SvK_C(&SvK, &mp_def);

    // gsNeoHookLogMaterial<real_t> NH(E_modulus, PoissonRatio2, &mp, &mp_def);
    //gsMaterialBase<real_t> base(&mp, &mp_def);

    gsVector<> pt(3);
    pt.col(0)<<0.125,0.375,0.5;

    gsMatrix<> pts(3,6);
    pts.col(0)<<0.125,0.375,0.5;
    pts.col(1)<<0.375,0.125,0.5;
    pts.col(2)<<0.125,0.25,0.5;
    pts.col(3)<<0.25,0.125,0.5;
    pts.col(4)<<0.25,0.25,0.5;
    pts.col(5)<<0.5,0.5,0.5;

    gsMatrix<> Fres, Eres, Sres;
    // strains
    SvK_E.piece(0).eval_into(pt,Sres);
    gsDebugVar(Sres);
    // strsses
    SvK_S.piece(0).eval_into(pt,Sres);
    gsDebugVar(Sres);
    // material matrix
    // SvK_C.piece(0).eval_into(pt,Sres); // problem evaluation of C
    // gsDebugVar(Sres);



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

    // gsDebugVar(disp);
    // gsDebugVar(disp_vec);

    gsMultiBasis<> dbasis(mp,true);
    
    space w = A.getSpace(dbasis,3); // 3 dimensional solution field?

    A.setIntegrationElements(dbasis); // necesito escribir la basis!

    geometryMap G = A.getMap(mp); // al mp o al mp_def? domain dimension 3

    solution u = A.getSolution(w, disp); // 3 rows one column

    w.setup(); //setup allocates the degree-of-freedom mapper

    auto lambda = E_modulus * PoissonRatio / ( ( 1. + PoissonRatio ) * ( 1. - 2. * PoissonRatio ) );
    auto mu     = E_modulus / ( 2. * ( 1. + PoissonRatio ) );

    gsMatrix<> I = gsMatrix<>::Identity(3,3); // es 3x3?

    auto phys_jacobian = ijac(u, G);

    // gsDebugVar(ev.eval(phys_jacobian,pt,0).cols());
    // gsDebugVar(ev.eval(phys_jacobian,pt,0).rows());

    // Linear strains
    //auto EE = 0.5*(phys_jacobian.cwisetr() + phys_jacobian); //?? in matrix??? strain
    // Green lagrange strain
    auto EE = 0.5*(phys_jacobian.cwisetr() + phys_jacobian + phys_jacobian.cwisetr() * phys_jacobian); //?? in matrix??? strain
    
    // gsDebugVar(ev.eval(phys_jacobian,pt,0));
    // gsDebugVar(ev.eval(EE,pt,0).cols());
    // gsDebugVar(ev.eval(EE,pt,0).rows());

    auto EE_eval = ev.eval(EE,pt,0); // Strain (small deformations) evaluated at pt (x,y,z)
    gsDebugVar(EE_eval);
    
    auto SS = lambda * EE.trace().val()*I + 2.0*mu*EE; // stress
    auto SS_eval = ev.eval(SS,pt,0); // Stress evaluated at pt (x,y,z)
    gsDebugVar(SS_eval);

    // ============================================================
    //

    // ------------- Elasticity check (with expression assembler) -------------
    // // gsDebugVar( mp.patch(0).coefs());
    // // gsDebugVar( mp_def.patch(0).coefs());

    // auto matrix_disp =  mp_def.patch(0).coefs()- mp.patch(0).coefs();
    // gsDebugVar(matrix_disp);
    
    // // for (index_t q = 0; q < matrix_disp.cols(); ++q)
    // // {
    // //     auto dispGrad = matrix_disp.jacobian(q)*(matrix_disp.jacobian(q).cramerInverse());
    // // }
    // // auto eps = (dispGrad + dispGrad.transpose())/2;
    // gsMatrix<> result;
    // //auto holahola =  gsCauchyStressFunction<index_t>::linearElastic(matrix_disp,result);
    // // gsWriteParaview(mp,"holahola",1,true);

    
    
    
    
    // gsExprAssembler<> A(1,1);

    // typedef gsExprAssembler<>::geometryMap geometryMap;
    // typedef gsExprAssembler<>::variable    variable;
    // typedef gsExprAssembler<>::space       space;
    // typedef gsExprAssembler<>::solution    solution;
    
    // gsMultiBasis<> dbasis(mp,true);
    
    // space w = A.getSpace(dbasis);

    // A.setIntegrationElements(dbasis); // necesito escribir la basis!
    // gsExprEvaluator<> ev(A);
    // A.assembleJacobian()
    // geometryMap G = A.getMap(mp); // al mp o al mp_def?

    

    // // Set the discretization space
    // space u_trial = A.getSpace(dbasis, 3); //?

    // auto phys_jacobian = ijac(u_trial, G);
    // auto epsilon = 0.5*(phys_jacobian.cwisetr() + phys_jacobian);
    // gsDebugVar(epsilon);

    // auto matrix = lambda * Itr + mu * Idev;

    // epsilon.eval_into(pt);
    // -------------------------------------------------------------------------


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