/** @file gsMaterialMatrix_test.cpp

    @brief Simple example for material matrix evaluations

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): H.M. Verhelst
*/

#include <gismo.h>
#include <gsElasticity/gsMaterialBase.h>
#include <gsElasticity/gsLinearMaterial.h>
#include <gsElasticity/gsNeoHookeLogMaterial.h>

using namespace gismo;

int main (int argc, char** argv)
{
    index_t testCase=0;
    gsCmdLine cmd(".");
    cmd.addInt("t","testCase","testCase",testCase);

    try { cmd.getValues(argc,argv); } catch (int rv) { return rv; }

    // Make geometry
    gsMultiPatch<> mp, mp_def;
    mp.addPatch( gsNurbsCreator<>::BSplineSquare(1) ); // degree
    mp.addAutoBoundaries();

    if (testCase==0)
    {
        mp_def = mp;
        mp_def.patch(0).coefs().col(0) *= 2;
        mp_def.patch(0).coefs().col(1) *= 1.5;

    }
    else if (testCase==1)
    {
        mp_def = mp;
        mp_def.patch(0).coefs()(1,1) += 0.5;
        mp_def.patch(0).coefs()(3,1) += 0.5;
    }
    gsDebugVar(mp_def.patch(0).coefs());

    real_t thickness = 1.0;
    real_t E_modulus = 210e9;
    real_t PoissonRatio = 0.3;
    real_t PoissonRatio2 = 0.5;
    real_t Ratio = 1;
    gsFunctionExpr<> t(std::to_string(thickness), 2);
    gsFunctionExpr<> E(std::to_string(E_modulus),2);
    gsFunctionExpr<> nu(std::to_string(PoissonRatio),2);
    gsFunctionExpr<> nu2(std::to_string(PoissonRatio2),2);
    gsConstantFunction<> ratio(Ratio,2);

    gsLinearMaterial<2,real_t> SvK(E_modulus, PoissonRatio, &mp, &mp_def);
    gsNeoHookeLogMaterial<2,real_t> NH(E_modulus, PoissonRatio, &mp, &mp_def);
    gsMaterialBase<real_t> base(&mp, &mp_def);



    gsVector<> pt(2);
    pt.col(0)<<0.125,0.375;

    gsMatrix<> pts(2,6);
    pts.col(0)<<0.125,0.375;
    pts.col(1)<<0.375,0.125;
    pts.col(2)<<0.125,0.25;
    pts.col(3)<<0.25,0.125;
    pts.col(4)<<0.25,0.25;
    pts.col(5)<<0.5,0.5;

    gsMatrix<> Fres, Eres, Sres;
    base.gradient_into(pt,Fres);
    base.strain_into(pt,Eres);

    SvK.eval_vector_into(pt,Sres,0);
    gsDebugVar(Sres);

    NH.eval_vector_into(pt,Sres,0);
    gsDebugVar(Sres);


    return 0;
}