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
#include <gsUtils/gsStopwatch.h>

using namespace gismo;

int main (int argc, char** argv)
{
    gsCmdLine cmd(".");

    try { cmd.getValues(argc,argv); } catch (int rv) { return rv; }

    // Make geometry
    gsMultiPatch<> mp, mp_def;
    mp.addPatch( gsNurbsCreator<>::BSplineCube(1) ); // degree
    mp.addAutoBoundaries();
    mp.uniformRefine(); // refines both
    mp_def = mp;
    mp_def.patch(0).coefs().col(0) *= 2; // deformation *2 in the x direction (col(0))

    real_t E_modulus = 210e9;
    real_t PoissonRatio = 0.3;
    gsFunctionExpr<> E(std::to_string(E_modulus),3);
    gsFunctionExpr<> nu(std::to_string(PoissonRatio),3);

    gsMatrix<> pts(3,6);
    pts.col(0)<<0.125,0.375,0.5;
    pts.col(1)<<0.375,0.125,0.5;
    pts.col(2)<<0.125,0.25,0.5;
    pts.col(3)<<0.25,0.125,0.5;
    pts.col(4)<<0.25,0.25,0.5;
    pts.col(5)<<0.5,0.5,0.5;

    gsMatrix<> Fres, Eres, Sres, Cres;

    gsLinearMaterial<real_t> SvK(E_modulus, PoissonRatio, 3);
    gsStopwatch time;
    // ============================================================
    time.restart();
    // Using precomputation of material matrices
    // Precompute the material matrix
    gsMaterialData<real_t> data;
    SvK.precompute(&mp,&mp_def,0,pts,data);
    SvK.eval_deformation_gradient_into(data,Fres);
    SvK.eval_strain_into(data,Eres);
    SvK.eval_stress_into(data,Sres);
    SvK.eval_matrix_into(data,Cres);
    gsInfo<<"Computation took "<<time.stop()<<" s\n";

    gsInfo<<"Deformation gradient\n";
    gsInfo<<Fres<<"\n";

    gsInfo<<"Strain\n";
    gsInfo<<Eres<<"\n";

    gsInfo<<"Stress\n";
    gsInfo<<Sres<<"\n";

    gsInfo<<"Material matrix\n";
    gsInfo<<Cres<<"\n";

    // ============================================================
    // Using gsMaterialEval
    time.restart();
    gsMaterialEval<real_t, gsMaterialOutput::F, true> SvK_F(&SvK, mp, mp_def);
    gsMaterialEval<real_t, gsMaterialOutput::E, true> SvK_E(&SvK, mp, mp_def);
    gsMaterialEval<real_t, gsMaterialOutput::S, true> SvK_S(&SvK, mp, mp_def);
    gsMaterialEval<real_t, gsMaterialOutput::C, true> SvK_C(&SvK, mp, mp_def);
    SvK_F.piece(0).eval_into(pts,Fres);
    SvK_E.piece(0).eval_into(pts,Eres);
    SvK_S.piece(0).eval_into(pts,Sres);
    SvK_C.piece(0).eval_into(pts,Cres);
    gsInfo<<"Computation took "<<time.stop()<<" s\n";

    // gsInfo<<"Deformation gradient\n";
    // gsInfo<<Fres<<"\n";

    gsInfo<<"Strain\n";
    gsInfo<<Eres<<"\n";

    gsInfo<<"Stress\n";
    gsInfo<<Sres<<"\n";

    gsInfo<<"Material matrix\n";
    gsInfo<<Cres<<"\n";

    return 0;
}