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
#include <gsElasticity/gsLinearDegradedMaterial.h>
#include <gsElasticity/gsLinearMaterial.h>
#include <gsElasticity/gsMaterialEval.h>
#include <gsElasticity/gsElasticityAssembler.h>
#include <gsElasticity/gsPhaseFieldAssembler.h>
#include <gsElasticity/gsPSOR.h>
#include <gsUtils/gsStopwatch.h>

using namespace gismo;
//! [Include namespace]

#define PRINT(w) std::setw(w)<<std::left

template <short_t dim, class T>
void solve(gsOptionList & materialParameters,
           gsOptionList & controlParameters,
           gsMultiPatch<T> & mp,
           gsMultiPatch<T> & damage,
           gsBoundaryConditions<T> & bc_u,
           gsBoundaryConditions<T> & bc_d,
           bool plot,
           index_t plotmod,
           std::string & outputdir);

int main(int argc, char *argv[])
{
    //! [Parse command line]
    bool plot = false;
    index_t numHRef = 4;
    index_t numElev = 1;
    index_t order = -1;
    index_t AT = -1;
    index_t plotmod = 1;
    index_t dimension = 2;
    std::string output;

    gsCmdLine cmd("Tutorial on solving a Linear Elasticity problem.");
    cmd.addInt("e", "numElev","Number of degree elevation steps to perform before solving",numElev);
    cmd.addInt("r", "numHRef","Number of Uniform h-refinement loops", numHRef);
    cmd.addInt("O", "order","Order of the basis functions", order);
    cmd.addInt("A", "AT","AT-1 or AT-2 model", AT);
    cmd.addInt("p", "plotmod","Modulo for plotting", plotmod);
    cmd.addInt("d", "dimension","Dimension of the problem", dimension);
    cmd.addSwitch("plot","Create a ParaView visualization file with the solution", plot);
    cmd.addString("o", "output", "Output directory", output);

    try { cmd.getValues(argc,argv); } catch (int rv) { return rv; }

    GISMO_ASSERT(order == 2 || order == 4, "Please specify the order of the model (2 or 4).");
    GISMO_ASSERT(AT == 1 || AT == 2, "Please specify the AT model (1 or 2).");
    GISMO_ASSERT(dimension == 2 || dimension == 3, "Please specify the dimension of the problem (2 or 3).");

    ///////////////////////////////////////////////////////////////////////////////////////
    //DEFINE PROBLEM PARAMETERS////////////////////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////
    if (output.empty())
        output = "./output/";

    std::string outputdir = output + gsFileManager::getNativePathSeparator();
    gsFileManager::mkdir(output);

    ///////////////////////////////////////////////////////////////////////////////////////
    // Create the geometry
    ///////////////////////////////////////////////////////////////////////////////////////
    gsMultiPatch<> mp;

    real_t L = 20.;
    real_t H = 1.;

    if (dimension == 2)
    {
        gsKnotVector<> kv_x(0,1,math::ceil(L/H)-1,2,1);
        gsKnotVector<> kv_y(0,1,math::ceil(H/L)-1,2,1);
        gsTensorBSplineBasis<2,real_t> tbasis(kv_x, kv_y);
        gsTensorBSpline<2,real_t> tb(tbasis,tbasis.anchors().transpose());
        tb.coefs().col(0) *= L;
        tb.coefs().col(1) *= H;
        mp.addPatch(tb);
    }
    else if (dimension == 3)
    {
        gsKnotVector<> kv_x(0,1,math::ceil(L/H)-1,2,1);
        gsKnotVector<> kv_y(0,1,math::ceil(H/L)-1,2,1);
        gsKnotVector<> kv_z(0,1,math::ceil(H/L)-1,2,1);
        gsTensorBSplineBasis<3,real_t> tbasis(kv_x, kv_y, kv_z);
        gsTensorBSpline<3,real_t> tb(tbasis,tbasis.anchors().transpose());
        tb.coefs().col(0) *= L;
        tb.coefs().col(1) *= H;
        tb.coefs().col(2) *= H;
        mp.addPatch(tb);
    }
    else
        GISMO_ERROR("Invalid dimension");

    if (plot) gsWriteParaview(mp,outputdir+"mp",10,true);

    mp.degreeIncrease(numElev);
    for (index_t i = 0; i<numHRef; ++i)
        mp.uniformRefine(1);

    //// Material parameters
    gsOptionList materialParameters;
    // Young's modulus [kN/mm^2]
    materialParameters.addReal("E", "Young's modulus", 100.*1e3);
    // Poisson's ratio [-]
    materialParameters.addReal("nu", "Poisson's ratio", 0.);
    // Toughness [kN/mm]
    materialParameters.addReal("Gc", "Toughness", 0.01*1e3);
    // Internal length [mm]
    materialParameters.addReal("l0", "Internal length", 0.125);
    // Order of the phase-field model
    materialParameters.addInt("order", "Order of the phase-field model", order);
    // AT model
    materialParameters.addInt("AT", "AT model", AT);

    //// Boundary control parameters
    gsOptionList controlParameters;
    // Maximum displacement [mm]
    controlParameters.addReal("uend", "Maximum displacement", 20e-2);
    // Initial displacement [mm]
    controlParameters.addReal("umin", "Initial displacement", 0.0);
    // Displacement step [mm]
    controlParameters.addReal("ustep", "Displacement step", 20e-2/2000);
    // Step transition [mm]
    controlParameters.addReal("utrans", "Step transition", 20e-2);
    // Step reduction factor [-]
    controlParameters.addReal("ured", "Step reduction factor", 1.0);
    // Maximum number of iterations
    controlParameters.addInt("maxIt", "Maximum number of iterations", 10000);
    // Maximum number of iterations for elasticity problem
    controlParameters.addInt("maxItEl", "Maximum number of iterations for elasticity problem", 100);
    // Maximum number of iterations for phase-field problem
    controlParameters.addInt("maxItPf", "Maximum number of iterations for phase-field problem", 1);
    // Tolerance for the elasticity problem
    controlParameters.addReal("tolEl", "Tolerance for the elasticity problem", 1e-6);
    // Tolerance for the phase-field problem
    controlParameters.addReal("tolPf", "Tolerance for the phase-field problem", 1e-6);
    // Staggered tolerance
    controlParameters.addReal("tol", "Tolerance for the staggered scheme", 1e-6);

    // Initialize the damage field
    gsMultiPatch<> damage;
    gsMatrix<> coefs(mp.basis(0).size(),1);
    coefs.setZero();
    // coefs.setRandom();
    // coefs.array() += 1.;
    // coefs.array() *= 1e-3;
    damage.addPatch(mp.basis(0).makeGeometry(give(coefs)));
    if (plot) gsWriteParaview(mp,damage,outputdir+"cini_approx",100000);

    // Boundary conditions
    gsBoundaryConditions<> bc_u;
    bc_u.setGeoMap(mp);

    gsBoundaryConditions<> bc_d;
    bc_d.addCondition(boundary::west,condition_type::dirichlet,0,0);
    bc_d.addCondition(boundary::east,condition_type::dirichlet,0,0);
    bc_d.setGeoMap(mp);

    ///////////////////////////////////////////////////////////////////////////////////////
    // Call the dimensional solver
    ///////////////////////////////////////////////////////////////////////////////////////
    switch (mp.domainDim())
    {
        case 2:
            solve<2>(materialParameters,controlParameters,mp,damage,bc_u,bc_d,plot,plotmod,outputdir);
            break;
        case 3:
            solve<3>(materialParameters,controlParameters,mp,damage,bc_u,bc_d,plot,plotmod,outputdir);
            break;
        default:
            GISMO_ERROR("Invalid domain dimension");
    }

    return 0;
} // end main

template <short_t dim, class T>
void solve(gsOptionList & materialParameters,
           gsOptionList & controlParameters,
           gsMultiPatch<T> & mp,
           gsMultiPatch<T> & damage,
           gsBoundaryConditions<T> & bc_u,
           gsBoundaryConditions<T> & bc_d,
           bool plot,
           index_t plotmod,
           std::string & outputdir)
{
    ////////////////////////////////////////////////////////////////////////////////////
    // Load parameters
    ////////////////////////////////////////////////////////////////////////////////////
    // Young's modulus [N/mm^2]
    T E = materialParameters.getReal("E");
    // Poisson's ratio [-]
    T nu = materialParameters.getReal("nu");
    // Toughness [N/mm]
    T Gc = materialParameters.getReal("Gc");
    // Internal length [mm]
    T l0 = materialParameters.getReal("l0");
    // Order of the phase-field model
    index_t order = materialParameters.getInt("order");
    // AT model
    index_t AT = materialParameters.getInt("AT");

    GISMO_ASSERT(order == 2 || order == 4, "Please specify the order of the model (2 or 4).");
    GISMO_ASSERT(AT == 1 || AT == 2, "Please specify the AT model (1 or 2).");

    ////////////////////////////////////////////////////////////////////////////////////
    // Boundary control parameters
    ////////////////////////////////////////////////////////////////////////////////////
    // Maximum displacement [mm]
    T uend = controlParameters.getReal("uend");
    // Displacement step [mm]
    T ustep = controlParameters.getReal("ustep");
    // Initial displacement [mm]
    T ucurr = controlParameters.getReal("umin");
    // Step transition [mm]
    T utrans = controlParameters.askReal("utrans",uend);
    // Step reduction factor [-]
    T ured = controlParameters.askReal("ured",1.);
    // Maximum number of iterations
    index_t maxIt = controlParameters.getInt("maxIt");
    // Maximum number of iterations for elasticity problem
    index_t maxItEl = controlParameters.getInt("maxItEl");
    // Maximum number of iterations for phase-field problem
    index_t maxItPf = controlParameters.getInt("maxItPf");
    // Tolerance for the elasticity problem
    T tolEl = controlParameters.getReal("tolEl");
    // Tolerance for the phase-field problem
    T tolPf = controlParameters.getReal("tolPf");
    // Staggered tolerance
    T tol = controlParameters.getReal("tol");
    // // Fixed side patch id
    // index_t fixedSidePatch = controlParameters.getInt("patchId");
    // // Fixed side id
    // index_t fixedSideId = controlParameters.getInt("side");
    // // Fixed side direction
    // index_t fixedSideDir = controlParameters.getInt("direction");
    // Function on the boundary
    std::string bcFunctionLeft = controlParameters.askString("function","-u");
    std::string bcFunctionRight = controlParameters.askString("function","u");

    ///////////////////////////////////////////////////////////////////////////////////////
    //PROBLEM SETUP////////////////////////////////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////

    // Construct the basis
    gsMultiBasis<T> mb(mp);
    gsInfo<<"The basis has size "<<mb.size()<<" and degree "<<mb.degree()<<"\n";
    for (size_t b=0; b!=mb.nBases(); b++)
        gsInfo<<"Basis "<<b<<":\n"<<mb.basis(b)<<"\n";

    // Boundary conditions
    gsFunctionExpr<T> displ_left(bcFunctionLeft,dim);
    gsFunctionExpr<T> displ_right(bcFunctionRight,dim);
    displ_left.set_u(ucurr);
    displ_right.set_u(ucurr);
    bc_u.addCondition(boundary::west,condition_type::dirichlet,&displ_left ,0);
    bc_u.addCondition(boundary::east,condition_type::dirichlet,&displ_right,0);
    bc_u.addCondition(boundary::south,condition_type::dirichlet,0,1); //vertical constraint
    if (dim==3)
        bc_u.addCondition(boundary::back,condition_type::dirichlet,0,2); //vertical constraint

    bc_u.setGeoMap(mp);

    bc_d.setGeoMap(mp);

    ///////////////////////////////////////////////////////////////////////////////////////
    //INITIALIZATION///////////////////////////////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////

    // Initialize the solution (deformed geometry)
    gsMultiPatch<T> mp_def = mp;
    gsMultiPatch<T> displacement = mp;
    for (index_t p=0; p<mp.nPatches(); ++p)
        displacement.patch(p).coefs().setZero();

    // Initialize the material
    gsLinearDegradedMaterial<T> material(E,nu,damage,dim);
    // Initialize the elasticity assembler
    gsVector<T> bodyForceVec(dim);
    bodyForceVec.setZero();
    gsConstantFunction<T> bodyForce(bodyForceVec,dim);
    gsElasticityAssembler<T> elAssembler(mp,mb,bc_u,bodyForce,&material);
    elAssembler.options().setReal("quA",1.0);
    elAssembler.options().setInt ("quB",0);
    elAssembler.options().setSwitch ("SmallStrains",true);
    std::vector<gsMatrix<T> > fixedDofs = elAssembler.allFixedDofs();
    // elAssembler.assemble();
    gsBoundaryConditions<T> bc_u_dummy;
    bc_u_dummy.setGeoMap(mp);
    gsElasticityAssembler<T> fullElAssembler(mp,mb,bc_u_dummy,bodyForce,&material);
    fullElAssembler.options().setReal("quA",1.0);
    fullElAssembler.options().setInt ("quB",0);
    fullElAssembler.options().setSwitch ("SmallStrains",true);
    std::vector<gsMatrix<T> > dummyFixedDofs = fullElAssembler.allFixedDofs();

    // Initialize the phase-field assembler
    gsPhaseFieldAssemblerBase<T> * pfAssembler;
    if      (order == 2 && AT == 1)
        pfAssembler = new gsPhaseFieldAssembler<T,PForder::Second,PFmode::AT1>(mp,mb,bc_d);
    else if (order == 4 && AT == 1)
    {
        pfAssembler = new gsPhaseFieldAssembler<T,PForder::Fourth,PFmode::AT1>(mp,mb,bc_d);
        pfAssembler->options().setReal("cw",4.44847);
    }
    else if (order == 2 && AT == 2)
        pfAssembler = new gsPhaseFieldAssembler<T,PForder::Second,PFmode::AT2>(mp,mb,bc_d);
    else if (order == 4 && AT == 2)
        pfAssembler = new gsPhaseFieldAssembler<T,PForder::Fourth,PFmode::AT2>(mp,mb,bc_d);
    else
        GISMO_ERROR("Invalid order and/or AT model");

    pfAssembler->options().setReal("l0",l0);
    pfAssembler->options().setReal("Gc",Gc);
    pfAssembler->options().setReal("ExprAssembler.quA",1.0);
    pfAssembler->options().setInt ("ExprAssembler.quB",0);
    pfAssembler->initialize();

    //////////////////////////////////////////////////////////////////////////
    // SOLVE THE PROBLEM/////////////////////////////////////////////////////
    //////////////////////////////////////////////////////////////////////////

    gsMatrix<T> u(elAssembler.numDofs(),1);
    u.setZero();

#ifdef GISMO_WITH_PARDISO
    typename gsSparseSolver<T>::PardisoLDLT solver;
#else
    typename gsSparseSolver<T>::CGDiagonal solver;
#endif
    gsSparseMatrix<T> K;
    gsMatrix<T> R, F;

    T elAssemblyTime = 0.0;
    T elSolverTime = 0.0;
    T pfAssemblyTime = 0.0;
    T pfSolverTime = 0.0;
    T iterationTime  = 0.0;

    gsMatrix<T> D, deltaD;
    gsSparseMatrix<T> Q, QPhi, QPsi;
    gsMatrix<T> q, qpsi;
    // Phase-field assembly can already be performed since some terms are independent of the solutions
    pfAssembler->assembleMatrix();
    pfAssembler->matrix_into(QPhi);
    pfAssembler->assembleVector();
    pfAssembler->rhs_into(q);

    index_t step = 0;

    gsParaviewCollection damageCollection(outputdir+"damage");
    gsParaviewCollection psiCollection(outputdir+"Psi");
    gsParaviewCollection displCollection(outputdir+"displacement");
    gsStopwatch smallClock, bigClock;

    // pfAssembler->assembleMatrix();
    // pfAssembler->matrix_into(QPhi);
    // pfAssembler->constructSolution(damage,D);
    // gsInfo<<"D_0 = "<<(0.5 * D.transpose() * QPhi * D).value()<<"\n";
    //     gsWriteParaview(mp,damage,"damage_ini",100000);

    std::ofstream file(outputdir+"results.txt");
    file<<"u,Fx,Fy,E_u,E_d\n";
    file.close();

    real_t Rnorm, Fnorm;
    Rnorm = Fnorm = 1;
    while (ucurr<=uend)
    {
        // Update the boundary conditions
        displ_left.set_u(ucurr);
        displ_right.set_u(ucurr);
        elAssembler.computeDirichletDofs(0); // NOTE: This computes the DDofs for **unknown** 1, which should be component 1. This is a bug in the gsElasticity assembler
        fixedDofs = elAssembler.allFixedDofs();
        elAssembler.setFixedDofs(fixedDofs);

        gsInfo<<"---------------------------------------------------------------------------------------------------------------------------\n";
        gsInfo<<"Load step "<<step<<": u = "<<ucurr<<"\n\n";

        deltaD.setZero();
        for (index_t it=0; it!=maxIt; ++it)
        {
            elAssemblyTime = elSolverTime = 0.0;
            pfAssemblyTime = pfSolverTime = 0.0;
            iterationTime  = 0.0;
            bigClock.restart();
            gsInfo<<" - Staggered iteration "<<it<<":\n";
            gsInfo<<"\t"<<PRINT(20)<<"* Elasticity:"<<PRINT(6)<<"It."<<PRINT(14)<<"||R||"<<PRINT(14)<<"||dU||/||U||"<<PRINT(20)<<"cum. assembly [s]"<<PRINT(20)<<"cum. solver [s]"<<"\n";

            material.setParameter(2,damage);
            // Pre-assemble the elasticity problem
            smallClock.restart();
            elAssembler.assemble(u,fixedDofs);
            K = elAssembler.matrix();
            F = elAssembler.rhs();
            Fnorm = (F.norm() == 0) ? 1 : F.norm();
            elAssemblyTime += smallClock.stop();
            for (index_t elIt=0; elIt!=maxItEl; ++elIt)
            {
                // Solve
                smallClock.restart();
                solver.compute(K);
                u = solver.solve(F);
                elSolverTime += smallClock.stop();

                smallClock.restart();
                elAssembler.assemble(u,fixedDofs);
                K = elAssembler.matrix();
                F = elAssembler.rhs();
                elAssemblyTime += smallClock.stop();
                Rnorm = (K*u - F).norm();
                gsInfo<<"\t"<<PRINT(20)<<""<<PRINT(6)<<elIt<<PRINT(14)<<Rnorm/Fnorm<<PRINT(14)<<u.norm()<<PRINT(20)<<elAssemblyTime<<PRINT(20)<<elSolverTime<<"\n";

                if (Rnorm/* /F.norm() */ < tolEl || u.norm() < 1e-12)
                    break;
                else if (elIt == maxItEl-1 && maxItEl != 1)
                    GISMO_ERROR("Elasticity problem did not converge.");
            }

            elAssembler.setFixedDofs(fixedDofs);
            elAssembler.constructSolution(u,fixedDofs,displacement);
            for (size_t p=0; p!=mp.nPatches(); ++p)
                mp_def.patch(p).coefs() = mp.patch(p).coefs() + displacement.patch(p).coefs();

            // Initialize the function for the elastic energy
            gsMaterialEval<T,gsMaterialOutput::Psi> Psi(&material,mp,mp_def);

            // ==================================================================================

            gsInfo<<"\t"<<PRINT(20)<<"* Phase-field:"<<PRINT(6)<<"It."<<PRINT(14)<<"||R||"<<PRINT(14)<<"||dD||/||D||"<<PRINT(20)<<"cum. assembly [s]"<<PRINT(20)<<"cum. solver [s]"<<"\n";

            // Phase-field problem
            smallClock.restart();
            // gsInfo<<"Assembling phase-field problem"<<"\n";
            pfAssembler->assemblePsiMatrix(Psi);
            pfAssembler->matrix_into(QPsi);
            pfAssembler->assemblePsiVector(Psi);
            pfAssembler->rhs_into(qpsi);
            Q = QPhi + QPsi;
            // Reconstruct the solution from the damage field
            pfAssembler->constructSolution(damage,D);
            pfAssemblyTime = smallClock.stop();
            // gsInfo<<". Done\n";

            // Initialize the PSOR solver
            smallClock.restart();
            gsPSOR<T> PSORsolver(Q);
            PSORsolver.options().setInt("MaxIterations",30000);
            PSORsolver.options().setSwitch("Verbose",false);
            PSORsolver.options().setReal("tolU",1e-4);
            PSORsolver.options().setReal("tolNeg",1e-9);
            PSORsolver.options().setReal("tolPos",1e-9);
            pfSolverTime = smallClock.stop();
            for (index_t pfIt=0; pfIt!=maxItPf; ++pfIt)
            {
                // Assemble
                smallClock.restart();
                R = Q * D - qpsi + q;
                pfAssemblyTime += smallClock.stop();

                // Solve
                smallClock.restart();
                // solver.compute(Q);
                // deltaD = solver.solve(-R);
                // gsDebugVar(deltaD.norm());
                PSORsolver.solve(R,deltaD); // deltaD = Q \ R
                pfSolverTime += smallClock.stop();
                D += deltaD;

                gsInfo<<"\t"<<PRINT(20)<<""<<PRINT(6)<<pfIt<<PRINT(14)<<R.norm()<<PRINT(14)<<deltaD.norm()/D.norm()<<PRINT(20)<<pfAssemblyTime<<PRINT(20)<<pfSolverTime<<"\n";;
                if (deltaD.norm()/D.norm() < tolPf || D.norm() < 1e-12)
                    break;
                else if (pfIt == maxItPf-1 && maxItPf != 1)
                    GISMO_ERROR("Phase-field problem did not converge.");
            }

            // Update damage spline
            pfAssembler->constructSolution(D,damage);

            smallClock.restart();
            material.setParameter(2,damage);
            elAssembler.assemble(u,fixedDofs);
            K = elAssembler.matrix();
            F = elAssembler.rhs();
            Fnorm = (F.norm() == 0) ? 1 : F.norm();
            elAssemblyTime += smallClock.stop();
            Rnorm = (K*u - F).norm();

            iterationTime = bigClock.stop();
            gsInfo<<"\t"<<PRINT(20)<<"* Finished"<<PRINT(6)<<""<<PRINT(14)<<"||R||"<<PRINT(14)<<"||R||/||F||"<<PRINT(14)<<"total [s]"<<PRINT(20)<<"elasticity [s]"           <<PRINT(20)<<"phase-field [s]"          <<"\n";
            gsInfo<<"\t"<<PRINT(20)<<""          <<PRINT(6)<<""<<PRINT(14)<<Rnorm<<PRINT(14)<<Rnorm/Fnorm<<PRINT(14)<<iterationTime<<PRINT(20)<<elAssemblyTime+elSolverTime<<PRINT(20)<<pfAssemblyTime+pfSolverTime<<"\n";
            if (Rnorm/Fnorm < tol)
                break;
            else if (it == maxIt-1)
                GISMO_ERROR("Staggered iterations problem did not converge.");
        }

        // =========================================================================
        // Compute resulting force and energies
        gsMatrix<T> ufull = displacement.patch(0).coefs().reshape(displacement.patch(0).coefs().size(),1);
        fullElAssembler.assemble(ufull,dummyFixedDofs);
        gsMatrix<T> Rfull = fullElAssembler.matrix()*ufull - fullElAssembler.rhs();
        // sum the reaction forces in Y direction
        gsDofMapper mapper(mb,dim);
        mapper.finalize();
        gsMatrix<index_t> boundary = mb.basis(0).boundary(boundary::east);
        T Fx = 0, Fy = 0;
        for (index_t k=0; k!=boundary.size(); k++)
        {
            Fx += Rfull(mapper.index(boundary(k,0),0,0),0); // DoF index, patch, component
            Fy += Rfull(mapper.index(boundary(k,0),0,1),0); // DoF index, patch, component
        }

        std::vector<T> stepData(5);
        stepData[0] = ucurr;
        stepData[1] = Fx;
        stepData[2] = Fy;
        stepData[3] = (0.5 * ufull.transpose() * fullElAssembler.matrix() * ufull).value();
        stepData[4] = (0.5 * D.transpose() * QPhi * D).value() + (D.transpose() * q).value();

        gsInfo<<"\n";
        gsInfo<<"Converged with ||R||/||F|| = "<<Rnorm/Fnorm<<" < "<<tol<<" ||D|| = "<<D.norm()<<" ||U|| = "<<u.norm()<<"\n";
        // gsInfo<<"----------------------------------------------------------------------------------------------------\n\n";

        // =========================================================================
        // PLOT
        if (plot && step%plotmod==0)
        {
            gsMatrix<> eval_geo, eval_damage, eval_psi, eval_displacement, pts, ab;
            gsVector<> a, b;
            ab = mp.patch(0).support();
            a  = ab.col(0);
            b  = ab.col(1);
            // Generate a point grid
            gsVector<unsigned> np(dim);
            np[0] = 1000;
            np.segment(1,dim-1).setConstant(2);
            pts = gsPointGrid(a,b,np);
            // Evaluate the geometry
            eval_geo = mp_def.patch(0).eval(pts);

            std::string filename;
            filename = "damage_" + util::to_string(step);
            eval_damage = damage.patch(0).eval(pts);
            gsWriteParaviewTPgrid(eval_geo,eval_damage,np.template cast<index_t>(),outputdir+filename);
            // gsWriteParaview(mp,damage,outputdir+filename,100000);
            filename += "0";
            damageCollection.addPart(filename,step,"Solution",0);

            gsMaterialEval<T,gsMaterialOutput::Psi> Psi(&material,mp,mp_def);
            filename = "Psi_"+util::to_string(step);
            eval_psi = Psi.piece(0).eval(pts);
            gsWriteParaviewTPgrid(eval_geo,eval_psi,np.template cast<index_t>(),outputdir+filename);
            // gsWriteParaview(mp,Psi,outputdir+filename,100000);
            filename += "0";
            psiCollection.addPart(filename,step,"Solution",0);

            filename = "displacement_"+util::to_string(step);
            eval_displacement = displacement.patch(0).eval(pts);
            gsWriteParaviewTPgrid(eval_geo,eval_displacement,np.template cast<index_t>(),outputdir+filename);
            // gsWriteParaview(mp,displacement,outputdir+filename,1000);
            filename += "0";
            displCollection.addPart(filename,step,"Solution",0);
        }

        // =========================================================================
        // Write data
        std::ofstream file(outputdir+"results.txt",std::ios::app);
        // for (size_t i = 0; i != data.size(); ++i)
        //     file<<data[i][0]<<","<<-data[i][1]<<","<<-data[i][2]<<","<<data[i][3]<<","<<data[i][4]<<"\n";
        file<<stepData[0]<<","<<-stepData[1]<<","<<-stepData[2]<<","<<stepData[3]<<","<<stepData[4]<<"\n";
        file.close();

        ucurr += (ucurr+ustep > utrans) ? ustep/ured : ustep;
        step++;
    }

    if (plot)
    {
        damageCollection.save();
        psiCollection.save();
        displCollection.save();
    }

    delete pfAssembler;
}