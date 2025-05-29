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
    index_t plotmod = 1;
    index_t numHRef = 0;
    index_t numElev = 0;
    std::string output;
    std::string damageInput;
    std::string geoInput;
    std::string parInput;
    std::string inputDir;
    index_t testCase = 0;

    gsCmdLine cmd("Tutorial on solving a Linear Elasticity problem.");
    cmd.addInt("e", "numElev","Number of degree elevation steps to perform before solving",numElev);
    cmd.addInt("r", "numHRef","Number of Uniform h-refinement loops", numHRef);
    cmd.addInt("p", "plotmod","Modulo for plotting", plotmod);
    cmd.addSwitch("plot","Create a ParaView visualization file with the solution", plot);
    cmd.addString("o", "output", "Output directory", output);
    cmd.addString("d", "damage", "Damage file", damageInput);
    cmd.addString("g", "geometry", "Geometry file", geoInput);
    cmd.addString("i", "parInput", "Input XML file for model parameters", parInput);
    cmd.addString("I", "inputDir", "Input directory", inputDir);

    try { cmd.getValues(argc,argv); } catch (int rv) { return rv; }

    inputDir = inputDir + gsFileManager::getNativePathSeparator();
    std::string parInputPath = (parInput.empty() ? inputDir + "parameters.xml" : parInput);
    std::string geoInputPath = (geoInput.empty() ? inputDir + "geometry.xml" : geoInput);
    std::string damageInputPath = (damageInput.empty() ? inputDir + "damage.xml" : damageInput);
    GISMO_ASSERT(gsFileManager::fileExists(parInputPath), "Input parameter file "<<parInputPath<<" not found.");
    GISMO_ASSERT(gsFileManager::fileExists(geoInputPath), "Input geometry file "<<geoInputPath<<" not found.");
    GISMO_ASSERT(gsFileManager::fileExists(damageInputPath), "Input damage file "<<damageInputPath<<" not found.");
    gsInfo << "Input parameter file "<< parInputPath <<"\n";
    gsInfo << "Input geometry file "<< geoInputPath <<"\n";
    gsInfo << "Input damage file "<< damageInputPath <<"\n";


    ///////////////////////////////////////////////////////////////////////////////////////
    //DEFINE PROBLEM PARAMETERS////////////////////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////
    if (output.empty())
        output = "./output/";

    std::string outputdir = output + gsFileManager::getNativePathSeparator();
    gsFileManager::mkdir(output);

    gsFileData<> fd_geo(geoInput.empty() ? inputDir + "geometry.xml" : geoInput);
    gsMultiPatch<> mp;
    fd_geo.getFirst(mp);
    if (numElev > 0)
        mp.degreeIncrease(numElev);
    for (index_t i = 0; i<numHRef; ++i)
        mp.uniformRefine(1);
    if (plot) gsWriteParaview(mp,outputdir+"mp",10,true);

    gsFileData<> fd_damage(damageInput.empty() ? inputDir + "damage.xml" : damageInput);
    gsMultiPatch<> damage;
    fd_damage.getFirst(damage);
    if (plot) gsWriteParaview(mp,damage,outputdir+"initial_damage",100000);

    gsFileData<> fd_pars(parInput.empty() ? inputDir + "parameters.xml" : parInput);
    GISMO_ASSERT(fd_pars.hasLabel("material"), "Material parameters not found in the input file.");
    GISMO_ASSERT(fd_pars.hasLabel("control"), "Displacement-Control parameters not found in the input file.");
    GISMO_ASSERT(fd_pars.hasLabel("BCs_u"), "Displacement boundary conditions not found in the input file.");
    GISMO_ASSERT(fd_pars.hasLabel("BCs_d"), "Phase-field boundary conditions not found in the input file.");

    //// Material parameters
    gsOptionList materialParameters;
    fd_pars.getLabel("material", materialParameters);

    //// Boundary control parameters
    gsOptionList controlParameters;
    fd_pars.getLabel("control", controlParameters);

    //// Boundary conditions
    gsBoundaryConditions<> bc_u;
    fd_pars.getLabel("BCs_u", bc_u);

    //// Boundary conditions
    gsBoundaryConditions<> bc_d;
    fd_pars.getLabel("BCs_d", bc_d);

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
    // Fixed side patch id
    index_t fixedSidePatch = controlParameters.getInt("patchId");
    // Fixed side id
    index_t fixedSideId = controlParameters.getInt("side");
    // Fixed side direction
    index_t fixedSideDir = controlParameters.getInt("direction");
    // Function on the boundary
    std::string bcFunction = controlParameters.askString("function","u");

    ///////////////////////////////////////////////////////////////////////////////////////
    //PROBLEM SETUP////////////////////////////////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////

    // Construct the basis
    gsMultiBasis<T> mb(mp);
    gsInfo<<"The basis has size "<<mb.size()<<" and degree "<<mb.degree()<<"\n";
    for (size_t b=0; b!=mb.nBases(); b++)
        gsInfo<<"Basis "<<b<<":\n"<<mb.basis(b)<<"\n";

    // Boundary conditions
    gsFunctionExpr<T> displ(bcFunction,dim);
    displ.set_u(ucurr);
    bc_u.addCondition(fixedSideId,condition_type::dirichlet,&displ,fixedSideDir);
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
    T elAssemblyTime = 0.0;
    T elSolverTime = 0.0;
    T pfAssemblyTime = 0.0;
    T pfSolverTime = 0.0;
    T iterationTime  = 0.0;

    gsMatrix<T> R;
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
    file<<"u,Fx,Fy,E_u,E_d,elAssemblyTime,elSolverTime,pfAssemblyTime,pfSolverTime\n";
    file.close();

    T Rnorm, Fnorm;
    Rnorm = Fnorm = 1;
    while (ucurr<=uend)
    {
        // Update the boundary conditions
        displ.set_u(ucurr);
        elAssembler.computeDirichletDofs(fixedSideDir); // NOTE: This computes the DDofs for **unknown** 1, which should be component 1. This is a bug in the gsElasticity assembler
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
            gsInfo<<"\t"<<PRINT(20)<<"* Elasticity:"<<PRINT(6)<<"It."<<PRINT(14)<<"||R||"<<PRINT(14)<<"||F||"<<PRINT(14)<<"||R||/||F||"<<PRINT(14)<<"||U||"<<PRINT(20)<<"cum. assembly [s]"<<PRINT(20)<<"cum. solver [s]"<<"\n";

            material.setParameter(2,damage);
            // Pre-assemble the elasticity problem
            smallClock.restart();
            elAssembler.assemble(u,fixedDofs);
            Fnorm = elAssembler.rhs().norm();
            Fnorm = (Fnorm == 0) ? 1 : Fnorm;
            elAssemblyTime += smallClock.stop();
            for (index_t elIt=0; elIt!=maxItEl; ++elIt)
            {
                // Solve
                smallClock.restart();
                solver.compute(elAssembler.matrix());
                u = solver.solve(elAssembler.rhs());
                elSolverTime += smallClock.stop();

                smallClock.restart();
                elAssembler.assemble(u,fixedDofs);
                Fnorm = elAssembler.rhs().norm();
                Fnorm = (Fnorm == 0) ? 1 : Fnorm;
                Rnorm = (elAssembler.matrix()*u - elAssembler.rhs()).norm();
                elAssemblyTime += smallClock.stop();
                gsInfo<<"\t"<<PRINT(20)<<""<<PRINT(6)<<elIt<<PRINT(14)<<Rnorm<<PRINT(14)<<Fnorm<<PRINT(14)<<Rnorm/Fnorm<<PRINT(14)<<u.norm()<<PRINT(20)<<elAssemblyTime<<PRINT(20)<<elSolverTime<<"\n";

                // Evaluate the geometry in the support
                gsMatrix<> supp(2,2);
                supp.col(0)<<0.0,0.48;
                supp.col(1)<<0.6,0.52;
                gsVector<unsigned> npts = uniformSampleCount<real_t>(supp.col(0), supp.col(1), 100000);
                gsMatrix<> points = gsPointGrid<real_t>(supp.col(0),supp.col(1),npts);
                gsMatrix<> eval_geo, eval_damage;
                mp.piece(0).eval_into(points, eval_geo);
                damage.piece(0).eval_into(points, eval_damage);
                gsWriteParaviewTPgrid(eval_geo,eval_damage,npts.template cast<index_t>(),outputdir+"damage_tmp");

                if (Rnorm/Fnorm < tolEl || u.norm() < 1e-12)
                    break;
                else if (elIt == maxItEl-1 && maxItEl != 1)
                    GISMO_ERROR("Elasticity problem did not converge.");
            }

            elAssembler.setFixedDofs(fixedDofs);
            elAssembler.constructSolution(u,fixedDofs,displacement);
            for (size_t p=0; p!=mp.nPatches(); ++p)
                mp_def.patch(p).coefs() = mp.patch(p).coefs() + displacement.patch(p).coefs();

            // Initialize the function for the elastic energy
            gsMaterialEval<T,gsMaterialOutput::Psi,true,true> Psi(&material,mp,mp_def);

            // ==================================================================================

            gsInfo<<"\t"<<PRINT(20)<<"* Phase-field:"<<PRINT(6)<<"It."<<PRINT(14)<<"||R||"<<PRINT(14)<<"||D||"<<PRINT(14)<<"||dD||/||D||"<<PRINT(20)<<"cum. assembly [s]"<<PRINT(20)<<"cum. solver [s]"<<"\n";

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
                gsInfo<<"\t"<<PRINT(20)<<""<<PRINT(6)<<pfIt<<PRINT(14)<<R.norm()<<PRINT(14)<<D.norm()<<PRINT(14)<<deltaD.norm()/D.norm()<<PRINT(20)<<pfAssemblyTime<<PRINT(20)<<pfSolverTime<<"\n";;
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
            Fnorm = elAssembler.rhs().norm();
            Fnorm = (Fnorm == 0) ? 1 : Fnorm;
            Rnorm = (elAssembler.matrix()*u - elAssembler.rhs()).norm();
            elAssemblyTime += smallClock.stop();

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
        gsDofMapper mapper(mb,dim);
        mapper.finalize();
        gsMatrix<index_t> boundary = mb.basis(0).boundary(fixedSideId);
        T Fx = 0, Fy = 0;
        for (index_t k=0; k!=boundary.size(); k++)
        {
            Fx += Rfull(mapper.index(boundary(k,0),0,0),0); // DoF index, patch, component
            Fy += Rfull(mapper.index(boundary(k,0),0,1),0); // DoF index, patch, component
        }

        std::vector<T> stepData(9);
        stepData[0] = ucurr;
        stepData[1] = Fx;
        stepData[2] = Fy;
        stepData[3] = (0.5 * ufull.transpose() * fullElAssembler.matrix() * ufull).value();
        stepData[4] = (0.5 * D.transpose() * QPhi * D).value() + (D.transpose() * q).value();
        stepData[5] = elAssemblyTime;
        stepData[6] = elSolverTime;
        stepData[7] = pfAssemblyTime;
        stepData[8] = pfSolverTime;
        gsInfo<<"\n";
        gsInfo<<"Converged with ||R||/||F|| = "<<Rnorm/Fnorm<<" < "<<tol<<" ||D|| = "<<D.norm()<<" ||U|| = "<<u.norm()<<"\n";
        // gsInfo<<"----------------------------------------------------------------------------------------------------\n\n";

        // =========================================================================
        // PLOT
        if (plot && step%plotmod==0)
        {
            std::string filename;
            filename = "damage_" + util::to_string(step);
            if (plot) gsWriteParaview(mp,damage,outputdir+filename,100000);
            // gsField<T> damage_step(zone,damage,false);
            // gsWriteParaview(damage_step,filename,1000);
            filename += "0";
            damageCollection.addPart(filename,step,"Solution",0);

            gsMaterialEval<T,gsMaterialOutput::Psi> Psi(&material,mp,mp_def);
            filename = "Psi_"+util::to_string(step);
            if (plot) gsWriteParaview(mp,Psi,outputdir+filename,100000);
            filename += "0";
            psiCollection.addPart(filename,step,"Solution",0);

            filename = "displacement_"+util::to_string(step);
            if (plot) gsWriteParaview(mp,displacement,outputdir+filename,1000);
            filename += "0";
            displCollection.addPart(filename,step,"Solution",0);
        }

        // =========================================================================
        // Write data
        std::ofstream file(outputdir+"results.txt",std::ios::app);
        // for (size_t i = 0; i != data.size(); ++i)
        //     file<<data[i][0]<<","<<-data[i][1]<<","<<-data[i][2]<<","<<data[i][3]<<","<<data[i][4]<<"\n";
        file<<stepData[0]<<","<<-stepData[1]<<","<<-stepData[2]<<","<<stepData[3]<<","<<stepData[4]<<","<<stepData[5]<<","<<stepData[6]<<","<<stepData[7]<<","<<stepData[8]<<"\n";
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

