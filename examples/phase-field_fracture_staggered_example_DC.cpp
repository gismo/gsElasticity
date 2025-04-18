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

int main(int argc, char *argv[])
{
    //! [Parse command line]
    bool plot = false;
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
    // Young's modulus [N/mm^2]
    real_t E = materialParameters.getReal("E");
    // Poisson's ratio [-]
    real_t nu = materialParameters.getReal("nu");
    // Toughness [N/mm]
    real_t Gc = materialParameters.getReal("Gc");
    // Internal length [mm]
    real_t l0 = materialParameters.getReal("l0");
    // Order of the phase-field model
    index_t order = materialParameters.getInt("order");
    // AT model
    index_t AT = materialParameters.getInt("AT");

    GISMO_ASSERT(order == 2 || order == 4, "Please specify the order of the model (2 or 4).");
    GISMO_ASSERT(AT == 1 || AT == 2, "Please specify the AT model (1 or 2).");

    //// Boundary control parameters
    gsOptionList controlParameters;
    fd_pars.getLabel("control", controlParameters);
    // Maximum displacement [mm]
    real_t uend = controlParameters.getReal("uend");
    // Displacement step [mm]
    real_t ustep = controlParameters.getReal("ustep");
    // Initial displacement [mm]
    real_t ucurr = controlParameters.getReal("umin");
    // Step transition [mm]
    real_t utrans = controlParameters.askReal("utrans",uend);
    // Step reduction factor [-]
    real_t ured = controlParameters.askReal("ured",1.);
    // Maximum number of iterations
    index_t maxIt = controlParameters.getInt("maxIt");
    // Tolerance
    real_t tol = controlParameters.getReal("tol");
    // Fixed side patch id
    index_t fixedSidePatch = controlParameters.getInt("patchId");
    // Fixed side id
    index_t fixedSideId = controlParameters.getInt("side");
    // Fixed side direction
    index_t fixedSideDir = controlParameters.getInt("direction");
    // Function on the boundary
    std::string bcFunction = controlParameters.askString("function","u");

    //// Boundary conditions
    gsBoundaryConditions<> bc_u;
    fd_pars.getLabel("BCs_u", bc_u);

    //// Boundary conditions
    gsBoundaryConditions<> bc_d;
    fd_pars.getLabel("BCs_d", bc_d);

    ///////////////////////////////////////////////////////////////////////////////////////
    //PROBLEM SETUP////////////////////////////////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////

    // Construct the basis
    gsMultiBasis<> mb(mp);
    gsInfo<<"The basis has size "<<mb.size()<<" and degree "<<mb.degree()<<"\n";
    for (size_t b=0; b!=mb.nBases(); b++)
        gsInfo<<"Basis "<<b<<":\n"<<mb.basis(b)<<"\n";

    // Boundary conditions
    gsFunctionExpr<> displ(bcFunction,2);
    displ.set_u(ucurr);
    bc_u.addCondition(fixedSideId,condition_type::dirichlet,&displ,fixedSideDir);
    bc_u.setGeoMap(mp);
    bc_d.setGeoMap(mp);

    ///////////////////////////////////////////////////////////////////////////////////////
    //INITIALIZATION///////////////////////////////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////

    // Initialize the solution (deformed geometry)
    gsMultiPatch<> mp_def = mp;
    gsMultiPatch<> displacement = mp;
    gsMultiPatch<> ddisplacement;
    for (index_t p=0; p<mp.nPatches(); ++p)
        displacement.patch(p).coefs().setZero();

    // Initialize the material
    gsLinearDegradedMaterial<real_t> material(E,nu,damage,2);
    // Initialize the elasticity assembler
    gsConstantFunction<> bodyForce(0.,0.,2);
    gsElasticityAssembler<real_t> elAssembler(mp,mb,bc_u,bodyForce,&material);
    std::vector<gsMatrix<> > fixedDofs = elAssembler.allFixedDofs();
    // elAssembler.assemble();
    gsBoundaryConditions<> bc_u_dummy;
    bc_u_dummy.setGeoMap(mp);
    gsElasticityAssembler<real_t> fullElAssembler(mp,mb,bc_u_dummy,bodyForce,&material);
    std::vector<gsMatrix<> > dummyFixedDofs = fullElAssembler.allFixedDofs();

    // Initialize the phase-field assembler
    gsPhaseFieldAssemblerBase<real_t> * pfAssembler;
    if      (order == 2 && AT == 1)
        pfAssembler = new gsPhaseFieldAssembler<real_t,PForder::Second,PFmode::AT1>(mp,mb,bc_d);
    else if (order == 4 && AT == 1)
    {
        pfAssembler = new gsPhaseFieldAssembler<real_t,PForder::Fourth,PFmode::AT1>(mp,mb,bc_d);
        pfAssembler->options().setReal("cw",4.44847);
    }
    else if (order == 2 && AT == 2)
        pfAssembler = new gsPhaseFieldAssembler<real_t,PForder::Second,PFmode::AT2>(mp,mb,bc_d);
    else if (order == 4 && AT == 2)
        pfAssembler = new gsPhaseFieldAssembler<real_t,PForder::Fourth,PFmode::AT2>(mp,mb,bc_d);
    else
        GISMO_ERROR("Invalid order and/or AT model");

    pfAssembler->options().setReal("l0",l0);
    pfAssembler->options().setReal("Gc",Gc);
    pfAssembler->initialize();

    //////////////////////////////////////////////////////////////////////////
    // SOLVE THE PROBLEM/////////////////////////////////////////////////////
    //////////////////////////////////////////////////////////////////////////

    gsMatrix<> u(elAssembler.numDofs(),1), du;
    u.setZero();


#ifdef GISMO_WITH_PARDISO
    gsSparseSolver<>::PardisoLDLT solver;
#else
    gsSparseSolver<>::CGDiagonal solver;
#endif
    gsSparseMatrix<> K;
    gsMatrix<> R;

    real_t elAssemblyTime = 0.0;
    real_t elSolverTime = 0.0;
    real_t pfAssemblyTime = 0.0;
    real_t pfSolverTime = 0.0;

    gsMatrix<> D, deltaD;
    gsSparseMatrix<> Q, QPhi, QPsi;
    gsMatrix<> q, qpsi;
    // Phase-field assembly can already be performed since some terms are independent of the solutions
    pfAssembler->assembleMatrix();
    pfAssembler->matrix_into(QPhi);
    pfAssembler->assembleVector();
    pfAssembler->rhs_into(q);

    index_t step = 0;

    gsParaviewCollection damageCollection("damage");
    gsParaviewCollection psiCollection("Psi");
    gsParaviewCollection displCollection("displacement");
    gsStopwatch clock;

    std::vector<std::vector<real_t>> data;

    // pfAssembler->assembleMatrix();
    // pfAssembler->matrix_into(QPhi);
    // pfAssembler->constructSolution(damage,D);
    // gsInfo<<"D_0 = "<<(0.5 * D.transpose() * QPhi * D).value()<<"\n";
    //     gsWriteParaview(mp,damage,"damage_ini",100000);

    std::ofstream file(outputdir+"results.txt");
    file<<"u,Fx,Fy,E_u,E_d\n";
    file.close();

    while (ucurr<=uend)
    {
        // Update the boundary conditions
        displ.set_u(ucurr);
        elAssembler.computeDirichletDofs(fixedSideDir); // NOTE: This computes the DDofs for **unknown** 1, which should be component 1. This is a bug in the gsElasticity assembler
        fixedDofs = elAssembler.allFixedDofs();
        elAssembler.setFixedDofs(fixedDofs);

        gsInfo<<"Load step "<<step<<": u = "<<ucurr<<"\n";
        gsInfo<<std::setw(20)<<"Iteration"<<std::setw(20)<<"||R||"<<std::setw(20)<<"||dU||/||U||"<<std::setw(20)<<"||dD||/||D||"<<std::setw(20)<<"el. assembly"<<std::setw(20)<<"el. solver"<<std::setw(20)<<"pf. assembly"<<std::setw(20)<<"pf. solver"<<"\n";

        du.setZero();
        deltaD.setZero();
        elAssemblyTime = elSolverTime = 0.0;
        pfAssemblyTime = pfSolverTime = 0.0;

        clock.restart();
        material.setParameter(2,damage);
        elAssembler.homogenizeFixedDofs(-1);
        elAssembler.assemble(u,fixedDofs);
        elAssemblyTime = clock.stop();
        for (index_t it=0; it!=maxIt; ++it)
        {
            // if (it!=0)
            // material.setParameter(2,damage);
            // elAssembler.homogenizeFixedDofs(-1);

            // clock.restart();
            // elAssembler.assemble(u,fixedDofs);
            // elAssemblyTime = clock.stop();
            K = elAssembler.matrix();
            R = elAssembler.rhs();
            clock.restart();
            solver.compute(K);
            du = solver.solve(R);
            elSolverTime = clock.stop();
            u += du;

            elAssembler.setFixedDofs(fixedDofs);
            elAssembler.constructSolution(u,fixedDofs,displacement);
            for (size_t p=0; p!=mp.nPatches(); ++p)
                mp_def.patch(p).coefs() = mp.patch(p).coefs() + displacement.patch(p).coefs();

            // Initialize the function for the elastic energy
            gsMaterialEval<real_t,gsMaterialOutput::Psi> Psi(&material,mp,mp_def);

            // ==================================================================================

            // Phase-field problem
            clock.restart();
            // gsInfo<<"Assembling phase-field problem"<<std::flush;
            pfAssembler->assemblePsiMatrix(Psi);
            pfAssembler->matrix_into(QPsi);
            pfAssembler->assemblePsiVector(Psi);
            pfAssembler->rhs_into(qpsi);
            Q = QPhi + QPsi;
            // gsInfo<<". Done\n";

            pfAssembler->constructSolution(damage,D);
            R = Q * D - qpsi + q;
            pfAssemblyTime = clock.stop();

            clock.restart();
            // solver.compute(Q);
            // deltaD = solver.solve(-R);
            // gsDebugVar(deltaD.norm());
            gsPSOR<real_t> PSORsolver(Q);
            PSORsolver.options().setInt("MaxIterations",30000);
            PSORsolver.options().setSwitch("Verbose",false);
            PSORsolver.options().setReal("tolU",1e-4);
            PSORsolver.options().setReal("tolNeg",1e-6);
            PSORsolver.options().setReal("tolPos",1e-6);
            PSORsolver.solve(R,deltaD); // deltaD = Q \ R
            pfSolverTime = clock.stop();

            D += deltaD;

            // Update damage spline
            pfAssembler->constructSolution(D,damage);

            clock.restart();
            material.setParameter(2,damage);
            elAssembler.homogenizeFixedDofs(-1);
            elAssembler.assemble(u,fixedDofs);
            R = elAssembler.rhs();
            elAssemblyTime = clock.stop();
            // if ((du.norm()/u.norm() < tol || u.norm() < 1e-12 ))//&& deltaD.norm()/D.norm() < tol)

            // Print
            gsInfo<<std::setw(20)<<it<<std::setw(20)<<R.norm()<<std::setw(20)<<du.norm()/u.norm()<<std::setw(20)<<deltaD.norm()/D.norm()<<std::setw(20)<<elAssemblyTime<<std::setw(20)<<elSolverTime<<std::setw(20)<<pfAssemblyTime<<std::setw(20)<<pfSolverTime<<"\n";

            if (R.norm() < 1e-7)
            {
                gsInfo<<"Converged\n";
                break;
            }
        }

        std::string filename;
        filename = "damage_" + util::to_string(step);
        gsWriteParaview(mp,damage,outputdir+filename,100000);
        // gsField<> damage_step(zone,damage,false);
        // gsWriteParaview(damage_step,filename,1000);
        filename += "0";
        damageCollection.addPart(filename,step,"Solution",0);

        gsMaterialEval<real_t,gsMaterialOutput::Psi> Psi(&material,mp,mp_def);
        filename = "Psi_"+util::to_string(step);
        gsWriteParaview(mp,Psi,outputdir+filename,100000);
        filename += "0";
        psiCollection.addPart(filename,step,"Solution",0);

        filename = "displacement_"+util::to_string(step);
        gsWriteParaview(mp,displacement,outputdir+filename,1000);
        filename += "0";
        displCollection.addPart(filename,step,"Solution",0);

        // =========================================================================
        // Compute resulting force and energies
        gsMatrix<> ufull = displacement.patch(0).coefs().reshape(displacement.patch(0).coefs().size(),1);
        fullElAssembler.assemble(ufull,dummyFixedDofs);
        gsMatrix<> Rfull = fullElAssembler.rhs();
        // sum the reaction forces in Y direction
        gsDofMapper mapper(mb,2);
        mapper.finalize();
        gsMatrix<index_t> boundary = mb.basis(0).boundary(boundary::north);
        real_t Fx = 0, Fy = 0;
        for (index_t k=0; k!=boundary.size(); k++)
        {
            Fx += Rfull(mapper.index(boundary(k,0),0,0),0); // DoF index, patch, component
            Fy += Rfull(mapper.index(boundary(k,0),0,1),0); // DoF index, patch, component
        }

        std::vector<real_t> stepData(5);
        stepData[0] = ucurr;
        stepData[1] = Fx;
        stepData[2] = Fy;
        stepData[3] = (0.5 * ufull.transpose() * fullElAssembler.matrix() * ufull).value();
        stepData[4] = (0.5 * D.transpose() * QPhi * D).value() + (D.transpose() * q).value();
        data.push_back(stepData);

        gsInfo<<"Rx = "<<-stepData[1]<<", "
              <<"Ry = "<<-stepData[2]<<", "
              <<"E_u = "<<stepData[3]<<", "
              <<"E_d = "<<stepData[4]<<"\n";

        std::ofstream file(outputdir+"results.txt",std::ios::app);
        // for (size_t i = 0; i != data.size(); ++i)
        //     file<<data[i][0]<<","<<-data[i][1]<<","<<-data[i][2]<<","<<data[i][3]<<","<<data[i][4]<<"\n";
        file<<stepData[0]<<","<<-stepData[1]<<","<<-stepData[2]<<","<<stepData[3]<<","<<stepData[4]<<"\n";
        file.close();

        ucurr += (ucurr+ustep > utrans) ? ustep/ured : ustep;

        step++;
    }

    damageCollection.save();
    psiCollection.save();
    displCollection.save();

    delete pfAssembler;
    return 0;
} // end main

