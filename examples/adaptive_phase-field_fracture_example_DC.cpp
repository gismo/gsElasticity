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

std::vector<real_t> labelElements(  const gsMultiPatch<> & geometry,
                                    const gsFunctionSet<>& damage,
                                    const gsMultiBasis<> & basis,
                                    const real_t         & lowerBound=0.0,
                                    const real_t         & upperBound=1.0);

void refineMesh    (      gsMultiBasis<>      & basis,
                    const std::vector<real_t> & vals,
                    const gsOptionList        & options = gsOptionList());

template <short_t dim, class T>
void solve(gsOptionList & materialParameters,
           gsOptionList & controlParameters,
           gsOptionList & mesherOptions,
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
    std::string parInput;
    std::string geoInput;
    std::string damageInput;
    std::string inputDir;

    gsCmdLine cmd("Tutorial on solving a Linear Elasticity problem.");
    cmd.addInt("e", "numElev","Number of degree elevation steps to perform before solving",numElev);
    cmd.addInt("r", "numHRef","Number of Uniform h-refinement loops", numHRef);
    cmd.addInt("p", "plotmod","Modulo for plotting", plotmod);
    cmd.addSwitch("plot","Create a ParaView visualization file with the solution", plot);
    cmd.addString("o", "output", "Output directory", output);
    cmd.addString("i", "parInput", "Input XML file", parInput);
    cmd.addString("g", "geometry", "Geometry file", geoInput);
    cmd.addString("d", "damage", "Damage file", damageInput);
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
    gsMultiPatch<> mp_ini;
    fd_geo.getFirst(mp_ini);
    if (numElev > 0)
        mp_ini.degreeIncrease(numElev);
    for (index_t i = 0; i<numHRef; ++i)
        mp_ini.uniformRefine(1);
    if (plot) gsWriteParaview(mp_ini,outputdir+"mp",10,true);

    gsFileData<> fd_damage(damageInput.empty() ? inputDir + "damage.xml" : damageInput);
    gsMultiPatch<> damage;
    fd_damage.getFirst(damage);
    if (plot) gsWriteParaview(mp_ini,damage,outputdir+"initial_damage",100000);

    gsFileData<> fd_pars(parInput.empty() ? inputDir + "parameters.xml" : parInput);
    GISMO_ASSERT(fd_pars.hasLabel("material"), "Material parameters not found in the input file.");
    GISMO_ASSERT(fd_pars.hasLabel("control"), "Displacement-Control parameters not found in the input file.");
    GISMO_ASSERT(fd_pars.hasLabel("meshing"), "Adaptive meshing parameters not found in the input file.");
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

    gsOptionList mesherOptions;
    fd_pars.getLabel("meshing", mesherOptions);

    ///////////////////////////////////////////////////////////////////////////////////////
    // Call the dimensional solver
    ///////////////////////////////////////////////////////////////////////////////////////
    switch (mp_ini.domainDim())
    {
        case 2:
            solve<2>(materialParameters,controlParameters,mesherOptions,mp_ini,damage,bc_u,bc_d,plot,plotmod,outputdir);
            break;
        case 3:
            solve<3>(materialParameters,controlParameters,mesherOptions,mp_ini,damage,bc_u,bc_d,plot,plotmod,outputdir);
            break;
        default:
            GISMO_ERROR("Invalid domain dimension");
    }

    return 0;
} // end main

std::vector<real_t> labelElements(  const gsMultiPatch<> & geometry,
                                    const gsFunctionSet<>& damage,
                                    const gsMultiBasis<> & basis,
                                    const real_t         & lowerBound,
                                    const real_t         & upperBound)
{
    GISMO_ASSERT(basis.nBases() == 1, "Labeling is only implemented for single basis meshes");
    gsBasis<real_t>::domainIter domIt  = basis.basis(0).domain()->beginAll();
    gsBasis<real_t>::domainIter domEnd = basis.basis(0).domain()->endAll();
    gsMatrix<real_t,2,2> corners;
    gsMatrix<real_t> points(2,5);
    gsMatrix<real_t> vals;
    std::vector<real_t> labels(basis.basis(0).numElements());
    for (; domIt<domEnd; ++domIt)
    {
        corners.col(0) = domIt.lowerCorner();
        corners.col(1) = domIt.upperCorner();
        gsGridIterator<real_t,CUBE,2> grid(corners,2);
        points.col(0) = domIt.centerPoint();
        points.block(0,1,2,4) = grid.toMatrix();
        damage.piece(0).eval_into(points,vals);
        labels[domIt.id()] = (vals.array() >= lowerBound && vals.array() <= upperBound).any();
    }
    return labels;
}

void refineMesh    (      gsMultiBasis<>      & basis,
                    const std::vector<real_t> & vals,
                    const gsOptionList        & options)
    {
    gsAdaptiveMeshing<2,real_t> mesher(basis);
    mesher.options().setSwitch("Admissible",true);
    mesher.options().setInt("MaxLevel",1);
    mesher.options().setInt("RefineRule",1);
    mesher.options().setInt("CoarsenRule",1);
    mesher.options().setReal("RefineParam",0.1);
    mesher.options().setReal("CoarsenParam",0.1);
    mesher.options().update(options,gsOptionList::ignoreIfUnknown);
    mesher.getOptions();

    gsHBoxContainer<2,real_t> refine, coarsen;
    // Mark the elements for refinement
    mesher.markRef_into(vals,refine);
    // Mesh adaptivity
    mesher.refine(refine);
}

template <short_t dim, class T>
void solve(gsOptionList & materialParameters,
           gsOptionList & controlParameters,
           gsOptionList & mesherOptions,
           gsMultiPatch<T> & mp_ini,
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

    // Convert to THB
    gsMultiPatch<> mp;
    for (index_t i = 0; i < mp_ini.nPatches(); ++i)
    {
        // Check if tensor basis
        if      ((dynamic_cast<const gsTensorBSpline<2,real_t> *>(&mp_ini.patch(i))))
        {
            // Create a THB spline basis
            const gsTensorBSpline<2,real_t> & tb = static_cast<const gsTensorBSpline<2,real_t> &>(mp_ini.patch(i));
            gsTHBSpline<2,real_t> thb(tb);
            mp.addPatch(memory::make_unique(thb.clone().release()));
        }
        else if ((dynamic_cast<const gsTHBSpline<2,real_t> *>(&mp_ini.patch(i))))
        {
            const gsTHBSpline<2,real_t> & thb = static_cast<const gsTHBSpline<2,real_t> &>(mp_ini.patch(i));
            mp.addPatch(memory::make_unique(thb.clone().release()));
        }
        else
            GISMO_ERROR("The basis is not a TB-spline basis or THB-spline basis.");
    }

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
    gsMultiPatch<> mp_def = mp;
    gsMultiPatch<> displacement = mp;
    gsMultiPatch<> ddisplacement;
    for (index_t p=0; p<mp.nPatches(); ++p)
        displacement.patch(p).coefs().setZero();

    // Initialize the material
    gsLinearDegradedMaterial<real_t> material(E,nu,damage,2);
    // Initialize the elasticity assembler
    gsConstantFunction<> bodyForce(0.,0.,2);

    gsBoundaryConditions<> bc_u_dummy;
    bc_u_dummy.setGeoMap(mp);

    // Initialize the phase-field assembler
    gsPhaseFieldAssemblerBase<real_t> * pfAssembler;

    //////////////////////////////////////////////////////////////////////////
    // INITIALIZE THE MESH ///////////////////////////////////////////////////
    //////////////////////////////////////////////////////////////////////////

    std::vector<real_t> elVals;
    // // REFINE MESH
    // for (index_t k=0; k!=mesherOptions.askInt("MaxLevel",1); ++k)
    // {
    //     elVals = labelElements(mp, damage, mb,0.1,1.0);
    //     gsInfo<<"Refine: "<<gsAsVector<real_t>(elVals).sum()<<"\n";
    //     if (gsAsVector<real_t>(elVals).sum())
    //         refineMesh(mb,elVals,mesherOptions);
    //     else
    //         break;
    // }
    // writeSingleCompMesh(mb.basis(0),mp.patch(0),outputdir+"mesh",1);

    //////////////////////////////////////////////////////////////////////////
    // SOLVE THE PROBLEM/////////////////////////////////////////////////////
    //////////////////////////////////////////////////////////////////////////

    gsMatrix<> u, du;

#ifdef GISMO_WITH_PARDISO
    gsSparseSolver<>::PardisoLDLT solver;
#else
    gsSparseSolver<>::CGDiagonal solver;
#endif
    gsSparseMatrix<> K;
    gsMatrix<> R, F;

    T elAssemblyTime = 0.0;
    T elSolverTime = 0.0;
    T pfAssemblyTime = 0.0;
    T pfSolverTime = 0.0;
    T iterationTime  = 0.0;

    gsMatrix<T> D, deltaD;
    gsSparseMatrix<T> Q, QPhi, QPsi;
    gsMatrix<T> q, qpsi;

    index_t step = 0;

    gsParaviewCollection damageCollection(outputdir+"damage");
    gsParaviewCollection psiCollection(outputdir+"Psi");
    gsParaviewCollection displCollection(outputdir+"displacement");
    gsParaviewCollection meshCollection(outputdir+"mesh");
    gsStopwatch smallClock, bigClock;

    std::vector<std::vector<T>> data;

    std::ofstream file(outputdir+"results.txt");
    file<<"u,Fx,Fy,E_u,E_d,basis_size,num_elements,refIt\n";
    file.close();

    std::vector<gsMatrix<> > fixedDofs;
    std::vector<gsMatrix<> > dummyFixedDofs;
    real_t Rnorm, Fnorm;
    Rnorm = Fnorm = 1;
    // gsMultiPatch<> displacement_old, damage_old;
    // displacement_old = displacement;
    // damage_old = damage;
    while (ucurr<=uend)
    {
        bool refined = true;
        index_t basis_size_old, basis_size;
        index_t refIt = 0;
        T basis_size_ratio;
        gsInfo<<"===========================================================================================================================\n";
        gsInfo<<"Load step "<<step<<": u = "<<ucurr<<"\n";
        // Refinement iterations
        for (refIt = 0; refIt!=mesherOptions.askInt("MaxRefIterations",5) && refined; refIt++)
        {
            basis_size = basis_size_old = mb.basis(0).size();
            gsInfo<<"---------------------------------------------------------------------------------------------------------------------------\n";
            gsInfo<<"Refinement iteration "<<refIt<<" (basis size: "<<basis_size<<"):\n";
            // CONSTRUCT ASSEMBLERS (since refreshing is not possible)
            gsElasticityAssembler<real_t> elAssembler(mp,mb,bc_u,bodyForce,&material);
            elAssembler.options().setReal("quA",1.0);
            elAssembler.options().setInt ("quB",0);
            elAssembler.options().setSwitch ("SmallStrains",true);
            fixedDofs = elAssembler.allFixedDofs();
            // Initialize u
            elAssembler.constructSolution(displacement,u);

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

            // Pre-assemble the phase-field operators that do not depend on the solution
            pfAssembler->options().setReal("l0",l0);
            pfAssembler->options().setReal("Gc",Gc);
            pfAssembler->options().setReal("ExprAssembler.quA",1.0);
            pfAssembler->options().setInt ("ExprAssembler.quB",0);
            pfAssembler->options().setInt("ExprAssembler.DirichletValues",dirichlet::l2Projection);
            pfAssembler->setSpaceBasis(mb);
            pfAssembler->initialize();

            // Construct Phase-Field solution vector from projected solution
            pfAssembler->constructSolution(damage,D);

            // Pre-assemble the phase-field operators that do not depend on the solution
            pfAssembler->assembleMatrix();
            pfAssembler->matrix_into(QPhi);
            pfAssembler->assembleVector();
            pfAssembler->rhs_into(q);

            // Update the boundary conditions
            displ.set_u(ucurr);
            elAssembler.computeDirichletDofs(fixedSideDir); // NOTE: This computes the DDofs for **unknown** 1, which should be component 1. This is a bug in the gsElasticity assembler
            fixedDofs = elAssembler.allFixedDofs();
            elAssembler.setFixedDofs(fixedDofs);

            deltaD.setZero(D.rows(),1);

            for (index_t it=0; it!=maxIt; ++it)
            {
                elAssemblyTime = elSolverTime = 0.0;
                pfAssemblyTime = pfSolverTime = 0.0;
                iterationTime  = 0.0;
                bigClock.restart();
                gsInfo<<"    --------------------------Staggered Iteration: "<<PRINT(4)<<it<<"--------------------------\n";
                gsInfo<<"    ---------------------------------ELASTICITY----------------------------------\n";
                gsInfo<<"    | "<<PRINT(6)<<"It."<<PRINT(14)<<"||R||"<<PRINT(14)<<"||dU||/||U||"<<PRINT(20)<<"cum. assembly [s]"<<PRINT(20)<<"cum. solver [s]"<<"|\n";

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
                    gsInfo<<"    | "<<PRINT(6)<<elIt<<PRINT(14)<<Rnorm/Fnorm<<PRINT(14)<<u.norm()<<PRINT(20)<<elAssemblyTime<<PRINT(20)<<elSolverTime<<"|\n";

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
                gsMaterialEval<T,gsMaterialOutput::Psi,true,true> Psi(&material,mp,mp_def);

                // ==================================================================================

                gsInfo<<"    ---------------------------------PHASE-FIELD---------------------------------\n";
                gsInfo<<"    | "<<PRINT(6)<<"It."<<PRINT(14)<<"||R||"<<PRINT(14)<<"||dD||/||D||"<<PRINT(20)<<"cum. assembly [s]"<<PRINT(20)<<"cum. solver [s]"<<"|\n";

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
                PSORsolver.options().setReal("tolNeg",1e-6);
                PSORsolver.options().setReal("tolPos",1e-6);
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

                    gsInfo<<"    | "<<PRINT(6)<<pfIt<<PRINT(14)<<R.norm()<<PRINT(14)<<deltaD.norm()/D.norm()<<PRINT(20)<<pfAssemblyTime<<PRINT(20)<<pfSolverTime<<"|\n";
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
                gsInfo<<"    ----------------------------------FINISHED-----------------------------------\n";
                gsInfo<<"    | "<<PRINT(6)<<"||R|| = "<<PRINT(14)<<Rnorm<<PRINT(20)<<"total time [s] = "<<PRINT(20)<<iterationTime<<"\n";
                gsInfo<<"    | "<<PRINT(6)<<"||F|| = "<<PRINT(14)<<Fnorm<<PRINT(20)<<"elasticity time [s] = "<<PRINT(20)<<elAssemblyTime+elSolverTime<<"\n";
                gsInfo<<"    | "<<PRINT(6)<<"||D|| = "<<PRINT(14)<<D.norm()<<PRINT(20)<<"phase-field time [s] = "<<PRINT(20)<<pfAssemblyTime+pfSolverTime<<"\n";
                gsInfo<<"    -----------------------------------------------------------------------------\n";
                if (Rnorm/* /Fnorm */ < tol)
                    break;
                else if (it == maxIt-1)
                    GISMO_ERROR("Staggered iterations problem did not converge.");
            }

            gsInfo<<"\n";
            gsInfo<<"Converged with ||R||/||F|| = "<<Rnorm/Fnorm<<" < "<<tol<<" ||D|| = "<<D.norm()<<" ||U|| = "<<u.norm()<<"\n";
            // gsInfo<<"----------------------------------------------------------------------------------------------------\n\n";

            // =========================================================================
            // REFINE MESH
            // All labelled elements are refined to the maximum level, step-by-step
            for (index_t i=0; i!=mesherOptions.askInt("MaxLevel",1); ++i)
            {
                elVals = labelElements(mp, damage, mb,0.1,1.0);
                if (gsAsVector<real_t>(elVals).sum() > 0)
                    refineMesh(mb,elVals,mesherOptions);

                basis_size = mb.basis(0).size();
                refined = basis_size > basis_size_old;
                if (!refined)
                    break;
            }
            basis_size_ratio = (T)basis_size/basis_size_old;
            gsInfo<<"Old mesh size: "<<basis_size_old<<", new mesh size: "<<basis_size<<", ratio = "<<basis_size_ratio<<"\n";
            refined &= basis_size_ratio > mesherOptions.askReal("SizeRatio",1.05);

            // =========================================================================
            // PROJECT SOLUTIONS
            gsMatrix<> projCoefs;
            // Geometry
            gsQuasiInterpolate<real_t>::localIntpl(mb.basis(0),mp.patch(0),projCoefs);
            mp.clear();
            mp.addPatch(mb.basis(0).makeGeometry(give(projCoefs)));
            mp_def = mp;
            // Displacement
            gsQuasiInterpolate<real_t>::localIntpl(mb.basis(0),displacement.patch(0),projCoefs);
            displacement.clear();
            displacement.addPatch(mb.basis(0).makeGeometry(give(projCoefs)));
            // gsQuasiInterpolate<real_t>::localIntpl(mb.basis(0),displacement_old.patch(0),projCoefs);
            // displacement_old.clear();
            // displacement_old.addPatch(mb.basis(0).makeGeometry(give(projCoefs)));
            // Damage
            gsQuasiInterpolate<real_t>::localIntpl(mb.basis(0),damage.patch(0),projCoefs);
            damage.clear();
            damage.addPatch(mb.basis(0).makeGeometry(give(projCoefs)));
            // gsQuasiInterpolate<real_t>::localIntpl(mb.basis(0),damage_old.patch(0),projCoefs);
            // damage_old.clear();
            // damage_old.addPatch(mb.basis(0).makeGeometry(give(projCoefs)));
        }

        // =========================================================================
        // PLOT
        if (plot && step%plotmod==0)
        {
            std::string filename;

            filename = "mesh_" + util::to_string(step);
            gsMesh<> mesh(mb.basis(0));
            writeSingleCompMesh(mb.basis(0),mp.patch(0),outputdir+filename,1);
            meshCollection.addPart(filename,step,"Mesh",0);

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
        }

        // =========================================================================

        // =========================================================================
        gsElasticityAssembler<real_t> fullElAssembler(mp,mb,bc_u_dummy,bodyForce,&material);
        fullElAssembler.options().setReal("quA",1.0);
        fullElAssembler.options().setInt ("quB",0);
        fullElAssembler.options().setSwitch ("SmallStrains",true);
        dummyFixedDofs = fullElAssembler.allFixedDofs();

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

        std::vector<real_t> stepData(8);
        stepData[0] = ucurr;
        stepData[1] = Fx;
        stepData[2] = Fy;
        stepData[3] = (0.5 * ufull.transpose() * fullElAssembler.matrix() * ufull).value();
        stepData[4] = (0.5 * D.transpose() * QPhi * D).value() + (D.transpose() * q).value();
        stepData[5] = mb.totalSize();
        stepData[6] = mb.totalElements();
        stepData[7] = refIt;;


        // Write data
        std::ofstream file(outputdir+"results.txt",std::ios::app);
        // for (size_t i = 0; i != data.size(); ++i)
        //     file<<data[i][0]<<","<<-data[i][1]<<","<<-data[i][2]<<","<<data[i][3]<<","<<data[i][4]<<"\n";
        file<<stepData[0]<<","<<-stepData[1]<<","<<-stepData[2]<<","<<stepData[3]<<","<<stepData[4]<<","<<stepData[5]<<","<<stepData[6]<<","<<stepData[7]<<"\n";
        file.close();

        // =========================================================================
        // INCREMENT STEP

        // displacement_old = displacement;
        // damage_old = damage;

        ucurr += (ucurr+ustep > utrans) ? ustep/ured : ustep;
        step++;
    }

    if (plot)
    {
        meshCollection.save();
        damageCollection.save();
        psiCollection.save();
        displCollection.save();
    }


    delete pfAssembler;
}