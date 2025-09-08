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
#include <gsElasticity/gsSolidAssembler.h>
#include <gsElasticity/gsPhaseFieldAssembler.h>
#include <gsElasticity/gsPSOR.h>
#include <gsUtils/gsStopwatch.h>

using namespace gismo;
//! [Include namespace]

#define PRINT(w) std::setw(w)<<std::left

template <class T>
struct times
{
    T elAssemblyTime;
    T elSolverTime;
    T pfAssemblyTime;
    T pfSolverTime;
    T projectionTime;

    void reset()
    {
        elAssemblyTime = 0;
        elSolverTime = 0;
        pfAssemblyTime = 0;
        pfSolverTime = 0;
        projectionTime = 0;
    }
};

template <class T>
typename gsMultiGridOp<T>::uPtr setupMultiGrid(const std::vector< gsSparseMatrix<T,RowMajor> > & transferMatrices,
                                              const gsSparseMatrix<T> & matrix,
                                              const gsOptionList & options)
{
    // Transfer matrices are consumed by the multigrid solver, so we need to make a copy
    std::vector< gsSparseMatrix<T,RowMajor> > myTransferMatrices = transferMatrices;

    // Setup the multigrid solver
    typename gsMultiGridOp<T>::uPtr mg = gsMultiGridOp<T>::make( matrix, myTransferMatrices );
    mg->setOptions( options );
    // Since we are solving a symmetric positive definite problem,we can use a Cholesky solver
    mg->setCoarseSolver( makeSparseCholeskySolver( mg->matrix(0) ) );


    // Parse smoother sequence and validate length
    std::vector<std::string> smoothers;
    std::istringstream ss(options.getString("SmootherSequence"));
    std::string smoother;
    while (std::getline(ss, smoother, ';'))
    {
        smoothers.push_back(smoother);
    }

    if (smoothers.size() < mg->numLevels() - 1)
    {
        gsWarn<<"WARNING: Number of smoothers ("<<smoothers.size()<<") is less than number of levels-1 ("<<(mg->numLevels()-1)<<"). "
                <<"Using Jacobi as fallback for the remaining levels."<<std::endl;
        while (smoothers.size() < mg->numLevels() - 1)
            smoothers.push_back("Jacobi");
    }

    // Setup smoothers with profiling info
    for (index_t i = 1; i < mg->numLevels(); ++i)
    {
        gsPreconditionerOp<>::Ptr smootherOp;

        // Get smoother type from sequence (level i-1 since we start from level 1)
        std::string smootherType = smoothers[i-1];

        if      (smootherType == "j" ||
                 smootherType == "J" ||
                 smootherType == "jacobi" ||
                 smootherType == "Jacobi")
        {
            // Jacobi smoother with damping
            smootherOp = makeJacobiOp(mg->matrix(i), options.askReal("JacobiDamping",0.8));
            // gsInfo << "Level " << i << ": Jacobi smoother (damping=" << options.askReal("JacobiDamping",0.8) << ")" << std::endl;
        }
        else if (smootherType == "gs" ||
                 smootherType == "G" ||
                 smootherType == "Gauss-Seidel" ||
                 smootherType == "gauss-seidel" ||
                 smootherType == "GaussSeidel" ||
                 smootherType == "gaussseidel")
        {
            // Symmetric Gauss-Seidel smoother
            smootherOp = makeSymmetricGaussSeidelOp(mg->matrix(i));
            // gsInfo << "Level " << i << ": Symmetric Gauss-Seidel smoother" << std::endl;
        }
        else
        {
            gsInfo << "WARNING: Unknown smoother type '" << smootherType << "' at level " << i
                << ". Using Jacobi as fallback." << std::endl;
            smootherOp = makeJacobiOp(mg->matrix(i), options.askReal("JacobiDamping",0.8));
        }

        smootherOp->setOptions(options);
        mg->setSmoother(i, smootherOp);

        // gsInfo << "Level " << i << ": " << mg->matrix(i).rows() << "x" << mg->matrix(i).cols()
            // << " matrix (" << mg->matrix(i).nonZeros() << " nnz)" << std::endl;
    }
    return mg;
}

template <short_t dim, class T>
void solve(gsOptionList & materialParameters,
           gsOptionList & controlParameters,
           gsOptionList & solverParameters,
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

    //// Solver parameters
    gsOptionList solverParameters;
    if (fd_pars.hasLabel("solver"))
        fd_pars.getLabel("solver", solverParameters);

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
            solve<2>(materialParameters,controlParameters,solverParameters,mp,damage,bc_u,bc_d,plot,plotmod,outputdir);
            break;
        case 3:
            solve<3>(materialParameters,controlParameters,solverParameters,mp,damage,bc_u,bc_d,plot,plotmod,outputdir);
            break;
        default:
            GISMO_ERROR("Invalid domain dimension");
    }

    return 0;
} // end main

template <short_t dim, class T>
void solve(gsOptionList & materialParameters,
           gsOptionList & controlParameters,
           gsOptionList & solverParameters,
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
    bc_u.addCondition(fixedSideId,condition_type::dirichlet,&displ,0,false,fixedSideDir);
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
    gsSolidAssembler<dim,T,gsLinearDegradedMaterial<T>> elAssembler(mp,mb,bc_u,&material);
    elAssembler.options().setReal("ExprAssembler.quA",1.0);
    elAssembler.options().setInt ("ExprAssembler.quB",0);
    elAssembler.initialize();
    gsBoundaryConditions<T> bc_u_dummy;
    bc_u_dummy.setGeoMap(mp);
    gsSolidAssembler<dim,T,gsLinearDegradedMaterial<T>> fullElAssembler(mp,mb,bc_u_dummy,&material);
    fullElAssembler.options().setReal("ExprAssembler.quA",1.0);
    fullElAssembler.options().setInt ("ExprAssembler.quB",0);
    fullElAssembler.initialize();

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

// #ifdef gsMUMPS_ENABLED
//     // Initialize MUMPS solver
//     gsEigen::MUMPSLDLT<gsSparseMatrix<T>,gsEigen::Lower> solver;
// #elif GISMO_WITH_PARDISO
//     typename gsSparseSolver<T>::PardisoLDLT solver;
// #else
    // typename gsSparseSolver<T>::CGDiagonal solver;
// #endif

    // Setup multigrid hierarchy (only once)
    std::vector< gsSparseMatrix<T,RowMajor> > transferMatrices;
    gsInfo<<solverParameters<<"\n";
    if (solverParameters.askSwitch("MultiGrid",false) && solverParameters.hasGroup("MG"))
    {
        gsGridHierarchy<>::buildByCoarsening(mb,dim,bc_u,solverParameters.getGroup("MG")).moveTransferMatricesTo(transferMatrices);
        gsInfo<<"Using Multi-Grid solver with hierarchy:\n";
        for (size_t i = transferMatrices.size(); i!= 0; i--)
            gsInfo << "Level " << i << ": " << transferMatrices[i-1].rows() << " -> " << transferMatrices[i-1].cols() << "\n";
    }

    times<T> stagTimes;
    times<T> stepTimes;
    stagTimes.reset();
    stepTimes.reset();
    gsStopwatch smallClock, bigClock;

    gsSparseMatrix<T> elMatrix;
    gsMatrix<T> elRhs;

    gsMatrix<T> R;
    gsMatrix<T> D, deltaD;
    gsSparseMatrix<T> Q, QPhi, QPsi;
    gsMatrix<T> q, qpsi;
    // Phase-field assembly can already be performed since some terms are independent of the solutions
    smallClock.restart();
    pfAssembler->assemblePhi();
    pfAssembler->matrix_into(QPhi);
    pfAssembler->rhs_into(q);
    T pfAssemblyTime0 = smallClock.stop();
    index_t step = 0;

    gsParaviewCollection damageCollection(outputdir+"damage");
    gsParaviewCollection psiCollection(outputdir+"Psi");
    gsParaviewCollection displCollection(outputdir+"displacement");

    std::ofstream file;
    file.open(outputdir+"results.txt");
    file<<"LoadStep,u,Fx,Fy,E_u,E_d,elAssemblyTime,elSolverTime,pfAssemblyTime,pfSolverTime,basis_size,totIt_el,totIt_pf,numIt_stag\n";
    file.close();
    file.open(outputdir+"iteration_results.txt");
    file<<"LoadStep,StagIt,u,Unorm,Dnorm,Rnorm,Fnorm,relRnorm,elAssemblyTime,elSolverTime,pfAssemblyTime,pfSolverTime,basis_size,numIt_el,numIt_pf\n";
    file.close();

    T Rnorm, Fnorm;
    Rnorm = Fnorm = 1;
    index_t numIt_el = 0, totIt_el = 0; //numIt: per staggered iteration, totIt: total number of iterations across all staggered iterations
    index_t numIt_pf = 0, totIt_pf = 0; //numIt: per staggered iteration, totIt: total number of iterations across all staggered iterations
    index_t numIt_stag = 0; // number of staggered iterations per load step
    index_t basis_size = mb.basis(0).size();
    while (ucurr<=uend)
    {
        numIt_stag = 0;
        totIt_el = totIt_pf = 0;
        stepTimes.reset();
        stagTimes.reset();

        // Update the boundary conditions
        displ.set_u(ucurr);
        elAssembler.initialize();

        gsInfo<<"---------------------------------------------------------------------------------------------------------------------------\n";
        gsInfo<<"Load step "<<step<<": u = "<<ucurr<<"\n\n";

        deltaD.setZero();
        index_t stagIt = 0;
        while(true)
        {
            stagTimes.reset();
            if (step==0 && stagIt==0) stagTimes.pfAssemblyTime = pfAssemblyTime0; // save the time of the first assembly

            bigClock.restart();
            gsInfo<<" - Staggered iteration "<<stagIt<<":\n";
            gsInfo<<"\t"<<PRINT(20)<<"* Elasticity:"<<PRINT(6)<<"It."<<PRINT(14)<<"||R||"<<PRINT(14)<<"||F||"<<PRINT(14)<<"||R||/||F||"<<PRINT(14)<<"||U||"<<PRINT(20)<<"cum. assembly [s]"<<PRINT(20)<<"cum. solver [s]"<<PRINT(20)<<"it. solver [s]"<<PRINT(20)<<"solver it."<<"\n";

            material.setParameter(2,damage);
            // Pre-assemble the elasticity problem
            smallClock.restart();
            elAssembler.assemble(u);
            stagTimes.elAssemblyTime += smallClock.stop();
            elAssembler.matrix_into(elMatrix);
            elAssembler.rhs_into(elRhs);
            Fnorm = elRhs.norm();
            Fnorm = (Fnorm == 0) ? 1 : Fnorm;
            index_t elIt = 0;
            while(true)
            {
                // Solve
                T itSolverTime = 0;
                index_t itSolverIterations = 0;
                smallClock.restart();
                if (solverParameters.askSwitch("MultiGrid",false) && solverParameters.hasGroup("MG"))
                {
                    typename gsSparseSolver<T>::CGCustom solver;
                    solver.preconditioner().set(setupMultiGrid<T>(transferMatrices,elMatrix,solverParameters.getGroup("MG")));
                    solver.compute(elMatrix);
                    u = solver.solveWithGuess(elRhs,u);
                    itSolverIterations = solver.iterations();
                }
                else
                {
                    typename gsSparseSolver<T>::CGDiagonal solver;
                    solver.compute(elMatrix);
                    u = solver.solve(elRhs);
                    itSolverIterations = solver.iterations();
                }
                itSolverTime = smallClock.stop();
                stagTimes.elSolverTime += itSolverTime;

                smallClock.restart();
                elAssembler.assemble(u);
                stagTimes.elAssemblyTime += smallClock.stop();
                elAssembler.matrix_into(elMatrix);
                elAssembler.rhs_into(elRhs);
                Fnorm = elRhs.norm();
                Fnorm = (Fnorm == 0) ? 1 : Fnorm;

                // Check convergence with the old matrix and rhs (saves one assembly)
                Rnorm = (elMatrix*u - elRhs).norm();
                gsInfo<<"\t"<<PRINT(20)<<""<<PRINT(6)<<elIt<<PRINT(14)<<Rnorm<<PRINT(14)<<Fnorm<<PRINT(14)<<Rnorm/Fnorm<<PRINT(14)<<u.norm()<<PRINT(20)<<stagTimes.elAssemblyTime<<PRINT(20)<<stagTimes.elSolverTime<<PRINT(20)<<itSolverTime<<PRINT(20)<<itSolverIterations<<"\n";

                if (Rnorm/Fnorm < tolEl || u.norm() < 1e-12 || maxItEl==1)
                    break;
                else if (elIt == maxItEl-1)
                    GISMO_ERROR("Elasticity problem did not converge.");
                else
                    elIt++;
            }
            numIt_el = elIt+1;
            totIt_el+= numIt_el;

            elAssembler.constructSolution(u,displacement);
            for (size_t p=0; p!=mp.nPatches(); ++p)
                mp_def.patch(p).coefs() = mp.patch(p).coefs() + displacement.patch(p).coefs();

            // Initialize the function for the elastic energy
            gsMaterialEval<T,gsMaterialOutput::Psi,true,true> Psi(&material,mp,mp_def);

            // ==================================================================================

            gsInfo<<"\t"<<PRINT(20)<<"* Phase-field:"<<PRINT(6)<<"It."<<PRINT(14)<<"||R||"<<PRINT(14)<<"||D||"<<PRINT(14)<<"||dD||/||D||"<<PRINT(20)<<"cum. assembly [s]"<<PRINT(20)<<"cum. solver [s]"<<"\n";

            // Phase-field problem
            // gsInfo<<"Assembling phase-field problem"<<"\n";
            smallClock.restart();
            pfAssembler->assemblePsi(Psi);
            stagTimes.pfAssemblyTime = smallClock.stop();
            pfAssembler->matrix_into(QPsi);
            pfAssembler->rhs_into(qpsi);
            if (qpsi.rows()==0) // qpsi is empty for AT2 models
                qpsi = gsMatrix<T>::Zero(QPsi.rows(),1);
            Q = QPhi + QPsi;

            // Reconstruct the solution from the damage field
            pfAssembler->constructSolution(damage,D);
            // gsInfo<<". Done\n";

            // Initialize the PSOR solver
            smallClock.restart();
            gsPSOR<T> PSORsolver(Q);
            PSORsolver.options().setInt("MaxIterations",30000);
            PSORsolver.options().setSwitch("Verbose",false);
            PSORsolver.options().setReal("tolU",1e-4);
            PSORsolver.options().setReal("tolNeg",1e-9);
            PSORsolver.options().setReal("tolPos",1e-9);
            stagTimes.pfSolverTime = smallClock.stop();
            index_t pfIt = 0;
            while(true)
            {
                // Assemble
                smallClock.restart();
                R = Q * D - qpsi + q;
                stagTimes.pfAssemblyTime += smallClock.stop();

                // Solve
                // solver.compute(Q);
                // deltaD = solver.solve(-R);
                // gsDebugVar(deltaD.norm());
                smallClock.restart();
                PSORsolver.solve(R,deltaD); // deltaD = Q \ R
                stagTimes.pfSolverTime += smallClock.stop();
                D += deltaD;
                gsInfo<<"\t"<<PRINT(20)<<""<<PRINT(6)<<pfIt<<PRINT(14)<<R.norm()<<PRINT(14)<<D.norm()<<PRINT(14)<<deltaD.norm()/D.norm()<<PRINT(20)<<stagTimes.pfAssemblyTime<<PRINT(20)<<stagTimes.pfSolverTime<<"\n";;
                if (deltaD.norm()/D.norm() < tolPf || D.norm() < 1e-12 || maxItPf == 1)
                    break;
                else if (pfIt == maxItPf-1)
                    GISMO_ERROR("Phase-field problem did not converge.");
                else
                    pfIt++;
            }
            numIt_pf = pfIt+1;
            totIt_pf+= numIt_pf;

            // Update damage spline
            pfAssembler->constructSolution(D,damage);

            material.setParameter(2,damage);
            smallClock.restart();
            elAssembler.assemble(u);
            stagTimes.elAssemblyTime += smallClock.stop();
            Fnorm = elAssembler.rhs().norm();
            Fnorm = (Fnorm == 0) ? 1 : Fnorm;
            Rnorm = (elAssembler.matrix()*u - elAssembler.rhs()).norm();

            gsInfo<<"\t"<<PRINT(20)<<"* Finished"<<PRINT(6)<<""<<PRINT(14)<<"||R||"<<PRINT(14)<<"||R||/||F||"<<PRINT(14)<<"total [s]"<<PRINT(20)<<"elasticity [s]"           <<PRINT(20)<<"phase-field [s]"          <<"\n";
            gsInfo<<"\t"<<PRINT(20)<<""          <<PRINT(6)<<""<<PRINT(14)<<Rnorm<<PRINT(14)<<Rnorm/Fnorm<<PRINT(14)<<bigClock.stop()<<PRINT(20)<<stagTimes.elAssemblyTime+stagTimes.elSolverTime<<PRINT(20)<<stagTimes.pfAssemblyTime+stagTimes.pfSolverTime<<"\n";

            file.open(outputdir+"iteration_results.txt",std::ios::app);
            file<<step<<","<<stagIt<<","<<ucurr<<","
                <<u.norm()<<","<<D.norm()<<","<<Rnorm<<","<<Fnorm<<","<<Rnorm/Fnorm<<","
                <<stagTimes.elAssemblyTime<<","<<stagTimes.elSolverTime<<","
                <<stagTimes.pfAssemblyTime<<","<<stagTimes.pfSolverTime<<","
                <<basis_size<<","
                <<numIt_el<<","<<numIt_pf<<"\n";
            file.close();

            stepTimes.elAssemblyTime += stagTimes.elAssemblyTime;
            stepTimes.elSolverTime += stagTimes.elSolverTime;
            stepTimes.pfAssemblyTime += stagTimes.pfAssemblyTime;
            stepTimes.pfSolverTime += stagTimes.pfSolverTime;

            if (Rnorm/Fnorm < tol)
                break;
            else if (stagIt == maxIt-1)
                GISMO_ERROR("Staggered iterations problem did not converge.");
            else
                stagIt++;
        }
        numIt_stag += stagIt+1;

        // =========================================================================
        // Compute resulting force and energies
        gsMatrix<T> ufull = displacement.patch(0).coefs().reshape(displacement.patch(0).coefs().size(),1);
        fullElAssembler.assemble(ufull);
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

        gsInfo<<"\n";
        gsInfo<<"Converged with ||R||/||F|| = "<<Rnorm/Fnorm<<" < "<<tol<<" ||D|| = "<<D.norm()<<" ||U|| = "<<u.norm()<<"\n";

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

            gsMaterialEval<T,gsMaterialOutput::Psi,true,true> Psi(&material,mp,mp_def);
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
        file.open(outputdir+"results.txt",std::ios::app);
        file<<step<<","<<ucurr<<","<<-Fx<<","<<-Fy<<","
            <<(0.5 * ufull.transpose() * fullElAssembler.matrix() * ufull).value()<<","
            <<(0.5 * D.transpose() * QPhi * D).value() + (D.transpose() * q).value()<<","
            <<stepTimes.elAssemblyTime<<","<<stepTimes.elSolverTime<<","
            <<stepTimes.pfAssemblyTime<<","<<stepTimes.pfSolverTime<<","
            <<basis_size<<","
            <<totIt_el<<","<<totIt_pf<<","
            <<numIt_stag<<"\n";
        file.close();

        ucurr += (ucurr+ustep > utrans) ? ustep/ured : ustep;
        ucurr = math::min(ucurr,uend);
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

