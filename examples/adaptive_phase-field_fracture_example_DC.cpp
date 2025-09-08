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
#include <gsHSplines/gsHElementHelper.h>
#include <gsHSplines/gsHElementMarker.h>

using namespace gismo;
//! [Include namespace]

#define PRINT(w) std::setw(w)<<std::left

std::vector<real_t> labelElements(  const gsMultiPatch<> & geometry,
                                    const gsFunctionSet<>& damage,
                                    const gsMultiBasis<> & basis,
                                    const real_t         & lowerBound=0.0,
                                    const real_t         & upperBound=1.0);

template <short_t dim,class T>
T refineMesh(      gsMultiBasis<T> & basis,
             const std::vector<T>  & vals,
             const gsOptionList    & options = gsOptionList());

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

    gsOptionList mesherOptions;
    fd_pars.getLabel("meshing", mesherOptions);

    ///////////////////////////////////////////////////////////////////////////////////////
    // Call the dimensional solver
    ///////////////////////////////////////////////////////////////////////////////////////
    switch (mp_ini.domainDim())
    {
        case 2:
            solve<2>(materialParameters,controlParameters,solverParameters,mesherOptions,mp_ini,damage,bc_u,bc_d,plot,plotmod,outputdir);
            break;
        case 3:
            solve<3>(materialParameters,controlParameters,solverParameters,mesherOptions,mp_ini,damage,bc_u,bc_d,plot,plotmod,outputdir);
            break;
        default:
            GISMO_ERROR("Invalid domain dimension");
    }

    return 0;
} // end main

template<short_t dim, class T>
std::vector<T> labelElements(  const gsMultiPatch<> & geometry,
                                    const gsFunctionSet<>& damage,
                                    const gsMultiBasis<> & basis,
                                    const T         & lowerBound,
                                    const T         & upperBound)
{
    GISMO_ASSERT(basis.nBases() == 1, "Labeling is only implemented for single basis meshes");
    // typename gsBasis<T>::domainIter domIt  = basis.basis(0).domain()->beginAll();
    typename gsBasis<T>::domainIter domEnd = basis.basis(0).domain()->endAll();
    std::vector<T> labels(basis.basis(0).numElements());
    gsVector<index_t,dim> np;
    np.setConstant(2);
    gsLobattoRule<T> rule(np); // equivalent to using gsPointGrid with np = 2
// #pragma omp parallel
//     {
#pragma omp parallel for
    for (typename gsBasis<T>::domainIter domIt  = basis.basis(0).domain()->beginAll(); domIt<domEnd; ++domIt)
    {
        gsMatrix<T> nodes;
        gsMatrix<T> vals;
        gsVector<T> weights;

        rule.mapTo(domIt.lowerCorner(), domIt.upperCorner(),nodes,weights);
        damage.piece(0).eval_into(nodes,vals);
        labels[domIt.id()] = (vals.array() >= lowerBound && vals.array() <= upperBound).any();
    }

    // }
    return labels;
}

template <short_t dim, class T>
T refineMesh(         gsMultiBasis<T>& basis,
                const std::vector<T> & vals,
                const gsOptionList   & options)
{
    typedef typename gsHElementHelper<dim,T>::HElementContainer HElementContainer;

    gsHElementMarker<dim,T> marker(basis.basis(0));
    marker.options().setSwitch("Admissible",true);
    marker.options().setInt("MaxLevel",1);
    marker.options().setInt("RefineRule",1);
    marker.options().setInt("CoarsenRule",1);
    marker.options().setReal("RefineParam",0.1);
    marker.options().setReal("CoarsenParam",0.1);
    marker.options().update(options,gsOptionList::ignoreIfUnknown);

    marker.setErrors(vals);
    HElementContainer markedRef = marker.markRef();

    T area = 0.0;
    gsMatrix<T> box;
    for (const auto & elem : markedRef)
    {
        box = marker.helper().toBox(elem);
        area += (box.col(1)-box.col(0)).prod();
    }

    std::vector<index_t> refBox = marker.toRefBoxes(markedRef);
    basis.basis(0).refineElements(refBox);
    return area;
}

template <short_t dim, class T>
void solve(gsOptionList & materialParameters,
           gsOptionList & controlParameters,
           gsOptionList & solverParameters,
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
    gsMultiPatch<T> mp;
    for (index_t i = 0; i < mp_ini.nPatches(); ++i)
    {
        // Check if tensor basis
        if      ((dynamic_cast<const gsTensorBSpline<dim,T> *>(&mp_ini.patch(i))))
        {
            // Create a THB spline basis
            const gsTensorBSpline<dim,T> & tb = static_cast<const gsTensorBSpline<dim,T> &>(mp_ini.patch(i));
            gsTHBSpline<dim,T> thb(tb);
            mp.addPatch(memory::make_unique(thb.clone().release()));
        }
        else if ((dynamic_cast<const gsTHBSpline<dim,T> *>(&mp_ini.patch(i))))
        {
            const gsTHBSpline<dim,T> & thb = static_cast<const gsTHBSpline<dim,T> &>(mp_ini.patch(i));
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
    bc_u.addCondition(fixedSideId,condition_type::dirichlet,&displ,0,false,fixedSideDir);
    bc_u.setGeoMap(mp);
    bc_d.setGeoMap(mp);

    ///////////////////////////////////////////////////////////////////////////////////////
    //INITIALIZATION///////////////////////////////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////

    // Initialize the solution (deformed geometry)
    gsMultiPatch<T> mp_def = mp;
    gsMultiPatch<T> displacement = mp;
    gsMultiPatch<T> ddisplacement;
    for (index_t p=0; p<mp.nPatches(); ++p)
        displacement.patch(p).coefs().setZero();

    // Initialize the material
    gsLinearDegradedMaterial<T> material(E,nu,damage,dim);
    // Initialize the elasticity assembler
    gsVector<T> bodyForceVec(dim);
    bodyForceVec.setZero();
    gsConstantFunction<T> bodyForce(bodyForceVec,dim);

    gsBoundaryConditions<T> bc_u_dummy;
    bc_u_dummy.setGeoMap(mp);

    // Initialize the phase-field assembler
    gsPhaseFieldAssemblerBase<T> * pfAssembler;

    //////////////////////////////////////////////////////////////////////////
    // INITIALIZE THE MESH ///////////////////////////////////////////////////
    //////////////////////////////////////////////////////////////////////////

    std::vector<T> elVals;
    // // REFINE MESH
    // for (index_t k=0; k!=mesherOptions.askInt("MaxLevel",1); ++k)
    // {
    //     elVals = labelElements(mp, damage, mb,0.1,1.0);
    //     gsInfo<<"Refine: "<<gsAsVector<T>(elVals).sum()<<"\n";
    //     if (gsAsVector<T>(elVals).sum())
    //         refineMesh(mb,elVals,mesherOptions);
    //     else
    //         break;
    // }
    // writeSingleCompMesh(mb.basis(0),mp.patch(0),outputdir+"mesh",1);

    //////////////////////////////////////////////////////////////////////////
    // SOLVE THE PROBLEM/////////////////////////////////////////////////////
    //////////////////////////////////////////////////////////////////////////

    gsMatrix<T> u, du;

#ifdef GISMO_WITH_PARDISO
    typename gsSparseSolver<T>::PardisoLDLT solver;
#else
    typename gsSparseSolver<T>::CGDiagonal solver;
#endif

    times<T> stagTimes;
    times<T> stepTimes;
    stagTimes.reset();
    stepTimes.reset();

    gsSparseMatrix<T> elMatrix;
    gsMatrix<T> elRhs;

    gsMatrix<T> R;
    gsMatrix<T> D, deltaD;
    gsSparseMatrix<T> Q, QPhi, QPsi;
    gsMatrix<T> q, qpsi;

    gsParaviewCollection damageCollection(outputdir+"damage");
    gsParaviewCollection psiCollection(outputdir+"Psi");
    gsParaviewCollection displCollection(outputdir+"displacement");
    gsParaviewCollection meshCollection(outputdir+"mesh");
    gsStopwatch smallClock, bigClock;

    std::vector<std::vector<T>> data;

    std::ofstream file;
    file.open(outputdir+"results.txt");
    file<<"LoadStep,u,Fx,Fy,E_u,E_d,elAssemblyTime,elSolverTime,pfAssemblyTime,pfSolverTime,projectionTime,basis_size,ref_area,totIt_el,totIt_pf,numIt_stag,numIt_ref\n";
    file.close();
    file.open(outputdir+"iteration_results.txt");
    file<<"LoadStep,RefIt,StagIt,u,Unorm,Dnorm,Rnorm,Fnorm,relRnorm,elAssemblyTime,elSolverTime,pfAssemblyTime,pfSolverTime,basis_size,ref_area,numIt_el,numIt_pf\n";
    file.close();

    T Rnorm, Fnorm;
    Rnorm = Fnorm = 1;
    // gsMultiPatch<T> displacement_old, damage_old;
    // displacement_old = displacement;
    // damage_old = damage;
    index_t numIt_el = 0, totIt_el = 0; //numIt: per staggered iteration, totIt: total number of iterations across all staggered iterations
    index_t numIt_pf = 0, totIt_pf = 0; //numIt: per staggered iteration, totIt: total number of iterations across all staggered iterations
    index_t numIt_ref = 0; // number of refinement iterations per load step
    index_t numIt_stag = 0; // number of staggered iterations per load step

    index_t step = 0;
    while (ucurr<=uend)
    {
        numIt_ref = numIt_stag = 0;
        totIt_el = totIt_pf = 0;
        stepTimes.reset();
        stagTimes.reset();

        bool refined = true;
        index_t basis_size_old, basis_size;
        T basis_size_ratio;
        T markedArea = 0., tmpArea = 0.;
        gsInfo<<"===========================================================================================================================\n";
        gsInfo<<"Load step "<<step<<": u = "<<ucurr<<"\n";
        // Refinement iterations
        index_t refIt = 0;
        while(true)
        // for (refIt = 0; refIt!=mesherOptions.askInt("MaxRefIterations",5) && refined; refIt++)
        {
            basis_size = basis_size_old = mb.basis(0).size();
            gsInfo<<"---------------------------------------------------------------------------------------------------------------------------\n";
            gsInfo<<"Refinement iteration "<<refIt<<" (basis size: "<<basis_size<<"):\n";
            // CONSTRUCT ASSEMBLERS (since refreshing is not possible)
            gsSolidAssembler<dim,T,gsLinearDegradedMaterial<T>> elAssembler(mp,mb,bc_u,&material);
            elAssembler.options().setReal("ExprAssembler.quA",1.0);
            elAssembler.options().setInt ("ExprAssembler.quB",0);
            elAssembler.options().setInt("ExprAssembler.DirichletValues",dirichlet::l2Projection);
            elAssembler.initialize();
            // Initialize u
            elAssembler.constructSolution(displacement,u);

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
            smallClock.restart();
            pfAssembler->assemblePhi();
            pfAssembler->matrix_into(QPhi);
            pfAssembler->rhs_into(q);
            T pfAssemblyTime0 = smallClock.stop();

            // Update the boundary conditions
            displ.set_u(ucurr);
            elAssembler.initialize();

            deltaD.setZero(D.rows(),1);
            index_t stagIt = 0;

            // Setup multigrid hierarchy (everytime the mesh changes)
            std::vector< gsSparseMatrix<T,RowMajor> > transferMatrices;
            gsInfo<<solverParameters<<"\n";
            if (solverParameters.askSwitch("MultiGrid",false) && solverParameters.hasGroup("MG"))
            {
                gsGridHierarchy<>::buildByHierarchicalCoarsening(mb,dim,bc_u,solverParameters.getGroup("MG")).moveTransferMatricesTo(transferMatrices);
                gsInfo<<"Using Multi-Grid solver with hierarchy:\n";
                for (size_t i = transferMatrices.size(); i!= 0; i--)
                    gsInfo << "Level " << i << ": " << transferMatrices[i-1].rows() << " -> " << transferMatrices[i-1].cols() << "\n";
            }


            while(true)
            {
                stagTimes.reset();
                if (stagIt==0) stagTimes.pfAssemblyTime = pfAssemblyTime0; // save the time of the first assembly

                bigClock.restart();
                gsInfo<<"    --------------------------Staggered Iteration: "<<PRINT(4)<<stagIt<<"--------------------------\n";
                gsInfo<<"    ---------------------------------ELASTICITY----------------------------------\n";
                gsInfo<<"    | "<<PRINT(6)<<"It."<<PRINT(14)<<"||R||"<<PRINT(14)<<"||F||"<<PRINT(14)<<"||R||/||F||"<<PRINT(14)<<"||U||"<<PRINT(20)<<"cum. assembly [s]"<<PRINT(20)<<"cum. solver [s]"<<PRINT(20)<<"it. solver [s]"<<PRINT(20)<<"solver it."<<"|\n";

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
                        solver.setMaxIterations(solverParameters.askInt("MaxIterations",100));
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

                    // Re-assemble
                    smallClock.restart();
                    elAssembler.assemble(u);
                    stagTimes.elAssemblyTime += smallClock.stop();
                    elAssembler.matrix_into(elMatrix);
                    elAssembler.rhs_into(elRhs);
                    Fnorm = elRhs.norm();
                    Fnorm = (Fnorm == 0) ? 1 : Fnorm;

                    // Check convergence with the old matrix and rhs (saves one assembly)
                    Rnorm = (elMatrix*u - elRhs).norm();
                    gsInfo<<"    | "<<PRINT(6)<<elIt<<PRINT(14)<<Rnorm<<PRINT(14)<<Fnorm<<PRINT(14)<<Rnorm/Fnorm<<PRINT(14)<<u.norm()<<PRINT(20)<<stagTimes.elAssemblyTime<<PRINT(20)<<stagTimes.elSolverTime<<PRINT(20)<<itSolverTime<<PRINT(20)<<itSolverIterations<<"|\n";

                    gsInfo<<"Terminating simulation for testing purposes\n";
                    return;

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

                gsInfo<<"    ---------------------------------PHASE-FIELD---------------------------------\n";
                gsInfo<<"    | "<<PRINT(6)<<"It."<<PRINT(14)<<"||R||"<<PRINT(14)<<"||dD||/||D||"<<PRINT(20)<<"cum. assembly [s]"<<PRINT(20)<<"cum. solver [s]"<<"|\n";

                // Phase-field problem
                // gsInfo<<"Assembling phase-field problem"<<"\n";
                smallClock.restart();
                pfAssembler->assemblePsi(Psi);
                stagTimes.pfAssemblyTime += smallClock.stop();
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
                    gsInfo<<"    | "<<PRINT(6)<<pfIt<<PRINT(14)<<R.norm()<<PRINT(14)<<deltaD.norm()/D.norm()<<PRINT(20)<<stagTimes.pfAssemblyTime<<PRINT(20)<<stagTimes.pfSolverTime<<"|\n";
                    if (deltaD.norm()/D.norm() < tolPf || D.norm() < 1e-12 || maxItPf==1)
                        break;
                    else if (pfIt == maxItPf-1 && maxItPf != 1)
                        GISMO_ERROR("Phase-field problem did not converge.");
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

                gsInfo<<"    ----------------------------------FINISHED-----------------------------------\n";
                gsInfo<<"    | "<<PRINT(6)<<"||R|| = "<<PRINT(14)<<Rnorm<<PRINT(20)<<"total time [s] = "<<PRINT(20)<<bigClock.stop()<<"\n";
                gsInfo<<"    | "<<PRINT(6)<<"||F|| = "<<PRINT(14)<<Fnorm<<PRINT(20)<<"elasticity time [s] = "<<PRINT(20)<<stagTimes.elAssemblyTime+stagTimes.elSolverTime<<"\n";
                gsInfo<<"    | "<<PRINT(6)<<"||D|| = "<<PRINT(14)<<D.norm()<<PRINT(20)<<"phase-field time [s] = "<<PRINT(20)<<stagTimes.pfAssemblyTime+stagTimes.pfSolverTime<<"\n";
                gsInfo<<"    -----------------------------------------------------------------------------\n";

                file.open(outputdir+"iteration_results.txt",std::ios::app);
                file<<step<<","<<refIt<<","<<stagIt<<","<<ucurr<<","
                    <<u.norm()<<","<<D.norm()<<","<<Rnorm<<","<<Fnorm<<","<<Rnorm/Fnorm<<","
                    <<stagTimes.elAssemblyTime<<","<<stagTimes.elSolverTime<<","
                    <<stagTimes.pfAssemblyTime<<","<<stagTimes.pfSolverTime<<","
                    <<basis_size<<","<<markedArea<<","
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

            gsInfo<<"\n";
            gsInfo<<"Converged with ||R||/||F|| = "<<Rnorm/Fnorm<<" < "<<tol<<" ||D|| = "<<D.norm()<<" ||U|| = "<<u.norm()<<"\n";
            // gsInfo<<"----------------------------------------------------------------------------------------------------\n\n";

            // =========================================================================
            // REFINE MESH
            // All labelled elements are refined to the maximum level, step-by-step
            for (index_t i=0; i!=mesherOptions.askInt("MaxLevel",1); ++i)
            {
                smallClock.restart();
                elVals = labelElements<dim,T>(mp, damage, mb,0.1,1.0);
                gsInfo<<"Labelling level "<<i<<" took "<<smallClock.stop()<<" seconds\n";
                if (gsAsVector<T>(elVals).sum() > 0)
                {
                    smallClock.restart();
                    tmpArea = refineMesh<dim,T>(mb,elVals,mesherOptions);
                    gsInfo<<"Refining mesh took "<<smallClock.stop()<<" seconds\n";
                }

                tmpArea /= (mb.basis(0).support().col(1)-mb.basis(0).support().col(0)).prod();
                markedArea = math::max(markedArea,tmpArea);
                basis_size = mb.basis(0).size();
                refined = basis_size > basis_size_old;
                if (!refined)
                    break;
            }
            gsInfo<<"Marked area: "<<markedArea<<"\n";
            refined &= markedArea > mesherOptions.askReal("SizeRatio",1.01);

            // basis_size_ratio = (T)basis_size/basis_size_old;
            // gsInfo<<"Old mesh size: "<<basis_size_old<<", new mesh size: "<<basis_size<<", ratio = "<<basis_size_ratio<<"\n";
            // refined &= basis_size_ratio > mesherOptions.askReal("SizeRatio",1.05);

            // =========================================================================
            // PROJECT SOLUTIONS
            smallClock.restart();
            gsMatrix<T> projCoefs;
            // Geometry
            gsQuasiInterpolate<T>::localIntpl(mb.basis(0),mp.patch(0),projCoefs);
            mp.clear();
            mp.addPatch(mb.basis(0).makeGeometry(give(projCoefs)));
            mp_def = mp;
            // Displacement
            gsQuasiInterpolate<T>::localIntpl(mb.basis(0),displacement.patch(0),projCoefs);
            displacement.clear();
            displacement.addPatch(mb.basis(0).makeGeometry(give(projCoefs)));
            // gsQuasiInterpolate<T>::localIntpl(mb.basis(0),displacement_old.patch(0),projCoefs);
            // displacement_old.clear();
            // displacement_old.addPatch(mb.basis(0).makeGeometry(give(projCoefs)));
            // Damage
            gsQuasiInterpolate<T>::localIntpl(mb.basis(0),damage.patch(0),projCoefs);
            damage.clear();
            damage.addPatch(mb.basis(0).makeGeometry(give(projCoefs)));
            // gsQuasiInterpolate<T>::localIntpl(mb.basis(0),damage_old.patch(0),projCoefs);
            // damage_old.clear();
            // damage_old.addPatch(mb.basis(0).makeGeometry(give(projCoefs)));
            T ptime = smallClock.stop();
            stepTimes.projectionTime += ptime;
            gsInfo<<"Projection took "<<ptime<<" seconds\n";

        // for (refIt = 0; refIt!=mesherOptions.askInt("MaxRefIterations",5) && refined; refIt++)

            if (!refined)
                break;
            else if (refIt==mesherOptions.askInt("MaxRefIterations",5)-1)
            {
                gsWarn<<"Maximum number of refinement itertions reached\n";
                break;
            }
            else
                refIt++;
        }
        numIt_ref = refIt+1;

        // =========================================================================
        // PLOT
        if (plot && step%plotmod==0)
        {
            std::string filename;

            filename = "mesh_" + util::to_string(step);
            gsMesh<T> mesh(mb.basis(0));
            writeSingleCompMesh(mb.basis(0),mp.patch(0),outputdir+filename,1);
            meshCollection.addPart(filename,step,"Mesh",0);

            filename = "damage_" + util::to_string(step);
            gsWriteParaview(mp,damage,outputdir+filename,100000);
            // gsField<T> damage_step(zone,damage,false);
            // gsWriteParaview(damage_step,filename,1000);
            filename += "0";
            damageCollection.addPart(filename,step,"Solution",0);

            gsMaterialEval<T,gsMaterialOutput::Psi> Psi(&material,mp,mp_def);
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
        gsSolidAssembler<dim,T,gsLinearDegradedMaterial<T>> fullElAssembler(mp,mb,bc_u_dummy,&material);
        fullElAssembler.options().setReal("ExprAssembler.quA",1.0);
        fullElAssembler.options().setInt ("ExprAssembler.quB",0);
        fullElAssembler.initialize();

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

        file.open(outputdir+"results.txt",std::ios::app);
        file<<step<<","<<ucurr<<","<<-Fx<<","<<-Fy<<","
            <<(0.5 * ufull.transpose() * fullElAssembler.matrix() * ufull).value()<<","
            <<(0.5 * D.transpose() * QPhi * D).value() + (D.transpose() * q).value()<<","
            <<stepTimes.elAssemblyTime<<","<<stepTimes.elSolverTime<<","
            <<stepTimes.pfAssemblyTime<<","<<stepTimes.pfSolverTime<<","
            <<stepTimes.projectionTime<<","
            <<basis_size<<","<<markedArea<<","
            <<totIt_el<<","<<totIt_pf<<","
            <<numIt_stag<<","<<numIt_ref<<"\n";
        file.close();

        // =========================================================================
        // INCREMENT STEP

        // displacement_old = displacement;
        // damage_old = damage;
        if (ucurr == uend || math::abs(ucurr-uend) < 1e-10)
            break;

        ucurr += (ucurr+ustep > utrans) ? ustep/ured : ustep;
        ucurr = math::min(ucurr,uend);
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