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
    gsInfo << "Output directory "<< outputdir <<"\n";
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
    gsInfo<<"damage size before refinement: "<<damage.patch(0).coefs().size()<<"\n";
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
    // Min time [s]
    controlParameters.addReal("tend", "Maximum time", 15e-7);
    // Max time [s]
    controlParameters.addReal("tmin", "Initial time", 0.0);
    // Time step [s]
    controlParameters.addReal("tstep", "Time step", 1e-7);
    // Step transition [s]
    controlParameters.addReal("ttrans", "Step transition", controlParameters.getReal("tend"));
    // Step reduction factor [-]
    controlParameters.addReal("tred", "Step reduction factor", 1.0);
    // Maximum number of iterations
    controlParameters.addInt("maxIt", "Maximum number of iterations", 10000);
    // Maximum number of iterations for elasticity problem
    controlParameters.addInt("maxItEl", "Maximum number of iterations for elasticity problem", 1);
    // Maximum number of iterations for phase-field problem
    controlParameters.addInt("maxItPf", "Maximum number of iterations for phase-field problem", 1);
    // Tolerance for the elasticity problem
    controlParameters.addReal("tolEl", "Tolerance for the elasticity problem", 1e-5);
    // Tolerance for the phase-field problem
    controlParameters.addReal("tolPf", "Tolerance for the phase-field problem", 1e-5);
    // Staggered tolerance
    controlParameters.addReal("tol", "Tolerance for the staggered scheme", 1e-5);

    //// Solver parameters
    gsOptionList solverParameters;
    if (fd_pars.hasLabel("solver"))
        fd_pars.getLabel("solver", solverParameters);

    // //// Boundary conditions
    // gsBoundaryConditions<> bc_u;
    // fd_pars.getLabel("BCs_u", bc_u);

    // //// Boundary conditions
    // gsBoundaryConditions<> bc_d;
    // fd_pars.getLabel("BCs_d", bc_d);

    short_t dim = mp_ini.domainDim();
    
    
    gsBoundaryConditions<> bc_u;

    // Boundary conditions
    real_t sigma = 1.0;
    std::vector<std::string> bcFunctionTop(dim);
    std::vector<std::string> bcFunctionBottom(dim);
    bcFunctionBottom[1] = "-u"; // the y component [1], the x component [0]
    bcFunctionTop[1] = "u";
    // for (short_t d=1; d<dim; ++d)
    bcFunctionTop[0] = bcFunctionBottom[0] = "0"; // the x component is zero
    gsFunctionExpr<> sigma_top(bcFunctionTop,dim);
    gsFunctionExpr<> sigma_bottom(bcFunctionBottom,dim);
    sigma_top.set_u(sigma);
    sigma_bottom.set_u(sigma);
    bc_u.addCondition(boundary::north,condition_type::neumann,&sigma_top );
    bc_u.addCondition(boundary::south,condition_type::neumann,&sigma_bottom);
    // bc_u.addCondition(boundary::west,condition_type::dirichlet,0,0,false,0); //vertical constraint
    // bc_u.addCondition(boundary::east,condition_type::dirichlet,0,0,false,0); //vertical constraint
    // bc_u.addCondition(boundary::west, condition_type::dirichlet, 0, 0, false, 0); // fix x
    // bc_u.addCondition(boundary::east, condition_type::dirichlet, 0, 0, false, 0); // fix x
    if (dim==3)
    {
        bc_u.addCondition(boundary::back,condition_type::dirichlet,0,0,false,2); //vertical constraint
        bc_u.addCondition(boundary::front,condition_type::dirichlet,0,0,false,2); //vertical constraint
    }

    gsBoundaryConditions<> bc_d;


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
// #pragma omp parallel for
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
    // Density [tonn/mm^3]
    T rho = materialParameters.getReal("rho");
    gsInfo<<"Material parameters: E="<<E<<", nu="<<nu<<", Gc="<<Gc<<", l0="<<l0<<", rho="<<rho<<"\n";
    // Order of the phase-field model
    index_t order = materialParameters.getInt("order");
    // AT model
    index_t AT = materialParameters.getInt("AT");

    GISMO_ASSERT(order == 2 || order == 4, "Please specify the order of the model (2 or 4).");
    GISMO_ASSERT(AT == 1 || AT == 2, "Please specify the AT model (1 or 2).");

    ////////////////////////////////////////////////////////////////////////////////////
    // Boundary control parameters
    ////////////////////////////////////////////////////////////////////////////////////
    // Min time [s]
    T tmin = controlParameters.getReal("tmin");
    // Max time [s]
    T tend = controlParameters.getReal("tend");
    // Time step [s]
    T tstep = controlParameters.getReal("tstep");
    // Step transition [s]
    T ttrans = controlParameters.askReal("ttrans",tend);
    // Step reduction factor [-]
    T tred = controlParameters.askReal("tred",1.);
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

    T tcurr_old = tmin;
    T tcurr = tmin + tstep;

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

    // // Boundary conditions
    // T sigma = 1.0;
    // std::vector<std::string> bcFunctionBottom(dim);
    // std::vector<std::string> bcFunctionTop(dim);
    // bcFunctionBottom[0] = "-u";
    // bcFunctionTop[0] = "u";
    // for (short_t d=1; d<dim; ++d)
    //     bcFunctionBottom[d] = bcFunctionTop[d] = "0";
    // gsFunctionExpr<T> sigma_bottom(bcFunctionBottom,dim);
    // gsFunctionExpr<T> sigma_top(bcFunctionTop,dim);
    // sigma_bottom.set_u(sigma);
    // sigma_top.set_u(sigma);
    // bc_u.addCondition(boundary::south,condition_type::neumann,&sigma_bottom );
    // bc_u.addCondition(boundary::north,condition_type::neumann,&sigma_top);
    // bc_u.addCondition(boundary::west,condition_type::dirichlet,0,0,false,1); //vertical constraint
    // bc_u.addCondition(boundary::east,condition_type::dirichlet,0,0,false,1); //vertical constraint
    // if (dim==3)
    // {
    //     bc_u.addCondition(boundary::back,condition_type::dirichlet,0,0,false,2); //vertical constraint
    //     bc_u.addCondition(boundary::front,condition_type::dirichlet,0,0,false,2); //vertical constraint
    // }

    bc_u.setGeoMap(mp);
    bc_d.setGeoMap(mp);
    ///////////////////////////////////////////////////////////////////////////////////////
    //INITIALIZATION///////////////////////////////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////

    // Initialize the solution (deformed geometry)
    gsMultiPatch<T> mp_def = mp;
    gsMultiPatch<T> displacement = mp;
    gsMultiPatch<T> velocity = mp;
    gsMultiPatch<T> acceleration = mp;
    gsInfo<<"Initializing displacement, velocity and acceleration to zero.\n";
    for (index_t p=0; p<mp.nPatches(); ++p)
    {
        displacement.patch(p).coefs().setZero();
        velocity.patch(p).coefs().setZero();
        acceleration.patch(p).coefs().setZero();
    }
    gsInfo<<"Worked!\n";


    // Initialize the material
    gsLinearDegradedMaterial<T> material(E,nu,rho,damage,dim);

    gsBoundaryConditions<T> bc_u_dummy;
    bc_u_dummy.setGeoMap(mp);

    // Initialize the phase-field assembler
    gsPhaseFieldAssemblerBase<T> * pfAssembler;

    //////////////////////////////////////////////////////////////////////////
    // INITIALIZE THE MESH ///////////////////////////////////////////////////
    //////////////////////////////////////////////////////////////////////////

    std::vector<T> elVals;

    //////////////////////////////////////////////////////////////////////////
    // SOLVE THE PROBLEM/////////////////////////////////////////////////////
    //////////////////////////////////////////////////////////////////////////
    
#ifdef GISMO_WITH_PARDISO
    typename gsSparseSolver<T>::PardisoLDLT solver;
    gsInfo<<"Using Pardiso direct solver\n";
#else
    typename gsSparseSolver<T>::CGDiagonal solver;
    gsInfo<<"Using CG diagonal preconditioned iterative solver\n";
#endif

    times<T> stagTimes;
    times<T> stepTimes;
    stagTimes.reset();
    stepTimes.reset();

    gsSparseMatrix<T> elMatrix;
    gsMatrix<T> elRhs;

    gsMatrix<T> R, Fext;
    gsMatrix<T> D, delta_D;
    gsSparseMatrix<T> Q, QPhi, QPsi;
    gsMatrix<T> q, qpsi;

    T elAssemblyTime = 0.0;
    T elSolverTime = 0.0;
    T pfAssemblyTime = 0.0;
    T pfSolverTime = 0.0;
    T iterationTime  = 0.0;

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

    // gsMultiPatch<T> displacement_old, damage_old;
    // displacement_old = displacement;
    // damage_old = damage;
    index_t numIt_el = 0, totIt_el = 0; //numIt: per staggered iteration, totIt: total number of iterations across all staggered iterations
    index_t numIt_pf = 0, totIt_pf = 0; //numIt: per staggered iteration, totIt: total number of iterations across all staggered iterations
    index_t numIt_ref = 0; // number of refinement iterations per load step
    index_t numIt_stag = 0; // number of staggered iterations per load step

    real_t gamma   = 0.5;
    real_t beta    = 0.25;
    real_t Rnorm, R0;
    real_t Unorm, U0;
    T Fnorm;
    Rnorm = Fnorm = 1;
    T dt;

    gsMatrix<T> u_new, u_old, udot_new, udot_old, uddot_new, uddot_old, delta_u;

    

    index_t step = 0;

        gsInfo<<"Plotting results for step "<<step<<"\n";
        std::string filename2;

        filename2 = "meshsfsdfs_" + util::to_string(step);
        gsMesh<T> mesh(mb.basis(0));
        writeSingleCompMesh(mb.basis(0),mp.patch(0),outputdir+filename2,1);
        meshCollection.addPart(filename2,step,"Mesh",0);
        gsInfo<<"Plotting results for step "<<step<<"\n";

    while (tcurr<=tend)
    {
        dt = tcurr - tcurr_old;
        numIt_ref = numIt_stag = 0;
        totIt_el = totIt_pf = 0;
        stepTimes.reset();
        stagTimes.reset();

        bool refined = true;
        index_t basis_size_old, basis_size;
        T basis_size_ratio;
        T markedArea = 0., tmpArea = 0.;
        gsInfo<<"===========================================================================================================================\n";
        gsInfo<<"Load step "<<step<<": t = "<<tcurr<<"\n";
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
            // Initialize solutions in the new mesh
            elAssembler.constructSolution(displacement,u_old);
            elAssembler.constructSolution(velocity,udot_old);
            elAssembler.constructSolution(acceleration,uddot_old);

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
            // displ.set_u(ucurr);
            // elAssembler.initialize();
            gsInfo<<"Load step "<<step<<": t = "<<tcurr<<"\n\n";

            delta_u.setZero();
            delta_D.setZero();

            // if (step == 0)
            // {
            //     gsInfo<<"Initial step: computing initial acceleration.\n";
            //     elAssembler.initialize();
            //     elAssembler.assembleMass();
            //     elAssembler.assemble(u_old);
            //     gsSparseMatrix<> M;
            //     elAssembler.matrix_into(M);
            //     solver.compute(M);
            //     elAssembler.rhs_into(Fext);
            //     uddot_old = solver.solve(Fext);
            //     gsInfo<<"Initialization done!\n";
            //     gsInfo<< "Fext max: "<< Fext.maxCoeff()<<" Fext min: "<<Fext.minCoeff()<<"\n";
            // }

            gsInfo<<"uddot max: "<< uddot_old.maxCoeff()<<" uddot min: "<<uddot_old.minCoeff()<<"\n";
            
            // Prediction step (IGA book Eqs. (6.44)-(6.46))
            udot_new = udot_old;
            uddot_new = (gamma-1)/gamma * uddot_old;
            u_new = u_old + dt * udot_old + 0.5*math::pow(dt,2) * ((1-2*beta) * uddot_old + 2*beta * uddot_new);

            gsInfo<<"udot_bnew max: "<< udot_new.maxCoeff()<<" udot_bmin: "<<udot_new.minCoeff()<<"\n";
            gsInfo<<"uddot_bnew max: "<< uddot_new.maxCoeff()<<" uddot_bmin: "<<uddot_new.minCoeff()<<"\n";
            gsInfo<<"u_bnew max: "<< u_new.maxCoeff()<<" u_bmin: "<<u_new.minCoeff()<<"\n";
            
            Unorm = u_new.norm();
            U0 = (Unorm > 0) ? Unorm : 1.0; // Avoid division by zero
            Rnorm = R0 = 1;

            delta_D.setZero(D.rows(),1);

            index_t stagIt = 0;

            while(true) // staggered loop 
            {
                stagTimes.reset();
                if (stagIt==0) stagTimes.pfAssemblyTime = pfAssemblyTime0; // save the time of the first assembly

                bigClock.restart();
                gsInfo<<"    --------------------------Staggered Iteration: "<<PRINT(4)<<stagIt<<"--------------------------\n";
                gsInfo<<"    ---------------------------------ELASTICITY----------------------------------\n";
                // gsInfo << "    | " << PRINT(6) << "It." << PRINT(14) << "||R||" << PRINT(14) << "||F||" << PRINT(14) << "||R||/||F||" << PRINT(14) << "||U||" << PRINT(20) << "cum. assembly [s]" << PRINT(20) << "cum. solver [s]" << PRINT(20) << "it. solver [s]" << PRINT(20) << "solver it." << "|\n";
                gsInfo <<  "    | " << PRINT(18) << "||R||" << PRINT(18) << "||R||/||R0||" << PRINT(18) << "||ΔU||/||U0||" << PRINT(18) << "||dA||/||A||" << PRINT(18) << "||U||" << PRINT(18) << "||V||" << PRINT(18) << "||A||" << PRINT(20) << "cum. assembly [s]" << PRINT(20) << "cum. solver [s]" << "\n";

                material.setParameter(2,damage);
                elAssembler.initialize();
                elAssembler.assembleMass();
                gsSparseMatrix<> M = elAssembler.matrix();
                // M.setZero();
                gsSparseMatrix<T> K;
                elAssembler.initialize();
                // Assemble the elasticity problem
                smallClock.restart();
                elAssembler.assemble(u_new);
                elAssemblyTime += smallClock.stop();
                elAssembler.matrix_into(K);
                elAssembler.rhs_into(Fext);


                

                gsInfo << "Fext max: "<< Fext.maxCoeff()<<" Fext min: "<<Fext.minCoeff()<<"\n";

                
                R = K * u_new - Fext;
                Rnorm = R.norm();
                if (stagIt == 0)
                    R0 = (Rnorm > 0) ? Rnorm : 1.0; // For exit criterion (eq. (18))
            
                smallClock.restart();
                solver.compute(K);
                delta_u = solver.solve(-R);
                elSolverTime += smallClock.stop();

                u_new += delta_u;

                Unorm = u_new.norm();
                T DeltaUnorm = delta_u.norm();
                
                Fnorm = Fext.norm();
                Fnorm = (Fnorm == 0) ? 1 : Fnorm;
                // K+= 1/(beta*math::pow(dt,2)) * M; // Eq. (21) - Greco et al. 2025
                
                // Print values of the Elasticity solve 
                gsInfo << "    | " << PRINT(18) << Rnorm << PRINT(18) << Rnorm/R0 << PRINT(18) << DeltaUnorm/U0 << PRINT(18) << (uddot_new-uddot_old).norm()/uddot_new.norm() << PRINT(18) << Unorm << PRINT(18) << udot_new.norm() << PRINT(18) << uddot_new.norm() << PRINT(20) << elAssemblyTime << PRINT(20) << elSolverTime << "\n";
                
                // Recompute the residual for the staggered check
                smallClock.restart();
                elAssembler.assemble(u_new);
                elAssemblyTime += smallClock.stop();
                elAssembler.matrix_into(K);
                elAssembler.rhs_into(Fext);
                R = K * u_new - Fext;
                Rnorm = R.norm();

                
                gsInfo<<"    * Staggered check:\n";
                gsInfo <<"    | " << PRINT(18) << "||R||" 
                    << PRINT(18) << "||R||/||R0||" 
                    << PRINT(18) << "||dU||" 
                    << PRINT(18) << "||dU||/||U0||" << "\n";

                gsInfo <<"    | " << PRINT(18)  << Rnorm 
                                    << PRINT(18) << Rnorm/R0 
                                    << PRINT(18) << DeltaUnorm 
                                    << PRINT(18) << DeltaUnorm/U0 << "\n";

                elAssembler.constructSolution(u_new,displacement);
                // elAssembler.constructSolution(udot_new,velocity);
                // elAssembler.constructSolution(uddot_new,acceleration);

                for (size_t p=0; p!=mp.nPatches(); ++p)
                    mp_def.patch(p).coefs() = mp.patch(p).coefs() + displacement.patch(p).coefs();

                // Initialize the function for the elastic energy
                gsMaterialEval<T,gsMaterialOutput::Psi> Psi(&material,mp,mp_def);

                // index_t maxIt = 2;
                if (Rnorm/R0 < 1e-5 && DeltaUnorm/U0 < 1e-4)
                    break;
                else if (stagIt == maxIt-1)
                    break;
                    // GISMO_ERROR("Staggered iterations problem did not converge.");
                else
                    stagIt++;

                // ==================================================================================
                // Update damage spline
                pfAssembler->constructSolution(D,damage);


            }
            numIt_stag += stagIt+1;

            gsInfo<<"\n";
            gsInfo<<"Converged with ||R||/||F|| = "<<Rnorm/Fnorm<<" < "<<tol<<" ||D|| = "<<D.norm()<<" ||U|| = "<<u_new.norm()<<"\n";
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
            // ========================================================================
            // New variables for dynamics
            // gsQuasiInterpolate<T>::localIntpl(mb.basis(0),displacement_old.patch(0),projCoefs);
            // displacement_old.clear();
            // displacement_old.addPatch(mb.basis(0).makeGeometry(give(projCoefs)));
            gsQuasiInterpolate<T>::localIntpl(mb.basis(0),velocity.patch(0),projCoefs);
            velocity.clear();
            velocity.addPatch(mb.basis(0).makeGeometry(give(projCoefs)));
            gsQuasiInterpolate<T>::localIntpl(mb.basis(0),acceleration.patch(0),projCoefs);
            acceleration.clear();
            acceleration.addPatch(mb.basis(0).makeGeometry(give(projCoefs)));
            // ========================================================================
            // Damage
            gsInfo<<"damage max before projection: "<<damage.patch(0).coefs().maxCoeff()<<"\n";
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
            gsInfo<<"refined: "<<refined<<"\n";
            if (!refined)
                break;
            else if (refIt==mesherOptions.askInt("MaxRefIterations",5)-1)
            {
                gsWarn<<"Maximum number of refinement iterations reached\n";
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

            gsInfo<<"Plotting results for step "<<step<<"\n";
            filename = "mesh_" + util::to_string(step);
            gsMesh<T> mesh(mb.basis(0));
            writeSingleCompMesh(mb.basis(0),mp.patch(0),outputdir+filename,1);
            meshCollection.addPart(filename,step,"Mesh",0);
            gsInfo<<"Plotting results for step "<<step<<"\n";

            filename = "displacement_"+util::to_string(step);
            gsInfo<<"disp max after projection: "<<displacement.patch(0).coefs().maxCoeff()<<"\n";
            gsInfo<<"disp min after projection: "<<displacement.patch(0).coefs().minCoeff()<<"\n";
            gsWriteParaview(mp,displacement,outputdir+filename,100000);
            filename += "0";
            displCollection.addPart(filename,step,"Solution",0);
        }

        u_old = u_new;
        udot_old = udot_new;
        uddot_old = uddot_new;

        tcurr_old = tcurr;
        tcurr += (tcurr+tstep > ttrans) ? tstep/tred : tstep;
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