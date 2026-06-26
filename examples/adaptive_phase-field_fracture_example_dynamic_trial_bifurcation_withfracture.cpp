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
           gsOptionList & mesherOptions,
           gsMultiPatch<T> & mp,
           gsMultiPatch<T> & damage,
           gsBoundaryConditions<T> & bc_u,
           gsBoundaryConditions<T> & bc_d,
           bool plot,
           index_t plotmod,
           bool plotMesh,
           std::string & outputdir);

int main(int argc, char *argv[])
{
    //! [Parse command line]
    bool plot = false;
    bool plotMesh = false;
    index_t numHRef = 0;
    index_t numElev = 0;
    index_t order = -1;
    index_t AT = -1;
    index_t plotmod = 1;
    index_t dimension = 2;
    std::string output;
    
    std::string parInput;
    std::string geoInput;
    std::string damageInput;
    std::string inputDir;

    gsCmdLine cmd("Tutorial on solving a Linear Elasticity problem.");
    cmd.addInt("e", "numElev","Number of degree elevation steps to perform before solving",numElev);
    cmd.addInt("r", "numHRef","Number of Uniform h-refinement loops", numHRef);
    cmd.addInt("O", "order","Order of the basis functions", order);
    cmd.addInt("A", "AT","AT-1 or AT-2 model", AT);
    cmd.addInt("p", "plotmod","Modulo for plotting", plotmod);
    // cmd.addInt("d", "dimension","Dimension of the problem", dimension);
    cmd.addSwitch("plot","Create a ParaView visualization file with the solution", plot);
    cmd.addSwitch("plotMesh","Create a ParaView visualization file with the adaptive mesh", plotMesh);
    cmd.addString("o", "output", "Output directory", output);
    cmd.addString("i", "parInput", "Input XML file", parInput);
    cmd.addString("g", "geometry", "Geometry file", geoInput);
    cmd.addString("d", "damage", "Damage file", damageInput);
    cmd.addString("I", "inputDir", "Input directory", inputDir);

    try { cmd.getValues(argc,argv); } catch (int rv) { return rv; }

    GISMO_ASSERT(order == 2 || order == 4, "Please specify the order of the model (2 or 4).");
    GISMO_ASSERT(AT == 1 || AT == 2, "Please specify the AT model (1 or 2).");
    // GISMO_ASSERT(dimension == 2 || dimension == 3, "Please specify the dimension of the problem (2 or 3).");

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

    // if (numElev > 0)
    //     mp_ini.degreeIncrease(numElev);
    
    //     gsInfo<< "c3\n";

    // for (index_t i = 0; i<numHRef; ++i)
    //     mp_ini.uniformRefine(1);

    if (plot) gsWriteParaview(mp_ini,outputdir+"mp",10,true);

    gsFileData<> fd_damage(damageInput.empty() ? inputDir + "damage.xml" : damageInput);
    gsMultiPatch<> damage;
    fd_damage.getFirst(damage);
    gsInfo<<"damage size before refinement: "<<damage.patch(0).coefs().size()<<"\n";
    if (plot) gsWriteParaview(mp_ini,damage,outputdir+"initial_damage",100000);


    gsFileData<> fd_pars(geoInput.empty() ? inputDir + "parameters.xml" : parInput);
    gsInfo << "Parameter file "<< parInput <<" read.\n";
    //// Material parameters
    gsOptionList materialParameters;
    fd_pars.getLabel("material", materialParameters);
    gsInfo<<"Material parameters:\n"<<materialParameters<<"\n";

    //// Boundary control parameters
    gsOptionList controlParameters;
    fd_pars.getLabel("control", controlParameters);
    gsInfo<<"Control parameters:\n"<<controlParameters<<"\n";

    gsOptionList mesherOptions;
    fd_pars.getLabel("meshing", mesherOptions);
    gsInfo<<"Mesher options:\n"<<mesherOptions<<"\n";

    short_t dim = mp_ini.domainDim();
    gsBoundaryConditions<> bc_u, bc_d;

    // Boundary conditions
    fd_pars.getLabel("BCs_u", bc_u);
    // fd_pars.getLabel("BCs_d", bc_d);

    gsInfo<< "Displacement BC: \n"<< bc_u;
    gsInfo<< "Damage BC: \n"<< bc_d; 

    ///////////////////////////////////////////////////////////////////////////////////////
    // Call the dimensional solver
    ///////////////////////////////////////////////////////////////////////////////////////
    switch (mp_ini.domainDim())
    {
        case 2:
            solve<2>(materialParameters,controlParameters,mesherOptions,mp_ini,damage,bc_u,bc_d,plot,plotmod,plotMesh,outputdir);
            break;
        case 3:
            solve<3>(materialParameters,controlParameters,mesherOptions,mp_ini,damage,bc_u,bc_d,plot,plotmod,plotMesh,outputdir);
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
           gsOptionList & mesherOptions,
           gsMultiPatch<T> & mp_ini,
           gsMultiPatch<T> & damage,
           gsBoundaryConditions<T> & bc_u,
           gsBoundaryConditions<T> & bc_d,
           bool plot,
           index_t plotmod,
           bool plotMesh,
           std::string & outputdir)
{

    
    ////////////////////////////////////////////////////////////////////////////////////
    // Load parameters (paper Borden et al. 2012)
    ////////////////////////////////////////////////////////////////////////////////////
    // Young's modulus [N/mm^2]
    T E = materialParameters.getReal("E");
    // Poisson's ratio [-]
    T nu = materialParameters.getReal("nu"); 
    // Toughness [N/mm]
    T Gc = materialParameters.getReal("Gc");
    // Internal length [mm]
    T l0 = materialParameters.getReal("l0");
    // Density [kg/mm^3]
    T rho = materialParameters.getReal("rho");
    // Order of the phase-field model
    index_t order = materialParameters.getInt("order");
    // AT model
    index_t AT = materialParameters.getInt("AT");

    gsInfo<< "Material parameters loaded:\n";
    gsInfo<< "E: " << E << "\n";
    gsInfo<< "nu: " << nu << "\n";
    gsInfo<< "Gc: " << Gc << "\n";
    gsInfo<< "l0: " << l0 << "\n";
    gsInfo<< "rho: " << rho << "\n";
    gsInfo<< "Phase-field order: " << order << "\n";
    gsInfo<< "AT model: " << AT << "\n";

    GISMO_ASSERT(order == 2 || order == 4, "Please specify the order of the model (2 or 4).");
    GISMO_ASSERT(AT == 1 || AT == 2, "Please specify the AT model (1 or 2).");

    ////////////////////////////////////////////////////////////////////////////////////
    // Boundary control parameters
    ////////////////////////////////////////////////////////////////////////////////////
    // Min time [s]
    T tmin = controlParameters.getReal("umin");
    // Max time [s]
    T tend = controlParameters.getReal("uend");
    // Time step [s]
    T tstep = controlParameters.getReal("ustep");
    // Step transition [s]
    T ttrans = controlParameters.askReal("utrans",tend);
    // Step reduction factor [-]
    T tred = controlParameters.askReal("ured",1.);
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


    gsInfo <<  "tolEl: " << tolEl << " tolPf: " << tolPf << " tol: " << tol <<"\n";
 
    // // Fixed side patch id
    // index_t fixedSidePatch = controlParameters.getInt("patchId");
    // // Fixed side id
    // index_t fixedSideId = controlParameters.getInt("side");
    // // Fixed side direction
    // index_t fixedSideDir = controlParameters.getInt("direction");
    // Function on the boundary
    // std::string bcFunctionLeft = controlParameters.askString("function","-u");
    // std::string bcFunctionRight = controlParameters.askString("function","u");

    T tcurr_old = tmin;
    T tcurr = tmin + tstep;

    ///////////////////////////////////////////////////////////////////////////////////////
    //PROBLEM SETUP////////////////////////////////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////

    // Convert to THB
    int maxRefLvl;
    gsMultiPatch<T> mp;
    bool adaptive_switch;
    for (index_t i = 0; i < mp_ini.nPatches(); ++i)
    {
        // Check if tensor basis
        if      ((dynamic_cast<const gsTensorBSpline<dim,T> *>(&mp_ini.patch(i))))
        {
            gsInfo << "TB basis: non-adaptive solution\n";
            adaptive_switch = false;
            // Create a THB spline basis
            const gsTensorBSpline<dim,T> & tb = static_cast<const gsTensorBSpline<dim,T> &>(mp_ini.patch(i));
            gsTHBSpline<dim,T> thb(tb);
            mp.addPatch(memory::make_unique(thb.clone().release()));
        }
        else if ((dynamic_cast<const gsTHBSpline<dim,T> *>(&mp_ini.patch(i))))
        {
            adaptive_switch = true;
            gsInfo << "THB-spline basis: adaptivity triggered\n";
            const gsTHBSpline<dim,T> & thb = static_cast<const gsTHBSpline<dim,T> &>(mp_ini.patch(i));
            mp.addPatch(memory::make_unique(thb.clone().release()));
            maxRefLvl = thb.basis().maxLevel();
        }
        else
            GISMO_ERROR("The basis is not a TB-spline basis or THB-spline basis.");
    }

    gsInfo<<mesherOptions<<"\n";


    // Construct the basis
    gsMultiBasis<T> mb(mp);
    gsInfo<<"The basis has size "<<mb.size()<<" and degree "<<mb.degree()<<"\n";
    for (size_t b=0; b!=mb.nBases(); b++)
        gsInfo<<"Basis "<<b<<":\n"<<mb.basis(b)<<"\n";
    
    // // Boundary conditions
    // gsFunctionExpr<T> displ(bcFunction,dim);
    // displ.set_u(ucurr);
    // bc_u.addCondition(fixedSidePatch,fixedSideId,condition_type::dirichlet,&displ,0,false,fixedSideDir);
    bc_u.setGeoMap(mp);
    bc_d.setGeoMap(mp);

    ///////////////////////////////////////////////////////////////////////////////////////
    //INITIALIZATION///////////////////////////////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////

    // Initialize the solution (deformed geometry)
    gsMultiPatch<T> mp_def            = mp;
    gsMultiPatch<T> displacement      = mp;
    gsMultiPatch<T> velocity          = mp;
    gsMultiPatch<T> acceleration      = mp;
    gsMultiPatch<T> displacement_old  = mp;
    gsMultiPatch<T> velocity_old      = mp;
    gsMultiPatch<T> acceleration_old  = mp;
    gsMultiPatch<T> damage_old        = damage;

    for (index_t p = 0; p < mp.nPatches(); ++p)
    {
        displacement.patch(p).coefs().setZero();
        velocity.patch(p).coefs().setZero();
        acceleration.patch(p).coefs().setZero();
        displacement_old.patch(p).coefs().setZero();
        velocity_old.patch(p).coefs().setZero();
        acceleration_old.patch(p).coefs().setZero();
    }

    // Initialize the material
    gsLinearDegradedMaterial<T> material(E,nu,rho,damage,dim);
    
    // Initialize the elasticity assembler
    gsSolidAssembler<dim,T,gsLinearDegradedMaterial<T>> elAssembler(mp,mb,bc_u,&material);
    elAssembler.options().setReal("ExprAssembler.quA",1.0);
    elAssembler.options().setInt ("ExprAssembler.quB",1);
    elAssembler.options().setInt("ExprAssembler.DirichletValues",dirichlet::l2Projection);
    elAssembler.initialize();
    elAssembler.assemble();

    gsBoundaryConditions<T> bc_u_dummy;
    bc_u_dummy.setGeoMap(mp);
    
    // Initialize the mass assembler
    elAssembler.assembleMass();
    gsSparseMatrix<> M = elAssembler.matrix();

    // Initialize the phase-field assembler
    gsPhaseFieldAssemblerBase<T> * pfAssembler;
    if      (order == 2 && AT == 1)
        pfAssembler = new gsPhaseFieldAssembler<T,PForder::Second,PFmode::AT1>(mp,mb,bc_d);
    else if (order == 4 && AT == 1)
    {
        pfAssembler = new gsPhaseFieldAssembler<T,PForder::Fourth,PFmode::AT1>(mp,mb,bc_d);
        pfAssembler->options().setReal("cw",3.1615);
        pfAssembler->options().setReal("chi",0.0625);
    }
    else if (order == 2 && AT == 2)
        pfAssembler = new gsPhaseFieldAssembler<T,PForder::Second,PFmode::AT2>(mp,mb,bc_d);
    else if (order == 4 && AT == 2)
    {
        pfAssembler = new gsPhaseFieldAssembler<T,PForder::Fourth,PFmode::AT2>(mp,mb,bc_d);
    }   
    else
        GISMO_ERROR("Invalid order and/or AT model");

    pfAssembler->options().setReal("l0",l0);
    pfAssembler->options().setReal("Gc",Gc);
    pfAssembler->initialize();

    gsMatrix<> D_new, D_old, delta_D;
    gsMatrix<> u_new,u_old,udot_new,udot_old,uddot_new,uddot_old,delta_u;
    u_new.setZero();
    u_old.setZero();
    udot_new.setZero();
    udot_old.setZero();
    uddot_new.setZero();
    uddot_old.setZero();

    std::vector<T> elVals;          

#ifdef GISMO_WITH_PARDISO
    typename gsSparseSolver<T>::PardisoLDLT solver;
#else
    typename gsSparseSolver<T>::CGDiagonal solver;
#endif

    times<T> stagTimes;
    times<T> stepTimes;
    stagTimes.reset();
    stepTimes.reset();

    gsSparseMatrix<T> K;
    gsMatrix<T> R, Fext;

    T elAssemblyTime = 0.0;
    T elSolverTime = 0.0;
    T pfAssemblyTime = 0.0;
    T pfSolverTime = 0.0;
    T iterationTime  = 0.0;
    T projTime = 0.0;
    T labelTime = 0.0;
    T refTime = 0.0;
    T totalTime = 0.0;

    T energy_D = 0.0, energy_E = 0.0;

    gsSparseMatrix<T> Q, QPhi, QPsi;
    gsMatrix<T> q, qpsi;
    // Phase-field assembly can already be performed since some terms are independent of the solutions
    pfAssembler->assemblePhi();
    pfAssembler->matrix_into(QPhi);
    pfAssembler->rhs_into(q);
    // Initialize the damage solution vector
    pfAssembler->constructSolution(damage,D_old);

    index_t step = 0;

    gsParaviewCollection damageCollection(outputdir+"damage");
    gsParaviewCollection psiCollection(outputdir+"Psi");
    gsParaviewCollection displCollection(outputdir+"displacement");
    gsParaviewCollection PwaveCollection(outputdir+"Pwave");
    gsParaviewCollection meshCollection(outputdir+"mesh");
    gsStopwatch smallClock, bigClock;

    std::ofstream file(outputdir+"results.txt");
    file<<"step,Fx,Fy,E_u,E_d\n";
    file.close();

    /* @todo: add option list from file for dynamic parameters */
    // gsOptionList dynamicParameters;

    index_t numIt_el = 0, totIt_el = 0; //numIt: per staggered iteration, totIt: total number of iterations across all staggered iterations
    index_t numIt_pf = 0, totIt_pf = 0; //numIt: per staggered iteration, totIt: total number of iterations across all staggered iterations
    index_t numIt_ref = 0; // number of refinement iterations per load step
    index_t numIt_stag = 0; // number of staggered iterations per load step

    real_t gamma   = 0.5;
    real_t beta    = 0.25;

    real_t Rnorm, R0;
    real_t Unorm, U0;
    T dt;
    index_t basis_size_old, basis_size; 
    index_t num_dofs_tot, num_el_tot;
    basis_size = basis_size_old = mb.basis(0).size();

    std::ofstream csvTotalTimes;
    csvTotalTimes.open(outputdir+"/output_times.csv");
    csvTotalTimes << "TimeStep,"<< "RefIt,"<< "StagIt,"<< "NumDOFs,"<< "numEl,"<< "E_E,"<< "E_D,"<< "T_init,"<< "T_EL_A,"<<"T_EL_S,"<<"T_PF_A,"<<"T_PF_S,"<<"T_Mark,"<<"T_Ref,"<<"T_Proj,"<<"T_total\n";
    T tTime = 0.0;

    gsStopwatch totalStop;
    totalStop.restart();
    while (tcurr<=tend)
    {
        dt = tcurr - tcurr_old;
        numIt_ref = numIt_stag = 0;
        totIt_el = totIt_pf = 0;

        bool refined = false;
        T basis_size_ratio;
        T markedArea = 0., tmpArea = 0.;

        gsInfo<<"---------------------------------------------------------------------------------------------------------------------------\n";
        gsInfo<<"Load step "<<step<<": t = "<<tcurr<<"\n";
        gsInfo<<"---------------------------------------------------------------------------------------------------------------------------\n";
        
        index_t refIt = 0;
        while(true) // refinement iteration loop
        {
            gsInfo<<"---------------------------------------------------------------------------------------------------------------------------\n";
            gsInfo<<"Refinement iteration: "<<refIt<<"\n";
            gsInfo<<"---------------------------------------------------------------------------------------------------------------------------\n";

            iterationTime = 0;
            elAssemblyTime = elSolverTime = 0.0;
            pfAssemblyTime = pfSolverTime = 0.0;
            bigClock.restart();
            material.setParameter(2,damage); 

            // Construct assembler with the new basis 
            gsSolidAssembler<dim,T,gsLinearDegradedMaterial<T>> elAssembler(mp,mb,bc_u,&material);
            elAssembler.options().setReal("ExprAssembler.quA",1.0);
            elAssembler.options().setInt ("ExprAssembler.quB",1);
            elAssembler.options().setInt("ExprAssembler.DirichletValues",dirichlet::l2Projection);
            elAssembler.initialize();
            // get new Mass matrix with the new mesh
            elAssembler.assembleMass();
            M = elAssembler.matrix();

            gsMatrix<> D_new, D_old, delta_D;
            gsMatrix<> u_new,u_old,udot_new,udot_old,uddot_new,uddot_old,delta_u;

            u_old.setZero(elAssembler.numDofs(),1);      
            udot_old.setZero(elAssembler.numDofs(),1);
            uddot_old.setZero(elAssembler.numDofs(),1);

            elAssembler.constructSolution(displacement_old,u_old);
            elAssembler.constructSolution(velocity_old,udot_old);
            elAssembler.constructSolution(acceleration_old,uddot_old);

            gsInfo<<"Phase-field assembler re-initialization\n";

            // ================== Initialize the phase-field assembler ==================
            gsPhaseFieldAssemblerBase<T> * pfAssembler;
            if      (order == 2 && AT == 1)
            {
                pfAssembler = new gsPhaseFieldAssembler<T,PForder::Second,PFmode::AT1>(mp,mb,bc_d);
                pfAssembler->options().setReal("cw",3.1615);
                pfAssembler->options().setReal("chi",1);
            }
            else if (order == 4 && AT == 1)
            {
                pfAssembler = new gsPhaseFieldAssembler<T,PForder::Fourth,PFmode::AT1>(mp,mb,bc_d);
                pfAssembler->options().setReal("cw",3.1615);
                pfAssembler->options().setReal("chi",0.0625);
            }
            else if (order == 2 && AT == 2)
            {
                pfAssembler = new gsPhaseFieldAssembler<T,PForder::Second,PFmode::AT2>(mp,mb,bc_d);
                pfAssembler->options().setReal("cw",3.1615);
                pfAssembler->options().setReal("chi",1);
            }
            else if (order == 4 && AT == 2)
            {
                pfAssembler = new gsPhaseFieldAssembler<T,PForder::Fourth,PFmode::AT2>(mp,mb,bc_d);
                pfAssembler->options().setReal("cw",4.4485);
                pfAssembler->options().setReal("chi",1);
            }
            else
                GISMO_ERROR("Invalid order and/or AT model");

            pfAssembler->options().setReal("l0",l0);
            pfAssembler->options().setReal("Gc",Gc);
            pfAssembler->initialize();
            pfAssembler->constructSolution(damage_old,D_old);
            
            // Assemble QPhi and q for energy computation (needed even if staggered loop exits early)
            pfAssembler->assemblePhi();
            pfAssembler->matrix_into(QPhi);
            pfAssembler->rhs_into(q);
            // ===========================================================================

            delta_u.setZero();
            delta_D.setZero();

            // Predictor step 
            udot_new = udot_old;
            uddot_new = (gamma-1)/gamma * uddot_old;
            u_new = u_old + dt * udot_old + 0.5*math::pow(dt,2) * ((1-2*beta) * uddot_old + 2*beta * uddot_new);

            Unorm = u_new.norm();
            U0 = (Unorm > 0) ? Unorm : 1.0; // Avoid division by zero
            Rnorm = R0 = 1;
            D_new = D_old;

            iterationTime += bigClock.stop();
            
            index_t stagIt = 0;
            while(true) // Staggered scheme
            {
                
                gsInfo<<"    --------------------------Staggered Iteration: "<<PRINT(4)<<stagIt<<"--------------------------\n";
                gsInfo<<"\t"<<PRINT(20)<<"* Elasticity:"<<PRINT(6)<<"It."<<PRINT(18)<<"||R||"<<PRINT(18)<<"||R||/||R0||"<<PRINT(18)<<"||ΔU||/||U0||"<<PRINT(18)<<"||dA||/||A||"<<PRINT(18)<<"||U||"<<PRINT(18)<<"||V||"<<PRINT(18)<<"||A||"<<PRINT(20)<<PRINT(20)<<"cum. assembly [s]"<<PRINT(20)<<"cum. solver [s]"<<"\n";

                // elAssemblyTime = elSolverTime = 0.0;
                // pfAssemblyTime = pfSolverTime = 0.0;
                bigClock.restart();

                material.setParameter(2,damage); // needed here too because it changes within the staggered iteration loop
                elAssembler.initialize();

                // ================================================ ELASTICITY ==============================================
                smallClock.restart();
                elAssembler.assemble(u_new);
                elAssemblyTime += smallClock.stop();
                elAssembler.matrix_into(K);
                elAssembler.rhs_into(Fext);

                R = M * uddot_new + K * u_new - Fext;
                Rnorm = R.norm();
                if (stagIt == 0)
                    R0 = (Rnorm > 0) ? Rnorm : 1.0; // For exit criterion (eq. (18))
                    
                K+= 1/(beta*math::pow(dt,2)) * M; // Eq. (21) - Greco et al. 2025

                smallClock.restart();
                solver.compute(K);
                delta_u = solver.solve(-R);
                elSolverTime += smallClock.stop();

                u_new += delta_u;
                udot_new  = (gamma/beta/dt) * ( u_new - u_old ) + (1-gamma/beta) * udot_old + (dt*(1-gamma/(2*beta))) * uddot_old;
                uddot_new = (1/(beta*math::pow(dt,2))) * ( u_new - u_old ) - (1/(beta*dt)) * udot_old - (1/(2*beta) - 1) * uddot_old;
                Unorm = u_new.norm();
                T DeltaUnorm = delta_u.norm();

                gsInfo<<"\t"<<PRINT(20)<<""<<PRINT(6)<<stagIt<<PRINT(18)<<Rnorm<<PRINT(18)<<Rnorm/R0<<PRINT(18)<<DeltaUnorm/U0<<PRINT(18)<<(uddot_new-uddot_old).norm()/uddot_new.norm()<<PRINT(18)<<Unorm<<PRINT(18)<<udot_new.norm()<<PRINT(18)<<uddot_new.norm()<<PRINT(20)<<elAssemblyTime<<PRINT(20)<<elSolverTime<<"\n";

                // Recompute the residual for the staggered check
                smallClock.restart(); // do i need to initialize the assembler?
                elAssembler.assemble(u_new);
                elAssemblyTime += smallClock.stop();
                elAssembler.matrix_into(K);
                elAssembler.rhs_into(Fext);
                R = M * uddot_new + K * u_new - Fext;
                Rnorm = R.norm();

                gsInfo << "\n";
                gsInfo << "    ****************************************************************************\n";
                gsInfo << "    *                          Staggered check                                 *\n";
                gsInfo << "    * " 
                    << PRINT(18) << "||R||" 
                    << PRINT(18) << "||R||/||R0||" 
                    << PRINT(18) << "||dU||" 
                    << PRINT(18) << "||dU||/||U0||" << " *\n";
                gsInfo << "    * "
                    << PRINT(18) << Rnorm 
                    << PRINT(18) << Rnorm/R0  
                    << PRINT(18) << DeltaUnorm  
                    << PRINT(18) << DeltaUnorm/U0 << " *\n";
                gsInfo << "    ****************************************************************************\n\n";
    
                // Update the fields before breaking the staggered loop 
                // (otherwise it does not update the displacements if staggered converges in 1 iteration)
                smallClock.restart();
                elAssembler.constructSolution(u_new,displacement);
                elAssembler.constructSolution(udot_new,velocity);
                elAssembler.constructSolution(uddot_new,acceleration);
                elAssemblyTime += smallClock.stop();

                for (size_t p=0; p!=mp.nPatches(); ++p)
                    mp_def.patch(p).coefs() = mp.patch(p).coefs() + displacement.patch(p).coefs();

                // Initialize the function for the elastic energy
                gsMaterialEval<T,gsMaterialOutput::Psi> Psi(&material,mp,mp_def);
                energy_E = 0.5 * (u_new.transpose() * K * u_new).value();                    

                // if (Rnorm/R0 < 1e-5 && DeltaUnorm/U0 < 1e-4)
                if (Rnorm/R0 < tol && DeltaUnorm/U0 < tolEl)
                    break;
                else if (stagIt == maxIt-1)
                    GISMO_ERROR("Staggered iterations problem did not converge.");
                stagIt++;

                // ================================================ PHASE-FIELD ==============================================
                gsInfo<<"\t"<<PRINT(20)<<"* Phase-Field:"<<PRINT(6)<<"It."<<PRINT(18)<<"||R||"<<PRINT(18)<<"||ΔD||"<<PRINT(18)<<"||ΔD||/||D||"<<PRINT(20)<<"cum. assembly [s]"<<PRINT(20)<<"cum. solver [s]"<<"\n";

                smallClock.restart();
                pfAssembler->assemblePsi(Psi);
                pfAssemblyTime += smallClock.stop();
                pfAssembler->matrix_into(QPsi);
                pfAssembler->rhs_into(qpsi);
                if (qpsi.rows()==0) // qpsi is empty for AT2 models
                    qpsi = gsMatrix<T>::Zero(QPsi.rows(),1);

                // QPhi and q changes size if the mesh is refined 
                pfAssembler->initialize();
                pfAssembler->assemblePhi();
                pfAssembler->matrix_into(QPhi);
                pfAssembler->rhs_into(q);   
                
                Q = QPhi + QPsi;

                // Reconstruct the solution from the damage field
                pfAssembler->constructSolution(damage,D_new);

                // Initialize the PSOR solver
                smallClock.restart();
                gsPSOR<T> PSORsolver(Q);
                PSORsolver.options().setInt("MaxIterations",30000);
                PSORsolver.options().setSwitch("Verbose",false);
                PSORsolver.options().setReal("tolU",tolEl);
                PSORsolver.options().setReal("tolNeg",tolPf);
                PSORsolver.options().setReal("tolPos",tolPf);
                // gsInfo<< PSORsolver.options() << "\n";
                pfSolverTime = smallClock.stop();
                index_t pfIt = 0;
                while(true)
                {
                    // Assemble
                    smallClock.restart();
                    R = Q * D_new - qpsi + q;
                    pfAssemblyTime += smallClock.stop();

                    smallClock.restart();
                    PSORsolver.solve(R,delta_D); // delta_D = Q \ R
                    pfSolverTime += smallClock.stop();
                    D_new += delta_D;

                    gsInfo<<"\t"<<PRINT(20)<<" "<<PRINT(6)<<pfIt<<PRINT(18)<<R.norm()<<PRINT(18)<<delta_D.norm()<<PRINT(18)<<delta_D.norm()/D_new.norm()<<PRINT(20)<<pfAssemblyTime<<PRINT(20)<<pfSolverTime<<"\n";

                    // if (delta_D.norm()/D_new.norm() < tolPf || D_new.norm() < 1e-12 || maxItPf==1)
                    if (delta_D.norm()/D_new.norm() < tolPf || D_new.norm() < 1e-12 || maxItPf==1)
                        break;
                    else if (pfIt == maxItPf-1 && maxItPf != 1)
                        GISMO_ERROR("Phase-field problem did not converge.");
                    pfIt++;
                }
                numIt_pf = pfIt+1;
                totIt_pf+= numIt_pf;
                
                // Update damage spline
                pfAssembler->constructSolution(D_new,damage);

            } // end staggered loop
            numIt_stag += stagIt+1;

            energy_D = (0.5 * D_new.transpose() * QPhi * D_new).value() + (D_new.transpose() * q).value();

            // Save size of the basis and element size
            num_dofs_tot = (dim+1)*mb.basis(0).size();
            num_el_tot   = mb.basis(0).numElements();

            // Mesh refinement 
            if (adaptive_switch)
            {
                for (index_t i=0; i!=maxRefLvl; ++i)
                {
                    smallClock.restart();
                    elVals = labelElements<dim,T>(mp, damage, mb,0.1,1.0);
                    labelTime = smallClock.stop();
                    gsInfo<<"Labelling level "<<i<<" took "<<smallClock.stop()<<" seconds\n";
                    if (gsAsVector<T>(elVals).sum() > 0)
                    {
                        smallClock.restart();
                        tmpArea = refineMesh<dim,T>(mb,elVals,mesherOptions);
                        refTime = smallClock.stop();
                        gsInfo<<"Refining mesh took "<<smallClock.stop()<<" seconds\n";
                    }

                    tmpArea /= (mb.basis(0).support().col(1)-mb.basis(0).support().col(0)).prod();
                    markedArea = math::max(markedArea,tmpArea);
                    basis_size = mb.basis(0).size();
                    refined = basis_size > basis_size_old;
                    if (!refined)
                        break;
                }
            }
            // gsInfo<<"Marked area: "<<markedArea<<"\n";
            // gsInfo<<"refined "<< refined<<"\n";

            ///????
            // refined &= markedArea > mesherOptions.askReal("SizeRatio",1.01);

            if (refined)
            {
                gsInfo << "The basis HAS BEEN REFINED: "<< basis_size_old << "=>"<< basis_size<<"\n";
                // ======================================================================================
                // Project the new and the old solution onto the new mesh (we move to the next time step)
                // ======================================================================================
                smallClock.restart();
                gsMatrix<T> projCoefs;
                // Geometry (undeformed)
                gsQuasiInterpolate<T>::localIntpl(mb.basis(0),mp.patch(0),projCoefs);
                mp.clear();
                mp.addPatch(mb.basis(0).makeGeometry(give(projCoefs)));
                // Geometry (deformed)
                gsQuasiInterpolate<T>::localIntpl(mb.basis(0),mp_def.patch(0),projCoefs);
                mp_def.clear();
                mp_def.addPatch(mb.basis(0).makeGeometry(give(projCoefs)));
                // Displacement (new)
                gsQuasiInterpolate<T>::localIntpl(mb.basis(0),displacement.patch(0),projCoefs);
                displacement.clear();
                displacement.addPatch(mb.basis(0).makeGeometry(give(projCoefs)));
                // Displacement (old)
                gsQuasiInterpolate<T>::localIntpl(mb.basis(0),displacement_old.patch(0),projCoefs);
                displacement_old.clear();
                displacement_old.addPatch(mb.basis(0).makeGeometry(give(projCoefs)));
                // ========================================================================
                // New variables for dynamics
                // Velocity (new)
                gsQuasiInterpolate<T>::localIntpl(mb.basis(0),velocity.patch(0),projCoefs);
                velocity.clear();
                velocity.addPatch(mb.basis(0).makeGeometry(give(projCoefs)));
                // Velocity (old)
                gsQuasiInterpolate<T>::localIntpl(mb.basis(0),velocity_old.patch(0),projCoefs);
                velocity_old.clear();
                velocity_old.addPatch(mb.basis(0).makeGeometry(give(projCoefs)));
                // Acceleration (new)
                gsQuasiInterpolate<T>::localIntpl(mb.basis(0),acceleration.patch(0),projCoefs);
                acceleration.clear();
                acceleration.addPatch(mb.basis(0).makeGeometry(give(projCoefs)));
                // Acceleration (old)
                gsQuasiInterpolate<T>::localIntpl(mb.basis(0),acceleration_old.patch(0),projCoefs);
                acceleration_old.clear();
                acceleration_old.addPatch(mb.basis(0).makeGeometry(give(projCoefs)));
                // ========================================================================
                // Damage
                gsQuasiInterpolate<T>::localIntpl(mb.basis(0),damage.patch(0),projCoefs);
                damage.clear();
                damage.addPatch(mb.basis(0).makeGeometry(give(projCoefs)));
                // Damage (old)
                gsQuasiInterpolate<T>::localIntpl(mb.basis(0),damage_old.patch(0),projCoefs);
                damage_old.clear();
                damage_old.addPatch(mb.basis(0).makeGeometry(give(projCoefs)));
                projTime = smallClock.stop();
                stepTimes.projectionTime += projTime;
                gsInfo<<"Projection took "<<projTime<<" seconds\n";

                // Update the basis size
                basis_size_old = basis_size;
            }
            else
            {
                gsInfo << "The basis has NOT been refined.\n";
                // break;
            }

            
            // gsInfo<<"antes del csv time TIME: "<<pfSolverTime<< "\n";
            // Update csv file data
            totalTime = iterationTime + elAssemblyTime + elSolverTime + pfAssemblyTime + pfSolverTime + labelTime + refTime + projTime;
            csvTotalTimes << step  << "," << refIt << ","<< stagIt+1 <<"," << num_dofs_tot << "," << num_el_tot <<"," << energy_E<<"," << energy_D << "," << iterationTime <<","<< elAssemblyTime  <<","<< elSolverTime<<","<< pfAssemblyTime << ","<< pfSolverTime << "," << labelTime << ","<< refTime << ","<<projTime << ","<< totalTime <<"\n";
            csvTotalTimes.flush(); 
            
            if (!refined) // to make sure it writes the results of the last refinement iteration
                break;

            refIt++;
        } // end refinement loop
        numIt_ref = refIt+1;

        // =========================================================================
        // // Compute resulting force and energies (i need to check the size of the assembler!)
        gsSolidAssembler<dim,T,gsLinearDegradedMaterial<T>> fullElAssembler(mp,mb,bc_u_dummy,&material);
        fullElAssembler.options().setReal("ExprAssembler.quA",1.0);
        fullElAssembler.options().setInt ("ExprAssembler.quB",1);
        // fullElAssembler.options().setInt("ExprAssembler.DirichletValues",dirichlet::l2Projection);
        fullElAssembler.initialize();
        gsMatrix<T> ufull = displacement.patch(0).coefs().reshape(displacement.patch(0).coefs().size(),1);
        fullElAssembler.assemble(ufull);
        gsMatrix<T> Rfull = fullElAssembler.rhs();
        // sum the reaction forces in Y direction
        gsDofMapper mapper(mb,dim);
        mapper.finalize();
        gsMatrix<index_t> boundary = mb.basis(0).boundary(boundary::east);
        T Fx = 0, Fy = 0;
        // for (index_t k=0; k!=boundary.size(); k++)
        // {
        //     Fx += Rfull(mapper.index(boundary(k,0),0,0),0); // DoF index, patch, component
        //     Fy += Rfull(mapper.index(boundary(k,0),0,1),0); // DoF index, patch, component
        // }

        std::vector<T> stepData(5);
        stepData[0] = step;
        stepData[1] = Fx;
        stepData[2] = Fy;
        stepData[3] = (0.5 * ufull.transpose() * fullElAssembler.matrix() * ufull).value();
        stepData[4] = energy_D;

        gsInfo<<"\n";
        gsInfo<<"Converged with ||R||/||R0|| = "<<Rnorm/R0<<" < "<<tol<<" ||dD||= "<<delta_D.norm()<<" ||dU||/U0 = "<<delta_u.norm()<<"\n";
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
            np[0] = 400;
            np[1] = 400;
            // if (dim==2)
            //     np[1] = 2;
            // else
            //     np[1] = np[2] = 2;
            pts = gsPointGrid(a,b,np);
            // Evaluate the geometry
            eval_geo = mp_def.patch(0).eval(pts);

            std::string filename;
            std::string subfolder = outputdir + "damage_pvd/";
            gsFileManager::mkdir(subfolder);
            filename = "damage_pvd/damage_" + util::to_string(step);
            eval_damage = damage.patch(0).eval(pts);
            gsWriteParaviewTPgrid(eval_geo,eval_damage,np.template cast<index_t>(),outputdir+filename);
            // gsWriteParaview(mp,damage,outputdir+filename,100000);
            damageCollection.addPart(filename,step,"Solution",0);

            subfolder.clear();
            filename.clear();
            subfolder = outputdir + "Psi_pvd/";
            gsFileManager::mkdir(subfolder);
            gsMaterialEval<T,gsMaterialOutput::Psi,true,true> Psi(&material,mp,mp_def);
            filename = "Psi_pvd/Psi_"+util::to_string(step);
            eval_psi = Psi.piece(0).eval(pts);
            gsWriteParaviewTPgrid(eval_geo,eval_psi,np.template cast<index_t>(),outputdir+filename);
            // gsWriteParaview(mp,Psi,outputdir+filename,100000);
            psiCollection.addPart(filename,step,"Solution",0);

            subfolder.clear();
            filename.clear();
            subfolder = outputdir + "displacement_pvd/";
            gsFileManager::mkdir(subfolder);
            filename = "displacement_pvd/displacement_"+util::to_string(step);
            eval_displacement = displacement.patch(0).eval(pts);
            gsInfo<<"MAX PARAVIEW: "<<eval_displacement.maxCoeff()<<"\n";
            gsWriteParaviewTPgrid(eval_geo,eval_displacement,np.template cast<index_t>(),outputdir+filename);
            // gsWriteParaview(mp,displacement,outputdir+filename,1000);
            displCollection.addPart(filename,step,"Solution",0);
            

            gsInfo << "plotMesh = " << plotMesh << "\n";
            gsInfo << "Output dir: " << outputdir << "\n";
            // gsInfo << "Creating: " << subfolder << "\n";
            // gsInfo << "mkdir returned: " << ok << "\n";


            // Plot mesh
            if (plotMesh) // to be polished
            {
                subfolder.clear();
                filename.clear();
                subfolder = outputdir + "mesh_pvd/";
                gsFileManager::mkdir(subfolder);
                filename = "mesh_pvd/mesh_"+util::to_string(step);
                gsWriteParaview(mp, outputdir + filename, 10, true); // (creates a pvd file at every time step...)
                meshCollection.addPart(filename + "_0_mesh.vtp", step, "Mesh", 0);
                // delete auto-generated volume grid and per-step PVD
                std::remove((outputdir + filename + "_0.vts").c_str());
                gsInfo<< outputdir + filename + "_0.pvd" <<"\n";
                std::remove((outputdir + filename + ".pvd").c_str());
            }
        }

        // =========================================================================
        // Write data
        file.open(outputdir+"results.txt",std::ios::app);
        // for (size_t i = 0; i != data.size(); ++i)
        //     file<<data[i][0]<<","<<-data[i][1]<<","<<-data[i][2]<<","<<data[i][3]<<","<<data[i][4]<<"\n";
        file<<stepData[0]<<","<<-stepData[1]<<","<<-stepData[2]<<","<<stepData[3]<<","<<stepData[4]<<"\n";
        file.close();

        // u_old = u_new;
        // udot_old = udot_new;
        // uddot_old = uddot_new;
        // D_old = D_new;

        displacement_old    = displacement;
        velocity_old        = velocity;
        acceleration_old    = acceleration;
        damage_old          = damage;    
 
        tcurr_old = tcurr;
        tcurr += (tcurr+tstep > ttrans) ? tstep/tred : tstep;
        step++;
    }

    if (plot)
    {
        damageCollection.save();
        psiCollection.save();
        displCollection.save();
        PwaveCollection.save();
        meshCollection.save();
    }
    
    tTime += totalStop.stop();
    csvTotalTimes<< "The total simulation time is: " << tTime<<"\n";
    csvTotalTimes.close();

    delete pfAssembler;
}