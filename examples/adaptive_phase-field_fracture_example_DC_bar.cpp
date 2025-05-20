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
    gsMultiPatch<> mp_ini;
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
        mp_ini.addPatch(tb);
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
        mp_ini.addPatch(tb);
    }
    else
        GISMO_ERROR("Invalid dimension");

    if (plot) gsWriteParaview(mp_ini,outputdir+"mp_ini",10,true);

    if (numElev > 0)
        mp_ini.degreeIncrease(numElev);
    for (index_t i = 0; i<numHRef; ++i)
        mp_ini.uniformRefine(1);

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
    controlParameters.addReal("tolEl", "Tolerance for the elasticity problem", 1e-9);
    // Tolerance for the phase-field problem
    controlParameters.addReal("tolPf", "Tolerance for the phase-field problem", 1e-9);
    // Staggered tolerance
    controlParameters.addReal("tol", "Tolerance for the staggered scheme", 1e-6);

    // Mesher options
    gsOptionList mesherOptions;
    // Admissible meshing
    mesherOptions.addSwitch("Admissible", "Admissible meshing",true);
    // Maximum level
    mesherOptions.addInt("MaxLevel", "Maximum level of refinement", 5);
    // Refinement rule
    mesherOptions.addInt("RefineRule", "Refinement rule", 1);
    // Coarsening rule
    mesherOptions.addInt("CoarsenRule", "Coarsening rule", 1);
    // Refinement parameter
    mesherOptions.addReal("RefineParam", "Refinement parameter", 0.1);
    // Coarsening parameter
    mesherOptions.addReal("CoarsenParam", "Coarsening parameter", 0.1);

    // Initialize the damage field
    gsMultiPatch<> damage;
    gsMatrix<> coefs(mp_ini.basis(0).size(),1);
    coefs.setZero();

    damage.addPatch(mp_ini.basis(0).makeGeometry(give(coefs)));

    // Boundary conditions
    gsBoundaryConditions<> bc_u;
    bc_u.setGeoMap(mp_ini);

    gsBoundaryConditions<> bc_d;
    // bc_d.addCondition(boundary::west,condition_type::dirichlet,0,0);
    // bc_d.addCondition(boundary::east,condition_type::dirichlet,0,0);
    bc_d.setGeoMap(mp_ini);

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

template<short_t dim, class T>
std::vector<T> labelElements(  const gsMultiPatch<> & geometry,
                                    const gsFunctionSet<>& damage,
                                    const gsMultiBasis<> & basis,
                                    const T         & lowerBound,
                                    const T         & upperBound)
{
    GISMO_ASSERT(basis.nBases() == 1, "Labeling is only implemented for single basis meshes");
    typename gsBasis<T>::domainIter domIt  = basis.basis(0).domain()->beginAll();
    typename gsBasis<T>::domainIter domEnd = basis.basis(0).domain()->endAll();
    gsMatrix<T> points;
    gsMatrix<T> vals;
    std::vector<T> labels(basis.basis(0).numElements());
    gsVector<unsigned,dim> np;
    np.setConstant(2);
    for (; domIt<domEnd; ++domIt)
    {
        points = gsPointGrid(domIt.lowerCorner(),domIt.upperCorner(),np);
        damage.piece(0).eval_into(points,vals);
        labels[domIt.id()] = (vals.array() >= lowerBound && vals.array() <= upperBound).any();
    }
    return labels;
}

template <short_t dim, class T>
void refineMesh    (      gsMultiBasis<>      & basis,
                    const std::vector<T> & vals,
                    const gsOptionList        & options)
    {
    gsAdaptiveMeshing<dim,T> mesher(basis);
    mesher.options().setSwitch("Admissible",true);
    mesher.options().setInt("MaxLevel",1);
    mesher.options().setInt("RefineRule",1);
    mesher.options().setInt("CoarsenRule",1);
    mesher.options().setReal("RefineParam",0.1);
    mesher.options().setReal("CoarsenParam",0.1);
    mesher.options().update(options,gsOptionList::ignoreIfUnknown);
    mesher.getOptions();

    gsHBoxContainer<dim,T> refine, coarsen;
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

    // Convert to THB
    gsMultiPatch<> mp;
    for (index_t i = 0; i < mp_ini.nPatches(); ++i)
    {
        // Check if tensor basis
        if      ((dynamic_cast<const gsTensorBSpline<2,real_t> *>(&mp_ini.patch(i))))
        {
            // Create a THB spline basis
            const gsTensorBSpline<dim,real_t> & tb = static_cast<const gsTensorBSpline<dim,real_t> &>(mp_ini.patch(i));
            gsTHBSpline<dim,real_t> thb(tb);
            mp.addPatch(memory::make_unique(thb.clone().release()));
        }
        else if ((dynamic_cast<const gsTHBSpline<dim,real_t> *>(&mp_ini.patch(i))))
        {
            const gsTHBSpline<dim,real_t> & thb = static_cast<const gsTHBSpline<dim,real_t> &>(mp_ini.patch(i));
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
    gsMultiPatch<> mp_def = mp;
    gsMultiPatch<> displacement = mp;
    for (index_t p=0; p<mp.nPatches(); ++p)
        displacement.patch(p).coefs().setZero();

    // Initialize the material
    gsLinearDegradedMaterial<T> material(E,nu,damage,dim);
    // Initialize the elasticity assembler
    gsConstantFunction<> bodyForce(0.,0.,dim);

    gsBoundaryConditions<> bc_u_dummy;
    bc_u_dummy.setGeoMap(mp);

    // Initialize the phase-field assembler
    gsPhaseFieldAssemblerBase<real_t> * pfAssembler = NULL;

    //////////////////////////////////////////////////////////////////////////
    // INITIALIZE THE MESH ///////////////////////////////////////////////////
    //////////////////////////////////////////////////////////////////////////

    std::vector<real_t> elVals;
    // REFINE MESH
    for (index_t k=0; k!=mesherOptions.askInt("MaxLevel",1); ++k)
    {
        elVals.resize(mb.basis(0).numElements());
        // Loop over the elements
        auto domIt = mb.basis(0).domain()->beginAll();
        auto domEnd = mb.basis(0).domain()->endAll();
        for (; domIt<domEnd; ++domIt)
        {
            if ((domIt.lowerCorner()(0,0) >= 0.5-1e-6 && domIt.lowerCorner()(0,0) <= 0.5+1e-6)
                ||
                (domIt.upperCorner()(0,0) >= 0.5-1e-6 && domIt.upperCorner()(0,0) <= 0.5+1e-6)
                )
                elVals[domIt.id()] = 1;
            else
                elVals[domIt.id()] = 0;
        }
        if (gsAsVector<real_t>(elVals).sum())
            refineMesh<dim,T>(mb,elVals,mesherOptions);
        else
            break;
    }
    writeSingleCompMesh(mb.basis(0),mp.patch(0),outputdir+"initial_mesh",1);

    //////////////////////////////////////////////////////////////////////////
    // PERTURB THE DAMAGE/////////////////////////////////////////////////////
    //////////////////////////////////////////////////////////////////////////

    gsMatrix<> anchors = mb.basis(0).anchors();
    gsMatrix<> coefs(mb.basis(0).size(),1);
    coefs.setZero();
    // Find the anchors closest to the middle of the domain
    // and set the damage to 1e-4
    gsMatrix<> tmp = anchors;
    gsVector<> mid(dim);
    gsMatrix<> norms;
    mid.setConstant(0.5);
    tmp.colwise() -= mid;
    norms = tmp.colwise().norm();
    real_t min = norms.minCoeff();
    gsDebugVar(norms);
    gsDebugVar(min);
    for (index_t k=0; k<coefs.rows(); ++k)
        if (gsClose(norms(0,k),min,1e-10))
        {
            coefs(k,0) = 1e-4;
            gsInfo<<"Damage initialized at "<<anchors(0,k)<<" "<<anchors(0,k)<<"\n";
        }

    damage.clear();
    damage.addPatch(mb.basis(0).makeGeometry(give(coefs)));

    //////////////////////////////////////////////////////////////////////////
    // SOLVE THE PROBLEM/////////////////////////////////////////////////////
    //////////////////////////////////////////////////////////////////////////

    gsMatrix<T> u;

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
    file<<"u,Fx,Fy,E_u,E_d,basis_size,num_elements\n";
    file.close();

    std::vector<gsMatrix<> > fixedDofs;
    std::vector<gsMatrix<> > dummyFixedDofs;
    real_t Rnorm, Fnorm;
    Rnorm = Fnorm = 1;
    bool refined = true;
    while (ucurr<=uend)
    {
        if (refined)
        {
            // PROJECT SOLUTIONS
            gsMatrix<> projCoefs;
            gsQuasiInterpolate<real_t>::localIntpl(mb.basis(0),mp.patch(0),projCoefs);
            mp.clear();
            mp.addPatch(mb.basis(0).makeGeometry(give(projCoefs)));
            mp_def = mp;

            gsQuasiInterpolate<real_t>::localIntpl(mb.basis(0),displacement.patch(0),projCoefs);
            displacement.clear();
            displacement.addPatch(mb.basis(0).makeGeometry(give(projCoefs)));

            gsQuasiInterpolate<real_t>::localIntpl(mb.basis(0),damage.patch(0),projCoefs);
            damage.clear();
            damage.addPatch(mb.basis(0).makeGeometry(give(projCoefs)));
        }
        // CONSTRUCT ASSEMBLERS (since refreshing is not possible)
        gsElasticityAssembler<real_t> elAssembler(mp,mb,bc_u,bodyForce,&material);
        elAssembler.options().setReal("quA",1.0);
        elAssembler.options().setInt ("quB",0);
        elAssembler.options().setSwitch ("SmallStrains",true);
        fixedDofs = elAssembler.allFixedDofs();
        elAssembler.constructSolution(displacement,u);

        gsElasticityAssembler<real_t> fullElAssembler(mp,mb,bc_u_dummy,bodyForce,&material);
        fullElAssembler.options().setReal("quA",1.0);
        fullElAssembler.options().setInt ("quB",0);
        fullElAssembler.options().setSwitch ("SmallStrains",true);
        dummyFixedDofs = fullElAssembler.allFixedDofs();

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
        displ_left.set_u(ucurr);
        displ_right.set_u(ucurr);
        elAssembler.computeDirichletDofs(0); // NOTE: This computes the DDofs for **unknown** 1, which should be component 1. This is a bug in the gsElasticity assembler
        elAssembler.computeDirichletDofs(1); // NOTE: This computes the DDofs for **unknown** 1, which should be component 1. This is a bug in the gsElasticity assembler
        fixedDofs = elAssembler.allFixedDofs();
        elAssembler.setFixedDofs(fixedDofs);

        gsInfo<<"---------------------------------------------------------------------------------------------------------------------------\n";
        gsInfo<<"Load step "<<step<<": u = "<<ucurr<<"\n\n";

        deltaD.setZero(D.rows(),1);

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
            gsMaterialEval<T,gsMaterialOutput::Psi,true,true> Psi(&material,mp,mp_def);

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
            if (Rnorm/* /Fnorm */ < tol)
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

        std::vector<real_t> stepData(7);
        stepData[0] = ucurr;
        stepData[1] = Fx;
        stepData[2] = Fy;
        stepData[3] = (0.5 * ufull.transpose() * fullElAssembler.matrix() * ufull).value();
        stepData[4] = (0.5 * D.transpose() * QPhi * D).value() + (D.transpose() * q).value();
        stepData[5] = mb.totalSize();
        stepData[6] = mb.totalElements();

        gsInfo<<"\n";
        gsInfo<<"Converged with ||R||/||F|| = "<<Rnorm/Fnorm<<" < "<<tol<<" ||D|| = "<<D.norm()<<" ||U|| = "<<u.norm()<<"\n";
        // gsInfo<<"----------------------------------------------------------------------------------------------------\n\n";

        // =========================================================================
        // PLOT
        if ((plot && step%plotmod==0) || refined)
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
            filename = "mesh_" + util::to_string(step);
            gsMesh<> mesh(mb.basis(0));
            writeSingleCompMesh(mb.basis(0),mp.patch(0),outputdir+filename,1);
            meshCollection.addPart(filename,step,"Mesh",0);

            filename = "damage_" + util::to_string(step);
            eval_damage = damage.patch(0).eval(pts);
            gsWriteParaviewTPgrid(eval_geo,eval_damage,np.template cast<index_t>(),outputdir+filename);
            // gsWriteParaview(mp,damage,outputdir+filename,100000);
            damageCollection.addPart(filename,step,"Solution",0);

            gsMaterialEval<T,gsMaterialOutput::Psi> Psi(&material,mp,mp_def);
            filename = "Psi_"+util::to_string(step);
            eval_psi = Psi.piece(0).eval(pts);
            gsWriteParaviewTPgrid(eval_geo,eval_psi,np.template cast<index_t>(),outputdir+filename);
            // gsWriteParaview(mp,Psi,outputdir+filename,100000);
            psiCollection.addPart(filename,step,"Solution",0);

            filename = "displacement_"+util::to_string(step);
            eval_displacement = displacement.patch(0).eval(pts);
            gsWriteParaviewTPgrid(eval_geo,eval_displacement,np.template cast<index_t>(),outputdir+filename);
            // gsWriteParaview(mp,displacement,outputdir+filename,1000);
            displCollection.addPart(filename,step,"Solution",0);
        }

        // =========================================================================
        // Write data
        std::ofstream file(outputdir+"results.txt",std::ios::app);
        // for (size_t i = 0; i != data.size(); ++i)
        //     file<<data[i][0]<<","<<-data[i][1]<<","<<-data[i][2]<<","<<data[i][3]<<","<<data[i][4]<<"\n";
        file<<stepData[0]<<","<<-stepData[1]<<","<<-stepData[2]<<","<<stepData[3]<<","<<stepData[4]<<","<<stepData[5]<<","<<stepData[6]<<"\n";
        file.close();

        // =========================================================================
        // REFINE MESH
        elVals = labelElements<dim,T>(mp, damage, mb,0.1,1.0);
        refined = gsAsVector<real_t>(elVals).sum() > 0;
        if (refined)
            refineMesh<dim,T>(mb,elVals,mesherOptions);

            // =========================================================================
        // INCREMENT STEP

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