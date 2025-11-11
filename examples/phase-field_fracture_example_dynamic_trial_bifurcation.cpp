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
    index_t numElev = 0;
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

    real_t L = 100.;
    real_t H = 40.;

    if (dimension == 2)
    {
        index_t nx = 99; // elements in x direction
        index_t ny = 39;  // elements in y direction
        gsKnotVector<> kv_x(0, 1, nx, 2, 1); 
        gsKnotVector<> kv_y(0, 1, ny, 2, 1);
        gsTensorBSplineBasis<2,real_t> tbasis(kv_x, kv_y);
        gsTensorBSpline<2,real_t> tb(tbasis, tbasis.anchors().transpose());
        tb.coefs().col(0) *= L;
        tb.coefs().col(1) *= H;
        mp.addPatch(tb);
    }

    // if (dimension == 2)
    // {
    //     gsKnotVector<> kv_x(0,1,math::ceil(L/H)-1,2,1);
    //     gsKnotVector<> kv_y(0,1,math::ceil(H/L)-1,2,1);
    //     gsTensorBSplineBasis<2,real_t> tbasis(kv_x, kv_y);
    //     gsTensorBSpline<2,real_t> tb(tbasis,tbasis.anchors().transpose());
    //     tb.coefs().col(0) *= L;
    //     tb.coefs().col(1) *= H;
    //     mp.addPatch(tb);
    // }
    // else if (dimension == 3)
    // {
    //     gsKnotVector<> kv_x(0,1,math::ceil(L/H)-1,2,1);
    //     gsKnotVector<> kv_y(0,1,math::ceil(H/L)-1,2,1);
    //     gsKnotVector<> kv_z(0,1,math::ceil(H/L)-1,2,1);
    //     gsTensorBSplineBasis<3,real_t> tbasis(kv_x, kv_y, kv_z);
    //     gsTensorBSpline<3,real_t> tb(tbasis,tbasis.anchors().transpose());
    //     tb.coefs().col(0) *= L;
    //     tb.coefs().col(1) *= H;
    //     tb.coefs().col(2) *= H;
    //     mp.addPatch(tb);
    // }
    // else
    //     GISMO_ERROR("Invalid dimension");

    if (plot) gsWriteParaview(mp,outputdir+"mp",10,true);

    // mp.degreeIncrease(numElev);
    // for (index_t i = 0; i<numHRef; ++i)
    //     mp.uniformRefine();

    //// Material parameters 
    gsOptionList materialParameters;
    // Young's modulus [kN/mm^2]
    materialParameters.addReal("E", "Young's modulus", 10e3);
    // Poisson's ratio [-]
    materialParameters.addReal("nu", "Poisson's ratio", 0.);
    // Toughness [kN/mm]
    materialParameters.addReal("Gc", "Toughness", 1e-3); // corrected
    // Internal length [mm]
    materialParameters.addReal("l0", "Internal length", 0.5);
    // Density [tonn/mm^3]
    materialParameters.addReal("rho", "Density", 2.46e-9);
    // Order of the phase-field model
    materialParameters.addInt("order", "Order of the phase-field model", order);
    // AT model
    materialParameters.addInt("AT", "AT model", AT);

    //// Boundary control parameters
    gsOptionList controlParameters;
    // Min time [s]
    controlParameters.addReal("tend", "Maximum time", 3.1e-06);
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
    // bc_d.addCondition(boundary::west,condition_type::dirichlet,0,0,false,0);
    // bc_d.addCondition(boundary::east,condition_type::dirichlet,0,0,false,0);
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
    // Load parameters (paper Borden et al. 2012)
    ////////////////////////////////////////////////////////////////////////////////////
    // Young's modulus [N/mm^2]
    T E = 32e3;
    // Poisson's ratio [-]
    T nu = 0.2;
    // Toughness [N/mm]
    T Gc = 3e-3;
    // Internal length [mm]
    T l0 = 0.5;
    // Density [tonn/mm^3]
    T rho = 2.45e-6;
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

    // Construct the basis
    gsMultiBasis<T> mb(mp);
    gsInfo<<"The basis has size "<<mb.size()<<" and degree "<<mb.degree()<<"\n";
    for (size_t b=0; b!=mb.nBases(); b++)
        gsInfo<<"Basis "<<b<<":\n"<<mb.basis(b)<<"\n";

    // Boundary conditions
    T sigma = 1.0;
    std::vector<std::string> bcFunctionTop(dim, "0");
    std::vector<std::string> bcFunctionBottom(dim, "0");
    bcFunctionTop[1] = "u";
    bcFunctionBottom[1] = "-u";
    gsFunctionExpr<T> sigma_top(bcFunctionTop,dim);
    gsFunctionExpr<T> sigma_bottom(bcFunctionBottom,dim);
    sigma_top.set_u(sigma);
    sigma_bottom.set_u(sigma);
    gsInfo<< sigma_top << "\n";
    gsInfo<< sigma_bottom << "\n";
    gsVector<> pt(2);
    pt<<100,40;
    gsInfo<<"sigma top at pt: "<<pt[0]<<","<<pt[1]<<":"<<sigma_top.eval(pt)<<"\n";
    // gsInfo<<"sigma top at "<< pt: "<<sigma_top.eval(pt)<<"\n";
    bc_u.addCondition(boundary::north,condition_type::neumann,&sigma_top );
    bc_u.addCondition(boundary::south,condition_type::neumann,&sigma_bottom);
    // bc_u.addCondition(boundary::south,condition_type::dirichlet,0,0,false,1); //vertical constraint
    // bc_u.addCondition(boundary::north,condition_type::dirichlet,0,0,false,1); //vertical constraint
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
    gsMultiPatch<T> velocity;
    for (index_t p=0; p<mp.nPatches(); ++p)
        displacement.patch(p).coefs().setZero();

    // Initialize the material
    gsLinearDegradedMaterial<T> material(E,nu,rho,damage,dim);

    // Initialize the elasticity assembler
    // Initialize the elasticity assembler
    gsSolidAssembler<dim,T,gsLinearDegradedMaterial<T>> elAssembler(mp,mb,bc_u,&material);
    elAssembler.options().setReal("ExprAssembler.quA",1.0);
    elAssembler.options().setInt ("ExprAssembler.quB",0);
    elAssembler.initialize();
    elAssembler.assemble();

    gsDebugVar(elAssembler.rhs().norm());



    gsBoundaryConditions<T> bc_u_dummy;
    bc_u_dummy.setGeoMap(mp);
    gsSolidAssembler<dim,T,gsLinearDegradedMaterial<T>> fullElAssembler(mp,mb,bc_u_dummy,&material);
    fullElAssembler.options().setReal("ExprAssembler.quA",1.0);
    fullElAssembler.options().setInt ("ExprAssembler.quB",0);
    fullElAssembler.initialize();

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
    pfAssembler->initialize();

    //////////////////////////////////////////////////////////////////////////
    // SOLVE THE PROBLEM/////////////////////////////////////////////////////
    //////////////////////////////////////////////////////////////////////////

    gsMatrix<> u_new(elAssembler.numDofs(),1),
               u_old(elAssembler.numDofs(),1),
               udot_new(elAssembler.numDofs(),1),
               udot_old(elAssembler.numDofs(),1),
               uddot_new(elAssembler.numDofs(),1),
               uddot_old(elAssembler.numDofs(),1),
               delta_u;
    u_new.setZero();
    u_old.setZero();
    udot_new.setZero();
    udot_old.setZero();
    uddot_new.setZero();
    uddot_old.setZero();

    gsMatrix<> D_new(pfAssembler->numDofs(),1),
               D_old(pfAssembler->numDofs(),1),
               delta_D(pfAssembler->numDofs(),1);

#ifdef GISMO_WITH_PARDISO
    typename gsSparseSolver<T>::PardisoLDLT solver;
#else
    typename gsSparseSolver<T>::CGDiagonal solver;
#endif

    gsSparseMatrix<T> K;
    gsMatrix<T> R, Fext;

    T elAssemblyTime = 0.0;
    T elSolverTime = 0.0;
    T pfAssemblyTime = 0.0;
    T pfSolverTime = 0.0;
    T iterationTime  = 0.0;

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
    gsStopwatch smallClock, bigClock;

    // pfAssembler->assembleMatrix();
    // pfAssembler->matrix_into(QPhi);
    // pfAssembler->constructSolution(damage,D);
    // gsInfo<<"D_0 = "<<(0.5 * D.transpose() * QPhi * D).value()<<"\n";
    //     gsWriteParaview(mp,damage,"damage_ini",100000);

    std::ofstream file(outputdir+"results.txt");
    file<<"u,Fx,Fy,E_u,E_d\n";
    file.close();

    /* @todo: add option list from file for dynamic parameters */
    // gsOptionList dynamicParameters;
    real_t gamma   = 0.5;
    real_t beta    = 0.25;

    real_t Rnorm, R0;
    real_t Unorm, U0;
    T dt;
    while (tcurr<=tend)
    {
        dt = tcurr - tcurr_old;
        // Update the boundary conditions
        // displ_left.set_u(tcurr);
        // displ_right.set_u(tcurr);
        // elAssembler.computeDirichletDofs(0); // NOTE: This computes the DDofs for **unknown** 1, which should be component 1. This is a bug in the gsElasticity assembler
        // fixedDofs = elAssembler.allFixedDofs();
        // elAssembler.setFixedDofs(fixedDofs);

        // sigma_left.set_u(tcurr/tend * sigma);
        // sigma_right.set_u(tcurr/tend * sigma);

        gsInfo<<"---------------------------------------------------------------------------------------------------------------------------\n";
        gsInfo<<"Load step "<<step<<": t = "<<tcurr<<"\n\n";

        delta_u.setZero();
        delta_D.setZero();

        // Prediction step (IGA book Eqs. (6.44)-(6.46))
        udot_new = udot_old;
        uddot_new = (gamma-1)/gamma * uddot_old;
        u_new = u_old + dt * udot_old + 0.5*math::pow(dt,2) * ((1-2*beta) * uddot_old + 2*beta * uddot_new);
        
        Unorm = u_new.norm();
        U0 = (Unorm > 0) ? Unorm : 1.0; // Avoid division by zero
        Rnorm = R0 = 1;

        D_new = D_old;

        gsInfo<< "INITIAL VALUES"
                <<": ||U|| = "<<Unorm
                <<", ||U||/||U0|| = "<<Unorm/U0
                <<"\n";

        gsInfo<<"\t"<<PRINT(20)<<"* Elasticity:"<<PRINT(6)<<"It."<<PRINT(18)<<"||R||"<<PRINT(18)<<"||R||/||R0||"<<PRINT(18)<<"||ΔU||/||U0||"<<PRINT(18)<<"||dA||/||A||"<<PRINT(18)<<"||U||"<<PRINT(18)<<"||V||"<<PRINT(18)<<"||A||"<<PRINT(20)<<PRINT(20)<<"cum. assembly [s]"<<PRINT(20)<<"cum. solver [s]"<<"\n";


        for (index_t it=0; it!=maxIt; ++it)
        {
            
            elAssemblyTime = elSolverTime = 0.0;
            pfAssemblyTime = pfSolverTime = 0.0;
            iterationTime  = 0.0;
            bigClock.restart();
            // gsInfo<<" - Staggered iteration "<<it<<":\n";

            // gsInfo<<"  D_new.maxCoeff()" << D_new.maxCoeff() << "\n";
            // gsInfo<<"  D_new.minCoeff()" << D_new.minCoeff() << "\n";

            material.setParameter(2,damage);
            elAssembler.initialize();

            // ================================================ ELASTICITY ==============================================
            smallClock.restart();
            elAssembler.assemble(u_new);
            elAssemblyTime += smallClock.stop();
            elAssembler.matrix_into(K);
            elAssembler.rhs_into(Fext);

            // // --- STATIC TEST (diagnostic) ---
            // elAssembler.assemble(u_new);
            // elAssembler.matrix_into(K);
            // elAssembler.rhs_into(Fext);

            // Try solve K * u_static = Fext (after BCs have been applied by assembler)
            solver.compute(K);
            if (solver.info() != 0) gsInfo<<"DEBUG static solver.compute failed: "<<solver.info()<<"\n";
            gsMatrix<T> u_static = solver.solve(Fext);
            // gsInfo << "DEBUG STATIC: ||K*u_static - Fext|| = " << (K*u_static - Fext).norm()
            //     << " ||u_static|| = " << u_static.norm() << "\n";

            R = M * uddot_new + K * u_new - Fext;
            Rnorm = R.norm();
            if (it == 0)
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

            // gsInfo << "MAX DELTA U = " << delta_u.maxCoeff() << "\n";


            gsInfo<<"\t"<<PRINT(20)<<""<<PRINT(6)<<it<<PRINT(18)<<Rnorm<<PRINT(18)<<Rnorm/R0<<PRINT(18)<<DeltaUnorm/U0<<PRINT(18)<<(uddot_new-uddot_old).norm()/uddot_new.norm()<<PRINT(18)<<Unorm<<PRINT(18)<<udot_new.norm()<<PRINT(18)<<uddot_new.norm()<<PRINT(20)<<elAssemblyTime<<PRINT(20)<<elSolverTime<<"\n";

            // Update the fields before breaking the staggered loop 
            // (otherwise it does not update the displacements if staggered converges in 1 iteration)
            elAssembler.constructSolution(u_new,displacement);

            for (size_t p=0; p!=mp.nPatches(); ++p)
                mp_def.patch(p).coefs() = mp.patch(p).coefs() + displacement.patch(p).coefs();


            // gsInfo<<"  ----------------------------------------\n";
            // gsInfo<< "maximum displacement: "<<displacement.patch(0).coefs().maxCoeff() << "\n";
            // gsInfo<< "minimum displacement: "<<displacement.patch(0).coefs().minCoeff() << "\n";


            // Initialize the function for the elastic energy
            gsMaterialEval<T,gsMaterialOutput::Psi> Psi(&material,mp,mp_def);

            // gsInfo << "DEBUG: u_new.maxCoeff() = " << u_new.maxCoeff()
            //     << ", displacement.patch(0).coefs().maxCoeff() = "
            //     << displacement.patch(0).coefs().maxCoeff() << "\n";

            // Staggered tolerance check for elasticity (according to Remark 1 in Greco et al. 2025)
            // if (Rnorm/R0 < 1e-5 && Unorm/U0 < 1e-4) // does not work!
            // gsInfo<< "  DeltaUnorm/U0: " << DeltaUnorm/U0 << "\n";
            // gsInfo<< "  Rnorm/R0 : " << Rnorm/R0  << "\n";
            // gsInfo<< "  DeltaUnorm: " << DeltaUnorm << "\n";
            // gsInfo<< "  Rnorm: " << Rnorm << "\n";

            if (Rnorm/R0 < 1e-5 && DeltaUnorm/U0 < 1e-4)
                break;
            else if (it == maxIt-1)
                GISMO_ERROR("Staggered iterations problem did not converge.");
        }

        // =========================================================================
        // Compute resulting force and energies
        gsMatrix<T> ufull = displacement.patch(0).coefs().reshape(displacement.patch(0).coefs().size(),1);
        fullElAssembler.assemble(ufull);
        gsMatrix<T> Rfull = fullElAssembler.rhs();
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
        stepData[0] = tcurr;
        stepData[1] = Fx;
        stepData[2] = Fy;
        stepData[3] = (0.5 * ufull.transpose() * fullElAssembler.matrix() * ufull).value();
        stepData[4] = (0.5 * D_new.transpose() * QPhi * D_new).value() + (D_new.transpose() * q).value();

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
            np[0] = 1000;
            if (dim==2)
                np[1] = 2;
            else
                np[1] = np[2] = 2;
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
        }

        // =========================================================================
        // Write data
        std::ofstream file(outputdir+"results.txt",std::ios::app);
        // for (size_t i = 0; i != data.size(); ++i)
        //     file<<data[i][0]<<","<<-data[i][1]<<","<<-data[i][2]<<","<<data[i][3]<<","<<data[i][4]<<"\n";
        file<<stepData[0]<<","<<-stepData[1]<<","<<-stepData[2]<<","<<stepData[3]<<","<<stepData[4]<<"\n";
        file.close();

        u_old = u_new;
        udot_old = udot_new;
        uddot_old = uddot_new;
        D_old = D_new;

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
    }

    delete pfAssembler;
}