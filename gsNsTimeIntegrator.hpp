/** @file gsNsTimeIntegrator.hpp

    @brief A class providing time integration for incompressible Navier-Stokes.

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s):
        O. Weeger    (2012 - 2015, TU Kaiserslautern),
        A.Shamanskiy (2016 - ...., TU Kaiserslautern)
*/

#pragma once

#include <gsElasticity/gsNsTimeIntegrator.h>

#include <gsElasticity/gsNsAssembler.h>
#include <gsElasticity/gsMassAssembler.h>
#include <gsElasticity/gsIterative.h>

namespace gismo
{

template <class T>
gsNsTimeIntegrator<T>::gsNsTimeIntegrator(gsNsAssembler<T> & stiffAssembler_,
                                          gsMassAssembler<T> & massAssembler_)
    : stiffAssembler(stiffAssembler_),
      massAssembler(massAssembler_)
{
    m_options = defaultOptions();
    m_ddof = stiffAssembler.allFixedDofs();
}

template <class T>
gsOptionList gsNsTimeIntegrator<T>::defaultOptions()
{
    gsOptionList opt = Base::defaultOptions();
    opt.addInt("Scheme","Time integration scheme",time_integration::implicit_nonlinear);
    opt.addReal("Theta","Time integration parametrer: 0 - explicit Euler, 1 - implicit Euler, 0.5 - Crank-Nicolson",0.5);
    opt.addInt("Verbosity","Amount of information printed to the terminal: none, some, all",solver_verbosity::none);
    opt.addReal("AbsTol","Absolute tolerance for the convergence cretiria",1e-10);
    opt.addReal("RelTol","Relative tolerance for the stopping criteria",1e-7);
    return opt;
}

template <class T>
void gsNsTimeIntegrator<T>::initialize()
{
    if (solVector.rows() != stiffAssembler.numDofs())
        solVector.setZero(stiffAssembler.numDofs(),1);

    if (m_options.getInt("Scheme") == time_integration::implicit_nonlinear)
        stiffAssembler.options().setInt("Assembly",ns_assembly::newton_next);
    if (m_options.getInt("Scheme") == time_integration::implicit_linear)
        stiffAssembler.options().setInt("Assembly",ns_assembly::ossen);

    stiffAssembler.assemble(solVector,m_ddof);

    oldSolVector = solVector;
    massAssembler.setFixedDofs(m_ddof);
    massAssembler.assemble();
    oldTimeStep = 1.;
}

template <class T>
void gsNsTimeIntegrator<T>::makeTimeStep(T timeStep)
{
    tStep = timeStep;
    if (m_options.getInt("Scheme") == time_integration::implicit_nonlinear)
        implicitNonlinear();
    if (m_options.getInt("Scheme") == time_integration::implicit_linear)
        implicitLinear();
}

template <class T>
void gsNsTimeIntegrator<T>::implicitLinear()
{
    T theta = m_options.getReal("Theta");
    index_t numDofsVel = massAssembler.numDofs();
    stiffAssembler.options().setInt("Assembly",ns_assembly::ossen);

    // rhs = M*u_n - dt*(1-theta)*A(u_n)*u_n + dt*(1-theta)*F_n + dt*theta*F_n+1
    // rhs: dt*(1-theta)*F_n
    m_system.rhs() = tStep*(1-theta)*stiffAssembler.rhs();
    // rhs: -dt*(1-theta)*A(u_n)*u_n
    m_system.rhs().middleRows(0,numDofsVel) -= tStep*(1-theta)*stiffAssembler.matrix().block(0,0,numDofsVel,numDofsVel) *
                                                               solVector.middleRows(0,numDofsVel);
    // rhs: M*u_n
    m_system.rhs().middleRows(0,numDofsVel) += massAssembler.matrix() *
                                               solVector.middleRows(0,numDofsVel);
    // rhs: M_FD*u_DDOFS_n
    m_system.rhs().middleRows(0,numDofsVel) -= massAssembler.rhs();

    stiffAssembler.assemble(solVector + tStep/oldTimeStep*(solVector-oldSolVector),
                            stiffAssembler.allFixedDofs());
    massAssembler.setFixedDofs(stiffAssembler.allFixedDofs());
    massAssembler.assemble(); // need to redo this - only need the rhs with elimated DDOFS
    // rhs: dt*theta*F_n+1
    m_system.rhs() += tStep*theta*stiffAssembler.rhs();
    // rhs: -M_FD*u_DDOFS_n+1
    m_system.rhs().middleRows(0,numDofsVel) += massAssembler.rhs();
    // matrix = M + dt*theta*A(u_exp)
    m_system.matrix() = tStep*stiffAssembler.matrix();
    // we need to modify the (0,0,numDofsVel,numDofsVel) block of the matrix.
    // unfortunately, eigen provides only a read-only block interfaces.
    // the following is an ugly way to overcome this
    gsSparseMatrix<T> tempVelocityBlock = m_system.matrix().block(0,0,numDofsVel,numDofsVel);
    tempVelocityBlock *= (theta-1);
    tempVelocityBlock += massAssembler.matrix();
    tempVelocityBlock.conservativeResize(stiffAssembler.numDofs(),numDofsVel);
    m_system.matrix().leftCols(numDofsVel) += tempVelocityBlock;

    m_system.matrix().makeCompressed();
    oldSolVector = solVector;
    oldTimeStep = tStep;
    m_ddof = stiffAssembler.allFixedDofs();

#ifdef GISMO_WITH_PARDISO
    gsSparseSolver<>::PardisoLU solver(m_system.matrix());
    solVector = solver.solve(m_system.rhs());
#else
    gsSparseSolver<>::LU solver(m_system.matrix());
    solVector = solver.solve(m_system.rhs());
#endif


    numIters = 1;
}

template <class T>
void gsNsTimeIntegrator<T>::makeTimeStepFSI2(T timeStep,gsMultiPatch<T> & velocityALE,
                      std::vector<std::pair<index_t,index_t> > & patches)
{
    T theta = m_options.getReal("Theta");
    index_t numDofsVel = massAssembler.numDofs();
    stiffAssembler.options().setInt("Assembly",ns_assembly::ossen);

    // rhs = M*u_n - dt*(1-theta)*A(u_n)*u_n + dt*(1-theta)*F_n + dt*theta*F_n+1
    // rhs: dt*(1-theta)*F_n
    m_system.rhs() = timeStep*(1-theta)*stiffAssembler.rhs();
    // rhs: -dt*(1-theta)*A(u_n)*u_n
    m_system.rhs().middleRows(0,numDofsVel) -= timeStep*(1-theta)*stiffAssembler.matrix().block(0,0,numDofsVel,numDofsVel) *
                                                                  solVector.middleRows(0,numDofsVel);
    // rhs: M*u_n
    m_system.rhs().middleRows(0,numDofsVel) += massAssembler.matrix() *
                                               solVector.middleRows(0,numDofsVel);
    // rhs: M_FD*u_DDOFS_n
    m_system.rhs().middleRows(0,numDofsVel) -= massAssembler.rhs();

    gsMultiPatch<T> velocity, pressure;
    stiffAssembler.constructSolution(solVector + timeStep/oldTimeStep*(solVector-oldSolVector),
                                     stiffAssembler.allFixedDofs(),velocity,pressure);
    for (auto &it : patches)
        velocity.patch(it.first).coefs() -= velocityALE.patch(it.second).coefs();
    stiffAssembler.assemble(velocity,pressure);
    massAssembler.setFixedDofs(stiffAssembler.allFixedDofs());
    massAssembler.assemble(); // need to redo this - only need the rhs with elimated DDOFS

    // rhs: dt*theta*F_n+1
    m_system.rhs() += timeStep*theta*stiffAssembler.rhs();
    // rhs: -M_FD*u_DDOFS_n+1
    m_system.rhs().middleRows(0,numDofsVel) += massAssembler.rhs();
    // matrix = M + dt*theta*A(u_exp)
    m_system.matrix() = timeStep*stiffAssembler.matrix();
    // we need to modify the (0,0,numDofsVel,numDofsVel) block of the matrix.
    // unfortunately, eigen provides only a read-only block interfaces.
    // the following is an ugly way to overcome this
    gsSparseMatrix<T> tempVelocityBlock = m_system.matrix().block(0,0,numDofsVel,numDofsVel);
    tempVelocityBlock *= (theta-1);
    tempVelocityBlock += massAssembler.matrix();
    tempVelocityBlock.conservativeResize(stiffAssembler.numDofs(),numDofsVel);
    m_system.matrix().leftCols(numDofsVel) += tempVelocityBlock;

    m_system.matrix().makeCompressed();
    oldSolVector = solVector;
    oldTimeStep = timeStep;
    m_ddof = stiffAssembler.allFixedDofs();

#ifdef GISMO_WITH_PARDISO
    gsSparseSolver<>::PardisoLU solver(m_system.matrix());
    solVector = solver.solve(m_system.rhs());
#else
    gsSparseSolver<>::LU solver(m_system.matrix());
    solVector = solver.solve(m_system.rhs());
#endif


    numIters = 1;
}

template <class T>
void gsNsTimeIntegrator<T>::implicitNonlinear()
{
    stiffAssembler.options().setInt("Assembly",ns_assembly::newton_next);
    T theta = m_options.getReal("Theta");
    index_t numDofsVel = massAssembler.numDofs();

    ALE = false;

    constRHS = tStep*(1-theta)*stiffAssembler.rhs();
    constRHS.middleRows(0,numDofsVel).noalias() -= tStep*(1-theta)*stiffAssembler.matrix().block(0,0,numDofsVel,numDofsVel)*solVector.middleRows(0,numDofsVel);
    constRHS.middleRows(0,numDofsVel).noalias() += massAssembler.matrix()*solVector.middleRows(0,numDofsVel);
    constRHS.middleRows(0,numDofsVel).noalias() -= massAssembler.rhs();
    massAssembler.setFixedDofs(stiffAssembler.allFixedDofs());
    massAssembler.assemble();
    constRHS.middleRows(0,numDofsVel).noalias() += massAssembler.rhs();

    gsIterative<T> solver(*this,solVector,m_ddof);
    solver.options().setInt("Verbosity",m_options.getInt("Verbosity"));
    solver.options().setInt("Solver",linear_solver::LU);
    solver.options().setInt("IterType",iteration_type::next);
    solver.options().setReal("AbsTol",m_options.getReal("AbsTol"));
    solver.options().setReal("RelTol",m_options.getReal("RelTol"));
    solver.solve();

    numIters = solver.numberIterations();

    solVector = solver.solution();
    m_ddof = stiffAssembler.allFixedDofs();
}

template <class T>
void gsNsTimeIntegrator<T>::makeTimeStepFSI3(T timeStep,gsMultiPatch<T> & velocityALE,
                                             std::vector<std::pair<index_t,index_t> > & patches)
{
    tStep = timeStep;
    stiffAssembler.options().setInt("Assembly",ns_assembly::newton_next);
    T theta = m_options.getReal("Theta");
    index_t numDofsVel = massAssembler.numDofs();

    ALE = true;
    aleVelocity.clear();
    for (index_t p = 0; p < velocityALE.nPatches();++p)
        aleVelocity.addPatch(velocityALE.patch(p).clone());
    alePatches.clear();
    for (auto it : patches)
        alePatches.push_back(it);

    constRHS = tStep*(1-theta)*stiffAssembler.rhs();
    constRHS.middleRows(0,numDofsVel).noalias() -= tStep*(1-theta)*stiffAssembler.matrix().block(0,0,numDofsVel,numDofsVel)*solVector.middleRows(0,numDofsVel);
    constRHS.middleRows(0,numDofsVel).noalias() += massAssembler.matrix()*solVector.middleRows(0,numDofsVel);
    constRHS.middleRows(0,numDofsVel).noalias() -= massAssembler.rhs();
    massAssembler.setFixedDofs(stiffAssembler.allFixedDofs());
    massAssembler.assemble();
    constRHS.middleRows(0,numDofsVel).noalias() += massAssembler.rhs();

    gsIterative<T> solver(*this,solVector,m_ddof);
    solver.options().setInt("Verbosity",m_options.getInt("Verbosity"));
    solver.options().setInt("Solver",linear_solver::LU);
    solver.options().setInt("IterType",iteration_type::next);
    solver.options().setReal("AbsTol",m_options.getReal("AbsTol"));
    solver.options().setReal("RelTol",m_options.getReal("RelTol"));
    solver.solve();


    numIters = solver.numberIterations();

    solVector = solver.solution();
    m_ddof = stiffAssembler.allFixedDofs();
}

template <class T>
bool gsNsTimeIntegrator<T>::assemble(const gsMatrix<T> & solutionVector,
                                     const std::vector<gsMatrix<T> > & fixedDoFs,
                                     bool assembleMatrix)
{
    T theta = m_options.getReal("Theta");
    if (!ALE)
        stiffAssembler.assemble(solutionVector,fixedDoFs,assembleMatrix);
    else
    {
        gsMultiPatch<T> velocity, pressure;
        stiffAssembler.constructSolution(solutionVector,fixedDoFs,velocity,pressure);
        for (auto &it : alePatches)
            velocity.patch(it.first).coefs() -= aleVelocity.patch(it.second).coefs();
        stiffAssembler.assemble(velocity,pressure);
    }

    m_system.matrix() = tStep*stiffAssembler.matrix();
    index_t numDofsVel = massAssembler.numDofs();
    gsSparseMatrix<T> tempVelocityBlock = m_system.matrix().block(0,0,numDofsVel,numDofsVel);
    tempVelocityBlock *= (theta-1);
    tempVelocityBlock += massAssembler.matrix();
    tempVelocityBlock.conservativeResize(stiffAssembler.numDofs(),numDofsVel);
    m_system.matrix().leftCols(numDofsVel) += tempVelocityBlock;
    m_system.matrix().makeCompressed();

    m_system.rhs() = tStep*theta*stiffAssembler.rhs() + constRHS;
    return true;
}

template <class T>
void gsNsTimeIntegrator<T>::saveState()
{
    velVecSaved = solVector;
    oldVecSaved = oldSolVector;
    massRhsSaved = massAssembler.rhs();
    stiffRhsSaved = stiffAssembler.rhs();
    stiffMatrixSaved = stiffAssembler.matrix();
    ddofsSaved = m_ddof;
}

template <class T>
void gsNsTimeIntegrator<T>::recoverState()
{
    GISMO_ASSERT(velVecSaved.rows() == solVector.rows(),"No state saved!");
    solVector = velVecSaved;
    oldSolVector = oldVecSaved;
    massAssembler.setRHS(massRhsSaved);
    stiffAssembler.setMatrix(stiffMatrixSaved);
    stiffAssembler.setRHS(stiffRhsSaved);
    m_ddof = ddofsSaved;
}



} // namespace ends
