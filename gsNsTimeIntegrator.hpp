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
    if (m_options.getInt("Scheme") == time_integration::explicit_lumped)
        control();
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
    stiffAssembler.assemble(solVector + tStep/oldTimeStep*(solVector-oldSolVector),
                            stiffAssembler.allFixedDofs());
    // rhs: dt*theta*F_n+1
    m_system.rhs() += tStep*theta*stiffAssembler.rhs();

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
}

template <class T>
void gsNsTimeIntegrator<T>::control()
{
    index_t numDofsVel = massAssembler.numDofs();
    stiffAssembler.options().setInt("Assembly",ns_assembly::ossen);
    massAssembler.assemble(true);

    // rhs = M*u_n - dt*(1-theta)*A(u_n)*u_n + dt*(1-theta)*F_n + dt*theta*F_n+1
    // rhs: dt*(1-theta)*F_n
    //m_system.rhs() = tStep*(1-theta)*stiffAssembler.rhs();
    // rhs: -dt*(1-theta)*A(u_n)*u_n
    //m_system.rhs().middleRows(0,numDofsVel) -= tStep*(1-theta)*stiffAssembler.matrix().block(0,0,numDofsVel,numDofsVel) *
    //                                                           solVector.middleRows(0,numDofsVel);
    // rhs: M*u_n
    //m_system.rhs().middleRows(0,numDofsVel) += massAssembler.matrix() *
    //                                           solVector.middleRows(0,numDofsVel);

    for (index_t d = 0; d < m_ddof.size(); ++d)
        m_ddof[d] += massAssembler.allFixedDofs()[d];
    stiffAssembler.assemble(solVector + tStep/oldTimeStep*(solVector-oldSolVector),
                            m_ddof);
    // rhs: dt*theta*F_n+1
    m_system.rhs() = tStep*stiffAssembler.rhs();
    gsInfo << "Arhs " << stiffAssembler.rhs().norm() << std::endl;

    gsInfo << "Mrhs " << massAssembler.rhs().norm() << std::endl;
    m_system.rhs().middleRows(0,numDofsVel) += massAssembler.matrix() *
                                                   solVector.middleRows(0,numDofsVel);
    m_system.rhs().middleRows(0,numDofsVel) += massAssembler.rhs();

    // matrix = M + dt*theta*A(u_exp)
    m_system.matrix() = tStep*stiffAssembler.matrix();
    // we need to modify the (0,0,numDofsVel,numDofsVel) block of the matrix.
    // unfortunately, eigen provides only a read-only block interfaces.
    // the following is an ugly way to overcome this
    gsSparseMatrix<T> tempVelocityBlock = m_system.matrix().block(0,0,numDofsVel,numDofsVel);
    tempVelocityBlock *= 0.;
    tempVelocityBlock += massAssembler.matrix();
    tempVelocityBlock.conservativeResize(stiffAssembler.numDofs(),numDofsVel);
    m_system.matrix().leftCols(numDofsVel) += tempVelocityBlock;

    m_system.matrix().makeCompressed();
    oldSolVector = solVector;
    oldTimeStep = tStep;

#ifdef GISMO_WITH_PARDISO
    gsSparseSolver<>::PardisoLU solver(m_system.matrix());
    solVector = solver.solve(m_system.rhs());
#else
    gsSparseSolver<>::LU solver(m_system.matrix());
    solVector = solver.solve(m_system.rhs());
#endif
}

template <class T>
void gsNsTimeIntegrator<T>::makeTimeStepFSI(T timeStep,
                                            gsMatrix<T> & solutionVector, gsMatrix<T> & solutionVectorOld,
                                            gsMultiPatch<T> & velocityALE,
                                            std::vector<std::pair<index_t,index_t> > & patches,
                                            gsSparseMatrix<T> & A_n, gsMatrix<T> & rhs_n)
{
    T theta = m_options.getReal("Theta");
    index_t numDofsVel = massAssembler.numDofs();
    stiffAssembler.options().setInt("Assembly",ns_assembly::ossen);

    // rhs = M*u_n - dt*(1-theta)*A(u_n-u_ALE)*u_n + dt*(1-theta)*F_n + dt*theta*F_n+1
    // rhs: dt*(1-theta)*F_n
    m_system.rhs() = timeStep*(1-theta)*rhs_n;
    // rhs: -dt*(1-theta)*A(u_n)*u_n
    m_system.rhs().middleRows(0,numDofsVel) -= timeStep*(1-theta)*A_n.block(0,0,numDofsVel,numDofsVel) *
                                                                solutionVector.middleRows(0,numDofsVel);
    // rhs: M*u_n
    m_system.rhs().middleRows(0,numDofsVel) += massAssembler.matrix() *
                                               solutionVector.middleRows(0,numDofsVel);
    gsMultiPatch<T> velocity, pressure;
    stiffAssembler.constructSolution(2*solutionVector -solutionVectorOld,
                                     stiffAssembler.allFixedDofs(),velocity,pressure);
    for (auto it : patches)
        velocity.patch(it.first).coefs() -= velocityALE.patch(it.second).coefs();
    //gsInfo << velocityALE.patch(0).coefs().norm() << " "
    //       <<velocityALE.patch(1).coefs().norm() << " "
    //      <<velocityALE.patch(2).coefs().norm() << std::endl;
    //velocity.patch(3).coefs() -= velocityALE.patch(0).coefs();
    stiffAssembler.assemble(velocity,pressure);
    // rhs: dt*theta*F_n+1
    m_system.rhs() += timeStep*theta*stiffAssembler.rhs();

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
    m_ddof = stiffAssembler.allFixedDofs();

#ifdef GISMO_WITH_PARDISO
    gsSparseSolver<>::PardisoLU solver(m_system.matrix());
    solVector = solver.solve(m_system.rhs());
#else
    gsSparseSolver<>::LU solver(m_system.matrix());
    solVector = solver.solve(m_system.rhs());
#endif
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

    gsMultiPatch<T> velocity, pressure;
    stiffAssembler.constructSolution(solVector,// + timeStep/oldTimeStep*(solVector-oldSolVector),
                                     stiffAssembler.allFixedDofs(),velocity,pressure);
    for (auto &it : patches)
        velocity.patch(it.first).coefs() -= velocityALE.patch(it.second).coefs();
    stiffAssembler.assemble(velocity,pressure);
    // rhs: dt*theta*F_n+1
    m_system.rhs() += timeStep*theta*stiffAssembler.rhs();

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
}

template <class T>
void gsNsTimeIntegrator<T>::implicitNonlinear()
{
    stiffAssembler.options().setInt("Assembly",ns_assembly::ossen);

    std::vector<gsMatrix<T> > tempDDofs = stiffAssembler.allFixedDofs();
    stiffAssembler.setFixedDofs(m_ddof);
    stiffAssembler.assemble(solVector,m_ddof);
    oldResidual = stiffAssembler.rhs();

    index_t numDofsVel = massAssembler.numDofs();
    oldResidual.middleRows(0,numDofsVel).noalias() -= stiffAssembler.matrix().block(0,0,numDofsVel,numDofsVel)*solVector.middleRows(0,numDofsVel);

    stiffAssembler.setFixedDofs(tempDDofs);
    gsIterative<T> solver(*this,solVector,m_ddof);
    solver.options().setInt("Verbosity",m_options.getInt("Verbosity"));
    solver.options().setInt("Solver",linear_solver::LU);
    solver.options().setInt("IterType",iteration_type::next);
    solver.options().setReal("AbsTol",m_options.getReal("AbsTol"));
    solver.options().setReal("RelTol",m_options.getReal("RelTol"));
    solver.solve();

    solVector = solver.solution();
    m_ddof = stiffAssembler.allFixedDofs();
}

template <class T>
bool gsNsTimeIntegrator<T>::assemble(const gsMatrix<T> & solutionVector,
                                     const std::vector<gsMatrix<T> > & fixedDoFs,
                                     bool assembleMatrix)
{
    T theta = m_options.getReal("Theta");
    stiffAssembler.options().setInt("Assembly",ns_assembly::newton_next);
    stiffAssembler.assemble(solutionVector,fixedDoFs,assembleMatrix);

    m_system.matrix() = tStep*stiffAssembler.matrix();
    index_t numDofsVel = massAssembler.numDofs();
    gsSparseMatrix<T> tempVelocityBlock = m_system.matrix().block(0,0,numDofsVel,numDofsVel);
    tempVelocityBlock *= (theta-1);
    tempVelocityBlock += massAssembler.matrix();
    tempVelocityBlock.conservativeResize(stiffAssembler.numDofs(),numDofsVel);
    m_system.matrix().leftCols(numDofsVel) += tempVelocityBlock;
    m_system.matrix().makeCompressed();

    m_system.rhs() = tStep*(stiffAssembler.rhs()*theta + (1-theta)*oldResidual);
    m_system.rhs().middleRows(0,numDofsVel).noalias() += massAssembler.matrix()*solVector.middleRows(0,numDofsVel);

    return true;
}



} // namespace ends
