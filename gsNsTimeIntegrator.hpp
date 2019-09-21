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
#include <gsElasticity/gsNewton.h>

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
    opt.addInt("Scheme","Time integration scheme",time_integration_NS::theta_scheme_linear);
    opt.addReal("Theta","Time integration parametrer: 0 - explicit Euler, 1 - implicit Euler, 0.5 - Crank-Nicolson",0.5);
    opt.addInt("Verbosity","Amount of information printed to the terminal: none, some, all",newton_verbosity::none);
    return opt;
}

template <class T>
void gsNsTimeIntegrator<T>::initialize()
{
    if (solVector.rows() != stiffAssembler.numDofs())
        solVector.setZero(stiffAssembler.numDofs(),1);

    if (m_options.getInt("Scheme") == time_integration_NS::theta_scheme)
        stiffAssembler.options().setInt("Iteration",iteration_type::newton);
    if (m_options.getInt("Scheme") == time_integration_NS::theta_scheme_linear)
        stiffAssembler.options().setInt("Iteration",iteration_type::picard);

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
    if (m_options.getInt("Scheme") == time_integration_NS::theta_scheme)
        implicitNonlinear();
    if (m_options.getInt("Scheme") == time_integration_NS::theta_scheme_linear)
        implicitLinear();
}

template <class T>
void gsNsTimeIntegrator<T>::implicitLinear()
{
    T theta = m_options.getReal("Theta");
    index_t numDofsVel = massAssembler.numDofs();

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
    m_ddof = stiffAssembler.allFixedDofs();

    gsSparseSolver<>::LU solver(m_system.matrix());
    solVector = solver.solve(m_system.rhs());
}

template <class T>
void gsNsTimeIntegrator<T>::implicitNonlinear()
{
    /*stiffAssembler.assemble(solVector,m_ddof,false);
    oldResidual = stiffAssembler.rhs();

    gsNewton<T> solver(*this,solVector,m_ddof);
    solver.options().setInt("Verbosity",m_options.getInt("Verbosity"));
    solver.options().setInt("Solver",linear_solver::LU);
    solver.solve();
    solVector = solver.solution();
    m_ddof = solver.allFixedDoFs();*/
}

template <class T>
bool gsNsTimeIntegrator<T>::assemble(const gsMatrix<T> & solutionVector,
                                     const std::vector<gsMatrix<T> > & fixedDoFs,
                                     bool assembleMatrix)
{
    T theta = m_options.getReal("Theta");
    stiffAssembler.options().setInt("Iteration",iteration_type::newton);
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
    m_system.rhs().middleRows(0,numDofsVel).noalias() += massAssembler.matrix()*((solVector-solutionVector).middleRows(0,numDofsVel));

    return true;
}



} // namespace ends
