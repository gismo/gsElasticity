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
}

template <class T>
gsOptionList gsNsTimeIntegrator<T>::defaultOptions()
{
    gsOptionList opt = Base::defaultOptions();
    opt.addInt("Scheme","Time integration scheme",time_integration::implicit_nonlinear);
    opt.addInt("Verbosity","Amount of information printed to the terminal: none, some, all",newton_verbosity::none);
    return opt;
}

template <class T>
void gsNsTimeIntegrator<T>::initialize()
{
    massAssembler.assemble();
    //stiffAssembler.assemble();

    if (solVector.rows() != stiffAssembler.numDofs())
        solVector.setZero(stiffAssembler.numDofs(),1);
}

template <class T>
void gsNsTimeIntegrator<T>::makeTimeStep(T timeStep)
{
    tStep = timeStep;
    if (m_options.getInt("Scheme") == time_integration_NS::implicit_euler_newton)
        solVector = implicitNewton();
    if (m_options.getInt("Scheme") == time_integration_NS::implicit_euler_oseen)
        solVector = implicitOseen();
    if (m_options.getInt("Scheme") == time_integration_NS::crank_nicolson_oseen)
        solVector = semiImplicitOseen();
}

template <class T>
gsMatrix<T> gsNsTimeIntegrator<T>::implicitNewton()
{
    gsNewton<T> solver(*this,solVector);
    solver.options().setInt("Verbosity",m_options.getInt("Verbosity"));
    solver.options().setInt("MaxIters",1);
    solver.solve();
    return solver.solution();
}

template <class T>
bool gsNsTimeIntegrator<T>::assemble(const gsMatrix<T> & solutionVector, bool assembleMatrix)
{
    if (!stiffAssembler.assemble(solutionVector))
        return false;
    Base::m_system.matrix() = tStep*stiffAssembler.matrix();
    index_t numDofsVel = massAssembler.numDofs();
    gsSparseMatrix<T> tempMassMatrix = massAssembler.matrix();
    tempMassMatrix.conservativeResize(stiffAssembler.numDofs(),numDofsVel);
    Base::m_system.matrix().leftCols(numDofsVel) += tempMassMatrix;
    Base::m_system.rhs() = tStep*stiffAssembler.rhs();
    Base::m_system.rhs().middleRows(0,numDofsVel).noalias() += massAssembler.matrix()*(solVector-solutionVector);

    return true;
}


template <class T>
gsMatrix<T> gsNsTimeIntegrator<T>::implicitOseen()
{
    stiffAssembler.options().setInt("Iteration",iteration_type::picard);
    stiffAssembler.assemble(solVector);
    m_system.matrix() = tStep*stiffAssembler.matrix();
    index_t numDofsVel = massAssembler.numDofs();
    gsSparseMatrix<T> tempMassMatrix = massAssembler.matrix();
    tempMassMatrix.conservativeResize(stiffAssembler.numDofs(),numDofsVel);
    m_system.matrix().leftCols(numDofsVel) += tempMassMatrix;
    m_system.matrix().makeCompressed();
    m_system.rhs() = tStep*stiffAssembler.rhs();
    m_system.rhs().middleRows(0,numDofsVel).noalias() += massAssembler.matrix()*solVector.middleRows(0,numDofsVel);// + massAssembler.rhs();

    gsSparseSolver<>::LU solver(m_system.matrix());
    return solver.solve(m_system.rhs());
}

template <class T>
gsMatrix<T> gsNsTimeIntegrator<T>::semiImplicitOseen()
{
    stiffAssembler.options().setInt("Iteration",iteration_type::picard);
    stiffAssembler.assemble(solVector);
    m_system.matrix() = 0.5*tStep*stiffAssembler.matrix();
    index_t numDofsVel = massAssembler.numDofs();
    gsSparseMatrix<T> tempMassMatrix = massAssembler.matrix();
    tempMassMatrix.conservativeResize(stiffAssembler.numDofs(),numDofsVel);
    m_system.matrix().leftCols(numDofsVel) += tempMassMatrix;
    m_system.matrix().makeCompressed();

    m_system.rhs() = tStep*(stiffAssembler.rhs()-0.5*stiffAssembler.matrix()*solVector);
    m_system.rhs().middleRows(0,numDofsVel).noalias() += massAssembler.matrix()*solVector.middleRows(0,numDofsVel);

    gsSparseSolver<>::LU solver(m_system.matrix());
    return solver.solve(m_system.rhs());
}

} // namespace ends
