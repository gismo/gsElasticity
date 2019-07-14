/** @file gsETimeIntegrator.hpp

    @brief A class providing time integration for dynamical elasticity.

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s):
        O. Weeger    (2012 - 2015, TU Kaiserslautern),
        A.Shamanskiy (2016 - ...., TU Kaiserslautern)
*/

#pragma once

#include <gsElasticity/gsElTimeIntegrator.h>

#include <gsElasticity/gsElasticityAssembler.h>
#include <gsElasticity/gsElMassAssembler.h>
#include <gsElasticity/gsElNewton.h>


namespace gismo
{

template <class T>
gsElTimeIntegrator<T>::gsElTimeIntegrator(gsElasticityAssembler<T> & stiffAssembler_,
                                          gsElMassAssembler<T> & massAssembler_)
    : stiffAssembler(stiffAssembler_),
      massAssembler(massAssembler_),
      timeStepNum(0)
{
    massAssembler.assemble();
    stiffAssembler.assemble();
    dispVector.setZero(stiffAssembler.numDofs(),1);
    velVector.setZero(stiffAssembler.numDofs(),1);
    gsSparseSolver<>::SimplicialLDLT solver(massAssembler.matrix());
    accVector = solver.solve(-1*stiffAssembler.matrix()*dispVector+stiffAssembler.rhs());

    m_options = defaultOptions();
}

template <class T>
gsOptionList gsElTimeIntegrator<T>::defaultOptions()
{
    gsOptionList opt = Base::defaultOptions();
    return opt;
}

template <class T>
void gsElTimeIntegrator<T>::makeTimeStep(T timeStep)
{
    T beta = 0.25;
    T gamma = 0.5;

    T alpha1 = 1./beta/pow(timeStep,2);
    T alpha2 = 1./beta/timeStep;
    T alpha3 = (1-2*beta)/2/beta;
    T alpha4 = gamma/beta/timeStep;
    T alpha5 = 1 - gamma/beta;
    T alpha6 = (1 - gamma/beta/2)*timeStep;

    m_matrix = alpha1*massAssembler.matrix() + stiffAssembler.matrix();
    m_rhs = massAssembler.matrix()*(alpha1*dispVector + alpha2*velVector + alpha3*accVector) + stiffAssembler.rhs();
    gsSparseSolver<>::SimplicialLDLT solver(m_matrix);
    gsMatrix<T> newDispVector = solver.solve(m_rhs);
    gsMatrix<T> tempVelVector = velVector;
    velVector = alpha4*(newDispVector - dispVector) + alpha5*tempVelVector + alpha6*accVector;
    accVector = alpha1*(newDispVector - dispVector) - alpha2*tempVelVector - alpha3*accVector;
    dispVector = newDispVector;
}

template <class T>
void gsElTimeIntegrator<T>::makeTimeStepNL(T timeStep)
{
    tStep = timeStep;
    T beta = 0.25;
    T gamma = 0.5;

    T alpha1 = 1./beta/pow(tStep,2);
    T alpha2 = 1./beta/tStep;
    T alpha3 = (1-2*beta)/2/beta;
    T alpha4 = gamma/beta/tStep;
    T alpha5 = 1 - gamma/beta;
    T alpha6 = (1 - gamma/beta/2)*tStep;

    gsElNewton<T> solver(*this,dispVector);
    solver.options().setInt("Verbosity",newton_verbosity::none);
    solver.solve();
    gsMatrix<T> tempVelVector = velVector;
    velVector = alpha4*(solver.solution() - dispVector) + alpha5*tempVelVector + alpha6*accVector;
    accVector = alpha1*(solver.solution() - dispVector) - alpha2*tempVelVector - alpha3*accVector;
    dispVector = solver.solution();
}

template <class T>
bool gsElTimeIntegrator<T>::assemble(const gsMatrix<T> & solutionVector)
{
    T beta = 0.25;
    T gamma = 0.5;

    T alpha1 = 1./beta/pow(tStep,2);
    T alpha2 = 1./beta/tStep;
    T alpha3 = (1-2*beta)/2/beta;
    T alpha4 = gamma/beta/tStep;
    T alpha5 = 1 - gamma/beta;
    T alpha6 = (1 - gamma/beta/2)*tStep;

    if (!stiffAssembler.assemble(solutionVector))
        return false;
    Base::m_system.matrix() = alpha1*massAssembler.matrix() + stiffAssembler.matrix();
    Base::m_system.matrix().makeCompressed();

    Base::m_system.rhs() = stiffAssembler.rhs() + massAssembler.matrix()*(alpha1*(dispVector-solutionVector) + alpha2*velVector + alpha3*accVector);
    return true;
}



} // namespace ends
