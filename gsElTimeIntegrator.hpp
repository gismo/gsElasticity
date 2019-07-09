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
}

template <class T>
void gsElTimeIntegrator<T>::makeTimeStep(T timeStep)
{
    T beta = 0.25;
    T gamma = 0.5;

    T alpha1 = 1./beta/pow(timeStep,2);
    T alpha2 = 1./beta/timeStep;
    T alpha3 = (1-2*beta)/2/beta;

    m_matrix = alpha1*massAssembler.matrix() + stiffAssembler.matrix();
    m_rhs = massAssembler.matrix()*(alpha1*dispVector + alpha2*velVector + alpha3*accVector) + stiffAssembler.rhs();
    gsSparseSolver<>::SimplicialLDLT solver(m_matrix);
    dispVector = solver.solve(m_rhs);
}

template <class T>
void gsElTimeIntegrator<T>::assemble()
{

}

} // namespace ends
