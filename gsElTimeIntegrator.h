/** @file gsElTimeIntegrator.h

    @brief Provides time integration for dynamical elasticity.

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s):
        A.Shamanskiy (2016 - ...., TU Kaiserslautern)
*/

#pragma once

#include <gsAssembler/gsAssembler.h>

namespace gismo
{

template <class T>
class gsElasticityAssembler;
template <class T>
class gsElMassAssembler;

template <class T>
class gsElTimeIntegrator
{
public:
    gsElTimeIntegrator(gsElasticityAssembler<T> & stiffAssembler_,
                       gsElMassAssembler<T> & massAssembler_);

    void makeTimeStep(T timeStep);

    const gsSparseMatrix<T> & matrix() const {return m_matrix; }

    const gsMatrix<T> & rhs() const {return m_rhs; }

    void assemble();



protected:
    /// assembler object that generates the static system
    gsElasticityAssembler<T> & stiffAssembler;
    /// assembler object that generates the mass matrix
    gsElMassAssembler<T> & massAssembler;
    /// Sparse matrix of the linear system to solve
    gsSparseMatrix<T> m_matrix;
    /// RHS vector of the linear system to solve
    gsMatrix<T> m_rhs;

    index_t timeStepNum;

    gsMatrix<T> dispVector;
    gsMatrix<T> velVector;
    gsMatrix<T> accVector;

};

}
