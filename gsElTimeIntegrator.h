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

#include <gsElasticity/gsElBaseAssembler.h>

#include <gsElasticity/gsElUtils.h>

namespace gismo
{

template <class T>
class gsElasticityAssembler;
template <class T>
class gsElMassAssembler;

/** @brief Time integation for equations of dynamic elasticity. Can perform explicit (with or without mass lamping)
 *         and implicit (linear and nonlinear) time steps.
*/
template <class T>
class gsElTimeIntegrator : public gsElBaseAssembler<T>
{
public:
    typedef gsElBaseAssembler<T> Base;
    /// constructor method. requires a gsElasticityAssembler for construction of the static linear system
    /// and a gsMassAssembler for the mass matrix
    gsElTimeIntegrator(gsElasticityAssembler<T> & stiffAssembler_,
                       gsElMassAssembler<T> & massAssembler_);

    /// @brief Returns the list of default options for assembly
    static gsOptionList defaultOptions();

    void makeTimeStep(T timeStep);

    void makeTimeStepNL(T timeStep);

    virtual bool assemble(const gsMatrix<T> & solutionVector);


    virtual int numDofs() const { return stiffAssembler.numDofs(); }

    const gsMatrix<T> & displacementVector() const {return dispVector;}


    virtual void setDirichletAssemblyScaling(T factor)
    {
        stiffAssembler.setDirichletAssemblyScaling(factor);
    }

    virtual void setDirichletConstructionScaling(T factor)
    {
        stiffAssembler.setDirichletConstructionScaling(factor);
    }

    virtual void setForceScaling(T factor)
    {
        stiffAssembler.setForceScaling(factor);
    }

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
    T tStep;

    gsMatrix<T> dispVector;
    gsMatrix<T> velVector;
    gsMatrix<T> accVector;

    using Base::m_options;
};

}
