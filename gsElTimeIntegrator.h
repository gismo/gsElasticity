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

    /// assemble the linear system for
    virtual bool assemble(const gsMatrix<T> & solutionVector, bool assembleMatrix = true);

    virtual int numDofs() const { return stiffAssembler.numDofs(); }

    /// returns  vector of displacement DoFs
    const gsMatrix<T> & displacementVector() const {return dispVector;}

    /// sets scaling of Dirichlet BC used for linear system assembly
    virtual void setDirichletAssemblyScaling(T factor) { stiffAssembler.setDirichletAssemblyScaling(factor); }
    /// sets scaling of Dirichlet BC used for construction of the solution as a gsMultiPatch object
    virtual void setDirichletConstructionScaling(T factor) { stiffAssembler.setDirichletConstructionScaling(factor); }
    /// set scaling of the force loading (volume and surface loading)
    virtual void setForceScaling(T factor) { stiffAssembler.setForceScaling(factor); }

protected:
    /// time integration scheme coefficients
    T alpha1() {return 1./m_options.getReal("Beta")/pow(tStep,2); }
    T alpha2() {return 1./m_options.getReal("Beta")/tStep; }
    T alpha3() {return (1-2*m_options.getReal("Beta"))/2/m_options.getReal("Beta"); }
    T alpha4() {return m_options.getReal("Gamma")/m_options.getReal("Beta")/tStep; }
    T alpha5() {return 1 - m_options.getReal("Gamma")/m_options.getReal("Beta"); }
    T alpha6() {return (1-m_options.getReal("Gamma")/m_options.getReal("Beta")/2)*tStep; }

protected:
    /// assembler object that generates the static system
    gsElasticityAssembler<T> & stiffAssembler;
    /// assembler object that generates the mass matrix
    gsElMassAssembler<T> & massAssembler;
    /// Sparse matrix of the linear system to solve
    gsSparseMatrix<T> m_matrix;
    /// RHS vector of the linear system to solve
    gsMatrix<T> m_rhs;

    /// time step length
    T tStep;
    /// vector of displacement DoFs
    gsMatrix<T> dispVector;
    /// vector of velocity DoFs
    gsMatrix<T> velVector;
    /// vector of acceleration DoFs
    gsMatrix<T> accVector;

    using Base::m_options;
};

}
