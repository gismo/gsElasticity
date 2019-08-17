/** @file gsNsTimeIntegrator.h

    @brief Provides time integration for incompressible Navier-Stokes equations.

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s):
        A.Shamanskiy (2016 - ...., TU Kaiserslautern)
*/

#pragma once

#include <gsElasticity/gsBaseAssembler.h>
#include <gsElasticity/gsBaseUtils.h>

namespace gismo
{

template <class T>
class gsNsAssembler;
template <class T>
class gsMassAssembler;

/** @brief Time integation for incompressible Navier-Stokes equations.
*/
template <class T>
class gsNsTimeIntegrator : public gsBaseAssembler<T>
{
public:
    typedef gsBaseAssembler<T> Base;
    /// constructor method. requires a gsNsAssembler for construction of the static linear system
    /// and a gsMassAssembler for the mass matrix
    gsNsTimeIntegrator(gsNsAssembler<T> & stiffAssembler_,
                       gsMassAssembler<T> & massAssembler_);

    /// @brief Returns the list of default options for assembly
    static gsOptionList defaultOptions();
    /// set intial conditions
    void setInitialSolution(const gsMatrix<T> & initialSolution) { solVector = initialSolution; }
    /// @brief Initialize the solver; execute before computing any time steps
    void initialize();
    /// make a time step according to a chosen scheme
    void makeTimeStep(T timeStep);
    /// assemble the linear system for the nonlinear solver
    virtual bool assemble(const gsMatrix<T> & solutionVector, bool assembleMatrix = true);

    virtual int numDofs() const { return stiffAssembler.numDofs(); }
    /// returns  vector of displacement DoFs
    const gsMatrix<T> & solutionVector() const {return solVector;}

    /// sets scaling of Dirichlet BC used for linear system assembly
    virtual void setDirichletAssemblyScaling(T factor) { stiffAssembler.setDirichletAssemblyScaling(factor); }
    /// sets scaling of Dirichlet BC used for construction of the solution as a gsMultiPatch object
    virtual void setDirichletConstructionScaling(T factor) { stiffAssembler.setDirichletConstructionScaling(factor); }
    /// set scaling of the force loading (volume and surface loading)
    virtual void setForceScaling(T factor) { stiffAssembler.setForceScaling(factor); }

protected:
    /// time integraton schemes
    gsMatrix<T> implicitNewton();
    gsMatrix<T> implicitOseen();

protected:
    /// assembler object that generates the static system
    gsNsAssembler<T> & stiffAssembler;
    /// assembler object that generates the mass matrix
    gsMassAssembler<T> & massAssembler;
    /// Sparse matrix of the linear system to solve
    gsSparseMatrix<T> m_matrix;
    /// RHS vector of the linear system to solve
    gsMatrix<T> m_rhs;

    /// time step length
    T tStep;
    /// vector of displacement DoFs
    gsMatrix<T> solVector;
    using Base::m_system;
    using Base::m_options;
};

}
