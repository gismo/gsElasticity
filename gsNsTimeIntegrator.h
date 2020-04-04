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
                       gsMassAssembler<T> & massAssembler_,
                       gsMultiPatch<T> * ALEvelocity = nullptr,
                       std::vector<std::pair<index_t,index_t> > * ALEpatches = nullptr);

    /// @brief Returns the list of default options for assembly
    static gsOptionList defaultOptions();
    /// set intial conditions
    void setSolutionVector(const gsMatrix<T> & solutionVector)
    {
        GISMO_ENSURE(solutionVector.rows() == stiffAssembler.numDofs(),"Wrong size of the solution vector: " + util::to_string(solutionVector.rows()) +
                     ". Must be: " + util::to_string(stiffAssembler.numDofs()));
        solVector = solutionVector;
        initialized = false;
    }
    /// set all fixed degrees of freedom
    virtual void setFixedDofs(const std::vector<gsMatrix<T> > & ddofs)
    {
        Base::setFixedDofs(ddofs);
        initialized = false;
    }
    /// make a time step according to a chosen scheme
    void makeTimeStep(T timeStep, bool ALE = false);
    /// assemble the linear system for the nonlinear solver
    virtual bool assemble(const gsMatrix<T> & solutionVector,
                          const std::vector<gsMatrix<T> > & fixedDoFs,
                          bool assembleMatrix = true);

    /// returns number of degrees of freedom
    virtual int numDofs() const { return stiffAssembler.numDofs(); }
    /// returns solution vector
    const gsMatrix<T> & solutionVector() const
    {
        GISMO_ENSURE(solVector.rows() == stiffAssembler.numDofs(),
                     "No initial conditions provided!");
        return solVector;
    }
    /// save solver state
    void saveState();
    /// recover solver state from the previously saved state
    void recoverState();
    /// number of iterations Newton's method required at the last time step; always 1 for IMEX
    index_t numberIterations() const { return numIters;}

protected:
    void initialize();
    /// time integraton schemes
    void implicitLinear();
    void implicitNonlinear();

protected:
    /// assembler object that generates the static system
    gsNsAssembler<T> & stiffAssembler;
    /// assembler object that generates the mass matrix
    gsMassAssembler<T> & massAssembler;
    /// initialization flag
    bool initialized;
    /// time step length
    T tStep;
    /// solution vector
    gsMatrix<T> solVector;
    using Base::m_system;
    using Base::m_options;
    using Base::m_ddof;

    /// IMEX stuff
    T oldTimeStep;
    gsMatrix<T> oldSolVector;

    /// Newton stuff
    gsMatrix<T> constRHS;
    index_t numIters;

    /// use ALE shift of the advective velocity
    bool flagALE;
    /// ALE velocity
    gsMultiPatch<T> * velocityALE;
    /// mapping between the geometry patches and the ALE velocity patches
    std::vector<std::pair<index_t,index_t> > * patchesALE;

    /// saved state
    bool hasSavedState;
    gsMatrix<T> velVecSaved;
    gsMatrix<T> oldVecSaved;
    gsMatrix<T> massRhsSaved;
    gsMatrix<T> stiffRhsSaved;
    gsSparseMatrix<T> stiffMatrixSaved;
    std::vector<gsMatrix<T> > ddofsSaved;
};

}

#ifndef GISMO_BUILD_LIB
#include GISMO_HPP_HEADER(gsNsTimeIntegrator.hpp)
#endif
