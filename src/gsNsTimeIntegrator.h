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

#include <gsElasticity/src/gsBaseAssembler.h>
#include <gsElasticity/src/gsBaseUtils.h>
#include <gsElasticity/src/gsNsAssembler.h>
#include <gsElasticity/src/gsMassAssembler.h>

namespace gismo
{

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
                       gsBoundaryInterface * interfaceALE2Flow = nullptr);

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
    using Base::setFixedDofs;
    virtual void setFixedDofs(const std::vector<gsMatrix<T> > & ddofs)
    {
        Base::setFixedDofs(ddofs);
        initialized = false;
    }

    /// make a time step according to a chosen scheme
    void makeTimeStep(T timeStep);

    /// assemble the linear system for the nonlinear solver
    using Base::assemble;
    virtual bool assemble(const gsMatrix<T> & solutionVector,
                          const std::vector<gsMatrix<T> > & fixedDoFs);

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

    /// construct the solution using the stiffness matrix assembler
    using Base::constructSolution;
    void constructSolution(gsMultiPatch<T> & velocity, gsMultiPatch<T> & pressure) const;

    /// assemblers' accessors
    gsBaseAssembler<T> & mAssembler();
    gsBaseAssembler<T> & assembler();

    /// get mapping between the flow domain patches and the ALE mapping patches (if only some patches of the flow domain are deformed)
    const gsBoundaryInterface & aleInterface() const {return *interface;}

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

    /// ALE velocity
    gsMultiPatch<T> * velocityALE;
    /// mapping between the geometry patches and the ALE velocity patches
    gsBoundaryInterface * interface;

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
