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

#include <gsElasticity/src/gsBaseAssembler.h>
#include <gsElasticity/src/gsBaseUtils.h>
#include <gsElasticity/src/gsElasticityAssembler.h>
#include <gsElasticity/src/gsMassAssembler.h>

namespace gismo
{

/** @brief Time integation for equations of dynamic elasticity with implicit schemes
*/
template <class T>
class gsElTimeIntegrator : public gsBaseAssembler<T>
{
public:
    typedef gsBaseAssembler<T> Base;
    /// constructor method. requires a gsElasticityAssembler for construction of the static linear system
    /// and a gsMassAssembler for the mass matrix
    gsElTimeIntegrator(gsElasticityAssembler<T> & stiffAssembler_,
                       gsMassAssembler<T> & massAssembler_);

    /// @brief Returns the list of default options for assembly
    static gsOptionList defaultOptions();

    /// set intial conditions
    void setDisplacementVector(const gsMatrix<T> & displacementVector)
    {
        GISMO_ENSURE(displacementVector.rows() == massAssembler.numDofs(),
                     "Wrong size of the displacement vector: " + util::to_string(displacementVector.rows()) +
                     ". Must be: " + util::to_string(massAssembler.numDofs()));
        solVector.middleRows(0,massAssembler.numDofs()) = displacementVector;
        initialized = false;
    }

    void setVelocityVector(const gsMatrix<T> & velocityVector)
    {
        GISMO_ENSURE(velocityVector.rows() == massAssembler.numDofs(),
                     "Wrong size of the velocity vector: " + util::to_string(velocityVector.rows()) +
                     ". Must be: " + util::to_string(massAssembler.numDofs()));
        velVector = velocityVector;
        initialized = false;
    }

    /// make a time step according to a chosen scheme
    void makeTimeStep(T timeStep);

    /// assemble the linear system for the nonlinear solver
    using Base::assemble;
    virtual bool assemble(const gsMatrix<T> & solutionVector,
                          const std::vector<gsMatrix<T> > & fixedDoFs);

    /// return the number of free degrees of freedom
    virtual int numDofs() const;

    /// returns complete solution vector (displacement + possibly pressure)
    const gsMatrix<T> & solutionVector() const { return solVector; }

    /// returns vector of displacement DoFs
    //const gsMatrix<T> & displacementVector() const { return solVector.middleRows(0,massAssembler.numDofs()); }

    /// returns vector of velocity DoFs
    const gsMatrix<T> & velocityVector() const { return velVector; }

    /// save solver state
    void saveState();

    /// recover solver state from saved state
    void recoverState();

    /// number of iterations Newton's method required at the last time step
    index_t numberIterations() const { return numIters;}

    /// construct displacement using the stiffness assembler
    void constructSolution(gsMultiPatch<T> & displacement) const;

    /// construct displacement and pressure (if applicable) using the stiffness assembler
    using Base::constructSolution;
    void constructSolution(gsMultiPatch<T> & displacement, gsMultiPatch<T> & pressure) const;

    /// assemblers' accessors
    gsBaseAssembler<T> & mAssembler() { return massAssembler; }
    gsBaseAssembler<T> & assembler() { return stiffAssembler; }

protected:
    void initialize();

    /// time integraton schemes
    gsMatrix<T> implicitLinear();
    gsMatrix<T> implicitNonlinear();

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
    gsMassAssembler<T> & massAssembler;
    /// initialization flag
    bool initialized;
    /// time step length
    T tStep;
    /// vector of displacement DoFs ( + possibly pressure)
    gsMatrix<T> solVector;
    /// vector of velocity DoFs
    gsMatrix<T> velVector;
    /// vector of acceleration DoFs
    gsMatrix<T> accVector;
    using Base::m_system;
    using Base::m_options;
    using Base::m_ddof;
    /// number of iterations Newton's method took to converge at the last time step
    index_t numIters;
    /// saved state
    bool hasSavedState;
    gsMatrix<T> solVecSaved;
    gsMatrix<T> velVecSaved;
    gsMatrix<T> accVecSaved;
    std::vector<gsMatrix<T> > ddofsSaved;
    /// temporary objects for memory efficiency
    gsMatrix<T> newSolVector, oldVelVector, dispVectorDiff;
    gsSparseMatrix<T> tempMassBlock;
};

}

#ifndef GISMO_BUILD_LIB
#include GISMO_HPP_HEADER(gsElTimeIntegrator.hpp)
#endif
