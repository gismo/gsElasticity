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

#include <gsElasticity/gsBaseAssembler.h>
#include <gsElasticity/gsBaseUtils.h>

namespace gismo
{

template <class T>
class gsElasticityAssembler;
template <class T>
class gsMassAssembler;

/** @brief Time integation for equations of dynamic elasticity. Can perform explicit (with or without mass lamping)
 *         and implicit (linear and nonlinear) time steps.
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
    void setInitialDisplacement(const gsMatrix<T> & initialDisplacement) { dispVector = initialDisplacement; }
    void setInitialVelocity(const gsMatrix<T> & initialVelocity) { velVector = initialVelocity; }
    /// @brief Initialize the solver; execute before computing any time steps
    void initialize();
    /// make a time step according to a chosen scheme
    void makeTimeStep(T timeStep);
    /// assemble the linear system for the nonlinear solver
    virtual bool assemble(const gsMatrix<T> & solutionVector,
                          const std::vector<gsMatrix<T> > & fixedDoFs,
                          bool assembleMatrix = true);
    /// return the number of free degrees of freedom
    virtual int numDofs() const { return stiffAssembler.numDofs(); }
    /// returns vector of displacement DoFs
    const gsMatrix<T> & displacementVector() const {return dispVector;}
    /// returns vector of velocity DoFs
    const gsMatrix<T> & velocityVector() const {return velVector;}

protected:
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
    using Base::m_ddof;
};

}
