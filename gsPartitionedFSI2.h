/** @file gsPartitionedFSI.h

    @brief Partitioned FSI solver.

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s):
        A.Shamanskiy (2016 - ...., TU Kaiserslautern)
*/

#pragma once

#include <gsIO/gsOptionList.h>

namespace gismo
{

template <class T>
class gsNsTimeIntegrator;
template <class T>
class gsElTimeIntegrator;
template <class T>
class gsALE;
template <class T>
class gsMultiPatch;

template <class T>
class gsPartitionedFSI2
{
public:

    gsPartitionedFSI2(gsNsTimeIntegrator<T> & nsSolver,
                     gsMultiPatch<T> & velocity, gsMultiPatch<T> & pressure,
                     gsALE<T> & aleSolver,
                     gsMultiPatch<T> & aleDisplacement, gsMultiPatch<T> & aleVelocity);

    /// default option list. used for initialization
    static gsOptionList defaultOptions();

    /// get options list to read or set parameters
    gsOptionList & options() { return m_options; }

    /// make the next time step
    bool makeTimeStep(T timeStep);

    /// form a residual vector
    void formVector(const gsMultiPatch<T> & disp, gsMatrix<T> & vector);

    /// perform Aitken relaxation step
    void aitken(gsMultiPatch<T> & dispA, gsMultiPatch<T> & dispB,
                gsMultiPatch<T> & dispB2, gsMultiPatch<T> & dispC);

    /// number of iterations the solver took to converge at the last time step
    index_t numberIterations() { return numIter; }
    /// amount of time consumed by each component at the last time step
    T timeNS() { return nsTime; }
    T timeALE() { return aleTime; }
    /// aitken relaxation parameter used to at the last time step
    T aitkenOmega() { return omega;}
    /// FSI interface residual norm
    T residualNormAbs() { return absResNorm;}
    /// FSI interface relative residual norm
    T residualNormRel() { return absResNorm/initResNorm; }

protected:
    /// component solvers
    gsNsTimeIntegrator<T> & m_nsSolver;
    gsMultiPatch<T> & m_velocity;
    gsMultiPatch<T> & m_pressure;
    gsALE<T> & m_aleSolver;
    gsMultiPatch<T> & m_ALEdisplacment;
    gsMultiPatch<T> & m_ALEvelocity;
    /// option list
    gsOptionList m_options;
    /// status variables
    index_t numIter; // number of iterations at the last time step
    bool converged; // convergence flag
    T nsTime, aleTime; // component computational times
    T omega; // aitken relaxation parameter
    T absResNorm, initResNorm; // residual norms for convergence cretirion

};

} // namespace ends

#ifndef GISMO_BUILD_LIB
#include GISMO_HPP_HEADER(gsPartitionedFSI2.hpp)
#endif
