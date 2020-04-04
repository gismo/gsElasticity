/** @file gsIterative.h

    @brief A class providing an iterative solver for nonlinear problems.

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s):
        O. Weeger    (2012 - 2015, TU Kaiserslautern),
        A.Shamanskiy (2016 - ...., TU Kaiserslautern)
*/

#pragma once

#include <gsIO/gsOptionList.h>
#include <gsElasticity/gsBaseUtils.h>
#include <functional>

namespace gismo
{

template <class T>
class gsBaseAssembler;
// TODO correct
/** @brief A general iterative solver for nonlinear problems.
 * An equation to solve is specified by an assembler class which
 * provides the following interfaces:
 *
 * int numDofs() const;
 * const gsSparseMatrix<T> & matrix() const;
 * const gsMatrix<T> & rhs() const;
 * void assemble(const gsMatrix<T> & solutionVector);
 * options().setReal("DirichletScaling",T);
 * options().setReal("ForceScaling",T);
 * .
 * Currently uses gsElBaseAssembler as a parent interface class. Potentially, can operate on gsAssembler.
*/
template <class T>
class gsIterative
{
public:
    typedef memory::shared_ptr<gsIterative> Ptr;
    typedef memory::unique_ptr<gsIterative> uPtr;
    /// constructor without an initial guess. Assumes a zero initial guess.
    gsIterative(gsBaseAssembler<T> & assembler_);
    /// constructor with an given initial free degrees of freedom.
    /// Fixed/Dirichlet degrees of freedom are taken from the assembler.
    /// fixed DoFs are given as a single vector arranged according to the function
    /// gsMatrix<> fixedDoFsAsVector() of gsBaseAssembler
    gsIterative(gsBaseAssembler<T> & assembler_,
             const gsMatrix<T> & initSolutionVector);
    /// constructor with an initial guess given as a combination of free and fixed/Dirichlet degrees of freedom.
    /// fixed DoFs are given as a single vector arranged according to the function
    /// gsMatrix<> fixedDoFsAsVector() of gsBaseAssembler
    gsIterative(gsBaseAssembler<T> & assembler_,
             const gsMatrix<T> & initSolutionVector,
             const std::vector<gsMatrix<T> > & initFixedDoFs);
    /// default option list. used for initialization
    static gsOptionList defaultOptions();
    /// get options list to read or set parameters
    gsOptionList & options() { return m_options; }
    /// solution procedure
    void solve();
    /// computes update or the next solution
    bool compute();
    /// returns the solution vector
    const gsMatrix<T> & solution() const { return solVector; }
    /// returns the fixed degrees of freedom
    const std::vector<gsMatrix<T> > & allFixedDofs() const { return fixedDoFs; }
    /// return solver status as a string
    std::string status();
    /// reset the solver state
    void reset();
    /// set all fixed degrees of freedom
    virtual void setFixedDofs(const std::vector<gsMatrix<T> > & ddofs);
    index_t numberIterations() const {return numIterations;}
    void setSolutionVector(const gsMatrix<T> & solutionVector) { solVector = solutionVector; }
    /// save solver state
    void saveState();
    /// recover solver state from saved state
    void recoverState();

protected:
    /// assembler object that generates the linear system
    gsBaseAssembler<T> & assembler;
    /// solution vector
    gsMatrix<T> solVector;
    /// current Dirichlet DoFs that the solution satisfies
    std::vector<gsMatrix<T> > fixedDoFs;
    /// ---- status variables ----- ///
    index_t numIterations; /// number of iterations performed
    solver_status m_status;  /// status of the solver (converged, interrupted, working)
    T residualNorm; /// norm of the residual vector
    T initResidualNorm; /// norm of the residual vector at the beginning of the loop
    T updateNorm; /// norm of the update vector
    T initUpdateNorm; /// norm of the update vector at the beginning of the loop
    /// option list
    gsOptionList m_options;

    gsMatrix<T> solVecSaved;
    std::vector<gsMatrix<T> > ddofsSaved;
};

} // namespace ends

#ifndef GISMO_BUILD_LIB
#include GISMO_HPP_HEADER(gsIterative.hpp)
#endif
