/** @file gsElNewton.h

    @brief A class providing a nonlinear solver based on Newton's method.
    Supports incremental loading (ILS = incremental loading step),
    step size adaptivity and damping for convergence control.

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
#include <gsElasticity/gsElUtils.h>
#include <functional>

namespace gismo
{

template <class T>
class gsElBaseAssembler;

/** @brief A general nonlinear solver based on Newton's method.
 * An equation to solve is specified by an assembler class which
 * provides the following interfaces:
 *
 * int numDofs() const;
 * const gsSparseMatrix<T> & matrix() const;
 * const gsMatrix<T> & rhs() const;
 * void assembler(const gsMatrix<T> & solutionVector);
 * options().setReal("DirichletScaling",T);
 * options().setReal("ForceScaling",T);
 * .
 * Currently uses gsElBaseAssembler as a parent interface class. Potentially, can operate on gsAssembler.
*/
template <class T>
class gsElNewton
{
public:
    /// constructor without an initial guess. Assumes a zero initial guess which probably does not
    /// satisfy Dirichlet BC. Uses incremental loading.
    gsElNewton(gsElBaseAssembler<T> & assembler_);
    /// constructor with an initial guess. Assumes that the initial guess satisfies DBC
    /// and further updates are 0 at the Dirichlet boundary. Does not use incremental loading.
    gsElNewton(gsElBaseAssembler<T> & assembler_, const gsMatrix<T> & initialSolVector);
    /// default option list. used for initialization
    static gsOptionList defaultOptions();
    /// get options list to read or set parameters
    gsOptionList & options() { return m_options; }
    /// solution procedure
    void solve();
    /// returns the solution vector
    const gsMatrix<T> & solution() const { return solVector; }
    /// return solver status as a string
    std::string status();
    /// reset the solver state
    void reset();

    void setPreProcessingFunction(std::function<void(const gsMatrix<T> &)> f)
    {
        preProcessingFunction = f;
    }

    void setPostProcessingFunction(std::function<void(const gsMatrix<T> &)> f)
    {
        postProcessingFunction = f;
    }

protected:
    /// computes update of the solution
    bool computeUpdate();
    /// solution procedure without an initial guess. Assumes 0 initial guess
    /// that does not satisfy Dirichlet BC. Uses incremental loading
    void solveNoGuess();
    /// solution procedure with an initial guess. Assumes that the initial guess
    /// satisfies Dirichlet BC. Does not use incremental loading.
    void solveWithGuess();

protected:
    /// assembler object that generates the linear system
    gsElBaseAssembler<T> & assembler;
    /// solution vector
    gsMatrix<T> solVector;
    bool initialGuess;
    /// ---- status variables ----- ///
    index_t numIterations; /// number of Newton's iterations performed at the current ILS
    index_t incStep; /// current incremental loading step
    newton_status m_status;  /// status of the solver (converged, interrupted, working)
    T residualNorm; /// norm of the residual vector
    T initResidualNorm; /// norm of the residual vector at the beginning of the loop
    T updateNorm; /// norm of the update vector
    T initUpdateNorm; /// norm of the update vector at the beginning of the loop
    /// option list
    gsOptionList m_options;

    std::function<void(const gsMatrix<T> &)> preProcessingFunction;
    std::function<void(const gsMatrix<T> &)> postProcessingFunction;
};

} // namespace ends
