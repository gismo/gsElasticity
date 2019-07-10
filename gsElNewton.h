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

namespace gismo
{

template <class T>
class gsElBaseAssembler;

struct newton_verbosity
{
    enum verbosity
    {
        none = 0,  /// no output
        some = 1,  /// only essential output
        all = 2  /// output everything
    };
};

template <class T>
class gsElNewton
{
public:
    /// constructor without an initial guess
    gsElNewton(gsElBaseAssembler<T> & assembler_);


    gsElNewton(gsElBaseAssembler<T> & assembler_, const gsMatrix<T> & initialSolVector);

    static gsOptionList defaultOptions();
    /// get options list to read or set parameters
    gsOptionList & options() { return m_options; }
    /// standard solution procedure; supports incremental loading
    void solve();
    /// returns the solution vector
    const gsMatrix<T> & solution() const { return solVector; }

protected:
    void computeUpdate(bool initUpdate);

    void printStatus();

protected:
    /// assembler object that generates the linear system
    gsElBaseAssembler<T> & assembler;
    /// solution vector
    gsMatrix<T> solVector;

    /// status variables
    index_t numIterations; /// number of Newton's iterations performed at the current ILS
    index_t incStep; /// current incremental loading step
    bool converged;  /// convergence status at the current ILS
    T residualNorm;
    T initResidualNorm; /// residual norm at the beginning of the current ILS
    T updateNorm;
    T initUpdateNorm; /// update vector norm at the beginning of the current ILS

    gsOptionList m_options;
};

} // namespace ends
