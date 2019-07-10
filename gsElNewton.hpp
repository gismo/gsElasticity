/** @file gsElNewton.hpp

    @brief Implementation of gsElNeton.

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s):
        O. Weeger    (2012 - 2015, TU Kaiserslautern),
        A.Shamanskiy (2016 - ...., TU Kaiserslautern)
*/

#pragma once

#include <gsElasticity/gsElNewton.h>

#include <gsElasticity/gsElBaseAssembler.h>

namespace gismo
{

template <class T>
gsElNewton<T>::gsElNewton(gsElBaseAssembler<T> & assembler_)
    : assembler(assembler_),
      m_options(defaultOptions())
{
    solVector.setZero(assembler.numDofs(),1);
}

template <class T>
gsElNewton<T>::gsElNewton(gsElBaseAssembler<T> & assembler_,
                          const gsMatrix<T> & initialSolVector)
    : assembler(assembler_),
      m_options(defaultOptions()),
      solVector(initialSolVector)
{

}

template <class T>
gsOptionList gsElNewton<T>::defaultOptions()
{
    gsOptionList opt;
    /// stopping creteria
    opt.addInt("MaxIter","Maximum number of iterations per loading step",50);
    opt.addReal("AbsTol","Absolute tolerance for the convergence cretiria",1e-12);
    opt.addReal("RelTol","Relative tolerance for the stopping criteria",1e-9);
    /// additional setting
    opt.addInt("Verbosity","Amount of information printed to the terminal: none, some, all",newton_verbosity::all);
    return opt;
}

template <class T>
void gsElNewton<T>::solve()
{
    // reset the status of Newton's method
    numIterations = 0;
    converged = false;

    computeUpdate(true);

    assembler.options().setReal("DirichletScaling",0.);

    while (!converged && numIterations < m_options.getInt("MaxIter"))
    {
        computeUpdate(false);
        if (residualNorm < m_options.getReal("AbsTol") ||
            updateNorm < m_options.getReal("AbsTol") ||
            residualNorm/initResidualNorm < m_options.getReal("RelTol") ||
            updateNorm/initUpdateNorm < m_options.getReal("RelTol"))
            converged = true;
    }
}

template <class T>
void gsElNewton<T>::computeUpdate(bool initUpdate)
{
    assembler.assemble(solVector);
    gsSparseSolver<>::SimplicialLDLT solver(assembler.matrix());
    gsVector<T> updateVector = solver.solve(assembler.rhs());

    updateNorm = updateVector.norm();
    residualNorm = assembler.rhs().norm();

    solVector += updateVector;

    if (initUpdate)
    {
        initUpdateNorm = updateNorm;
        initResidualNorm = residualNorm;
    }

    numIterations++;
    printStatus();
}

template <class T>
void gsElNewton<T>::printStatus()
{
    if (m_options.getInt("Verbosity") == newton_verbosity::all)
        gsInfo << "Iteration: " << numIterations
               << ", resAbs: " << residualNorm
               << ", resRel: " << residualNorm/initResidualNorm
               << ", updAbs: " << updateNorm
               << ", updRel: " << updateNorm/initUpdateNorm << std::endl;
}

} // namespace ends
