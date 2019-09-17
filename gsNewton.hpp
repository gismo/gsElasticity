/** @file gsNewton.hpp

    @brief Implementation of gsElNewton.

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s):
        O. Weeger    (2012 - 2015, TU Kaiserslautern),
        A.Shamanskiy (2016 - ...., TU Kaiserslautern)
*/

#pragma once

#include <gsElasticity/gsNewton.h>

#include <gsElasticity/gsBaseAssembler.h>

#include <sstream>

namespace gismo
{

template <class T>
gsNewton<T>::gsNewton(gsBaseAssembler<T> & assembler_)
    : assembler(assembler_),
      m_options(defaultOptions())
{
    solVector.setZero(assembler.numDofs(),1);
    fixedDoFs = assembler.allFixedDofs();
    for (index_t d = 0; d < index_t(fixedDoFs.size()); ++d)
        fixedDoFs[d].setZero();
    reset();
}

template <class T>
gsNewton<T>::gsNewton(gsBaseAssembler<T> & assembler_,
                      const gsMatrix<T> & initFreeDoFs)
    : assembler(assembler_),
      solVector(initFreeDoFs),
      m_options(defaultOptions())
{
    fixedDoFs = assembler.allFixedDofs();
    assembler.homogenizeFixedDofs(-1);
    reset();
}

template <class T>
gsNewton<T>::gsNewton(gsBaseAssembler<T> & assembler_,
                      const gsMatrix<T> & initFreeDoFs,
                      const std::vector<gsMatrix<T> > & initFixedDoFs)
    : assembler(assembler_),
      solVector(initFreeDoFs),
      fixedDoFs(initFixedDoFs),
      m_options(defaultOptions())
{
    reset();
}

template <class T>
void gsNewton<T>::reset()
{
    m_status = newton_status::working;
    numIterations = 0;
    residualNorm = 0.;
    initResidualNorm = 1.;
    updateNorm = 0.;
    initUpdateNorm = 1.;
}

template <class T>
gsOptionList gsNewton<T>::defaultOptions()
{
    gsOptionList opt;
    /// linear solver
    opt.addInt("Solver","Linear solver to use",linear_solver::SimplicialLDLT);
    /// stopping creteria
    opt.addInt("MaxIters","Maximum number of iterations per loop",50);
    opt.addReal("AbsTol","Absolute tolerance for the convergence cretiria",1e-12);
    opt.addReal("RelTol","Relative tolerance for the stopping criteria",1e-9);
    /// additional setting
    opt.addInt("Verbosity","Amount of information printed to the terminal: none, some, all",newton_verbosity::none);
    return opt;
}

template <class T>
void gsNewton<T>::solve()
{
    while (m_status == newton_status::working)
    {
        if (!computeUpdate())
        {
            m_status = newton_status::bad_solution;
            goto abort;
        }
        if (m_options.getInt("Verbosity") == newton_verbosity::all)
            gsInfo << status() << std::endl;
        if (residualNorm < m_options.getReal("AbsTol") ||
            updateNorm < m_options.getReal("AbsTol") ||
            residualNorm/initResidualNorm < m_options.getReal("RelTol") ||
            updateNorm/initUpdateNorm < m_options.getReal("RelTol"))
            m_status = newton_status::converged;
        else if (numIterations == m_options.getInt("MaxIters"))
            m_status = newton_status::interrupted;
    }

    abort:;
    if (m_options.getInt("Verbosity") != newton_verbosity::none)
        gsInfo << status() << std::endl;
} 

template <class T>
bool gsNewton<T>::computeUpdate()
{
    if (numIterations == 1) // set Dirichlet BC to zero after the first iteration
        assembler.homogenizeFixedDofs(-1);

    if (!assembler.assemble(solVector,fixedDoFs))
        return false;

    gsVector<T> updateVector;
    if (m_options.getInt("Solver") == linear_solver::BiCGSTABILUT)
    {
        gsSparseSolver<>::BiCGSTABILUT solver(assembler.matrix());
        updateVector = solver.solve(assembler.rhs());
    }
    if (m_options.getInt("Solver") == linear_solver::CGDiagonal)
    {
        gsSparseSolver<>::CGDiagonal solver(assembler.matrix());
        updateVector = solver.solve(assembler.rhs());
    }
    if (m_options.getInt("Solver") == linear_solver::LU)
    {
        gsSparseSolver<>::LU solver(assembler.matrix());
        updateVector = solver.solve(assembler.rhs());
    }
    if (m_options.getInt("Solver") == linear_solver::SimplicialLDLT)
    {
        gsSparseSolver<>::SimplicialLDLT solver(assembler.matrix());
        updateVector = solver.solve(assembler.rhs());
    }

    updateNorm = updateVector.norm();
    residualNorm = assembler.rhs().norm();

    solVector += updateVector;
    // update fixed degrees fo freedom at the first iteration only (they are zero afterwards)
    if (numIterations == 0)
        for (index_t d = 0; d < index_t(fixedDoFs.size()); ++d)
            fixedDoFs[d] += assembler.fixedDofs(d);

    if (numIterations == 0)
    {
        initUpdateNorm = updateNorm;
        initResidualNorm = residualNorm;
    }

    numIterations++;

    return true;
}

template <class T>
std::string gsNewton<T>::status()
{
    std::string statusString;
    if (m_status == newton_status::converged)
        statusString = "Newton's method converged after " +
                 util::to_string(numIterations) + " iteration(s).";
    else if (m_status == newton_status::interrupted)
        statusString = "Newton's method was interrupted after " +
                util::to_string(numIterations) + " iteration(s).";
    else if (m_status == newton_status::bad_solution)
        statusString = "Newton's method was interrupted after " +
                util::to_string(numIterations) + " iteration(s) due to an invalid solution";
    else if (m_status == newton_status::working)
        statusString = "It: " + util::to_string(numIterations) +
                 ", updAbs: " + util::to_string(updateNorm) +
                 ", updRel: " + util::to_string(updateNorm/initUpdateNorm) +
                 ", resAbs: " + util::to_string(residualNorm) +
                 ", resRel: " + util::to_string(residualNorm/initResidualNorm);
    return statusString;
}

template <class T>
void gsNewton<T>::setFixedDofs(const std::vector<gsMatrix<T> > & ddofs)
{
    GISMO_ENSURE(ddofs.size() == fixedDoFs.size(), "Wrong size of the container with fixed DoFs: " + util::to_string(ddofs.size()) +
                 ". Must be: " + util::to_string(fixedDoFs.size()));

    for (short_t d = 0; d < index_t(fixedDoFs.size()); ++d)
    {
        GISMO_ENSURE(fixedDoFs[d].rows() == ddofs[d].rows(),"Wrong number of fixed DoFs for " + util::to_string(d) + "component: " +
                     util::to_string(ddofs[d].rows()) + ". Must be: " + util::to_string(fixedDoFs[d].rows()));
        fixedDoFs[d] = ddofs[d];
    }
}


} // namespace ends
