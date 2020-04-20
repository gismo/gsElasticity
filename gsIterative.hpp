/** @file gsIterative.hpp

    @brief Implementation of gsElIterative.

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s):
        O. Weeger    (2012 - 2015, TU Kaiserslautern),
        A.Shamanskiy (2016 - ...., TU Kaiserslautern)
*/

#pragma once

#include <gsElasticity/gsIterative.h>

#include <gsElasticity/gsBaseAssembler.h>

#include <sstream>

namespace gismo
{

template <class T>
gsIterative<T>::gsIterative(gsBaseAssembler<T> & assembler_)
    : assembler(assembler_),
      m_options(defaultOptions())
{
    solVector.setZero(assembler.numDofs(),1);
    fixedDoFs = assembler.allFixedDofs();
    for (index_t d = 0; d < (index_t)(fixedDoFs.size()); ++d)
        fixedDoFs[d].setZero();
    reset();
}

template <class T>
gsIterative<T>::gsIterative(gsBaseAssembler<T> & assembler_,
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
gsIterative<T>::gsIterative(gsBaseAssembler<T> & assembler_,
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
void gsIterative<T>::reset()
{
    m_status = solver_status::working;
    numIterations = 0;
    residualNorm = 0.;
    initResidualNorm = 1.;
    updateNorm = 0.;
    initUpdateNorm = 1.;
}

template <class T>
gsOptionList gsIterative<T>::defaultOptions()
{
    gsOptionList opt;
    /// linear solver
    opt.addInt("Solver","Linear solver to use",linear_solver::LU);
    /// stopping creteria
    opt.addInt("MaxIters","Maximum number of iterations per loop",50);
    opt.addReal("AbsTol","Absolute tolerance for the convergence cretiria",1e-12);
    opt.addReal("RelTol","Relative tolerance for the stopping criteria",1e-9);
    /// additional setting
    opt.addInt("Verbosity","Amount of information printed to the terminal: none, some, all",solver_verbosity::none);
    opt.addInt("IterType","Type of iteration: update or next/full",iteration_type::update);
    return opt;
}

template <class T>
void gsIterative<T>::solve()
{
    while (m_status == solver_status::working)
    {
        if (!compute())
        {
            m_status = solver_status::bad_solution;
            goto abort;
        }
        if (m_options.getInt("Verbosity") == solver_verbosity::all)
            gsInfo << status() << std::endl;
        if (residualNorm < m_options.getReal("AbsTol") ||
            updateNorm < m_options.getReal("AbsTol") ||
            residualNorm/initResidualNorm < m_options.getReal("RelTol") ||
            updateNorm/initUpdateNorm < m_options.getReal("RelTol"))
            m_status = solver_status::converged;
        else if (numIterations == m_options.getInt("MaxIters"))
            m_status = solver_status::interrupted;
    }

    abort:;
    if (m_options.getInt("Verbosity") != solver_verbosity::none)
        gsInfo << status() << std::endl;
}

template <class T>
bool gsIterative<T>::compute()
{
    // update mode: set Dirichlet BC to zero after the first iteration
    if (numIterations == 1 && m_options.getInt("IterType") == iteration_type::update)
        assembler.homogenizeFixedDofs(-1);

    if (!assembler.assemble(solVector,fixedDoFs))
        return false;

    gsVector<T> solutionVector;
    if (m_options.getInt("Solver") == linear_solver::LU)
    {
#ifdef GISMO_WITH_PARDISO
        gsSparseSolver<>::PardisoLU solver(assembler.matrix());
        solutionVector = solver.solve(assembler.rhs());
#else
        gsSparseSolver<>::LU solver(assembler.matrix());
        solutionVector = solver.solve(assembler.rhs());
#endif
    }
    if (m_options.getInt("Solver") == linear_solver::LDLT)
    {
#ifdef GISMO_WITH_PARDISO
        gsSparseSolver<>::PardisoLDLT solver(assembler.matrix());
        solutionVector = solver.solve(assembler.rhs());
#else
        gsSparseSolver<>::SimplicialLDLT solver(assembler.matrix());
        solutionVector = solver.solve(assembler.rhs());
#endif
    }
    if (m_options.getInt("Solver") == linear_solver::BiCGSTABDiagonal)
    {
        gsSparseSolver<>::BiCGSTABDiagonal solver(assembler.matrix());
        solutionVector = solver.solve(assembler.rhs());
    }
    if (m_options.getInt("Solver") == linear_solver::CGDiagonal)
    {
        gsSparseSolver<>::CGDiagonal solver(assembler.matrix());
        solutionVector = solver.solve(assembler.rhs());
    }

    if (m_options.getInt("IterType") == iteration_type::update)
    {
        updateNorm = solutionVector.norm();
        residualNorm = assembler.rhs().norm();
        solVector += solutionVector;
        // update fixed degrees fo freedom at the first iteration only (they are zero afterwards)
        if (numIterations == 0)
            for (index_t d = 0; d < (index_t)(fixedDoFs.size()); ++d)
                fixedDoFs[d] += assembler.fixedDofs(d);
    }
    else if (m_options.getInt("IterType") == iteration_type::next)
    {
        updateNorm = (solutionVector-solVector).norm();
        residualNorm = 1.; // residual is not defined
        solVector = solutionVector;
        // copy the fixed degrees of freedom
        if (numIterations == 0)
            fixedDoFs = assembler.allFixedDofs();
    }

    if (numIterations == 0)
    {
        initUpdateNorm = updateNorm;
        initResidualNorm = residualNorm;
    }
    numIterations++;

    return true;
}

template <class T>
std::string gsIterative<T>::status()
{
    std::string statusString;
    if (m_status == solver_status::converged)
        statusString = "Iterative solver converged after " +
                 util::to_string(numIterations) + " iteration(s).";
    else if (m_status == solver_status::interrupted)
        statusString = "Iterative solver was interrupted after " +
                util::to_string(numIterations) + " iteration(s).";
    else if (m_status == solver_status::bad_solution)
        statusString = "Iterative solver was interrupted after " +
                util::to_string(numIterations) + " iteration(s) due to an invalid solution";
    else if (m_status == solver_status::working)
        statusString = "It: " + util::to_string(numIterations) +
                 ", updAbs: " + util::to_string(updateNorm) +
                 ", updRel: " + util::to_string(updateNorm/initUpdateNorm) +
                 ", resAbs: " + util::to_string(residualNorm) +
                 ", resRel: " + util::to_string(residualNorm/initResidualNorm);
    return statusString;
}

template <class T>
void gsIterative<T>::setFixedDofs(const std::vector<gsMatrix<T> > & ddofs)
{
    GISMO_ENSURE(ddofs.size() >= fixedDoFs.size(), "Wrong size of the container with fixed DoFs: " + util::to_string(ddofs.size()) +
                 ". Must be at least: " + util::to_string(fixedDoFs.size()));

    for (short_t d = 0; d < (short_t)(fixedDoFs.size()); ++d)
    {
        GISMO_ENSURE(fixedDoFs[d].rows() == ddofs[d].rows(),"Wrong number of fixed DoFs for " + util::to_string(d) + "component: " +
                     util::to_string(ddofs[d].rows()) + ". Must be: " + util::to_string(fixedDoFs[d].rows()));
        fixedDoFs[d] = ddofs[d];
    }
}

template <class T>
void gsIterative<T>::saveState()
{
    solVecSaved = solVector;
    ddofsSaved = fixedDoFs;
}

template <class T>
void gsIterative<T>::recoverState()
{
    GISMO_ASSERT(solVecSaved.rows() == solVector.rows(),"No state saved!");
    solVector = solVecSaved;
    fixedDoFs = ddofsSaved;
}

} // namespace ends
