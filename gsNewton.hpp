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
      initialGuess(false),
      m_options(defaultOptions()),
      preProcessingFunction([](const gsMatrix<T> &){}),
      postProcessingFunction([](const gsMatrix<T> &){})
{
    solVector.setZero(assembler.numDofs(),1);
    reset();
}

template <class T>
gsNewton<T>::gsNewton(gsBaseAssembler<T> & assembler_,
                          const gsMatrix<T> & initialSolVector)
    : assembler(assembler_),
      solVector(initialSolVector),
      initialGuess(true),
      m_options(defaultOptions()),
      preProcessingFunction([](const gsMatrix<T> &){}),
      postProcessingFunction([](const gsMatrix<T> &){})
{
    reset();
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
    /// incremental loading
    opt.addInt("NumIncSteps","Number of incremental loading steps",1);
    opt.addInt("MaxItersInter","Maximum number of Newton's iterations at intermediate loading steps. Same as MaxIter if < 1.",0);
    /// additional setting
    opt.addInt("Verbosity","Amount of information printed to the terminal: none, some, all",newton_verbosity::none);
    return opt;
}

template <class T>
void gsNewton<T>::solve()
{
    if (initialGuess)
        solveWithGuess();
    else
        solveNoGuess();
}

template <class T>
void gsNewton<T>::solveNoGuess()
{
    T stepSize = 1./m_options.getInt("NumIncSteps");
    for (index_t s = 0; s < m_options.getInt("NumIncSteps"); ++s)
    {
        if (m_options.getInt("Verbosity") != newton_verbosity::none)
            gsInfo << "Load: " << s*stepSize*100 << "% -> " << (s+1)*stepSize*100 << "%\n";
        reset();

        //scale load
        assembler.setForceScaling((s+1)*stepSize);
        // current solution satisfies a fraction s/numIncSteps of the Dirichlet BC
        assembler.setDirichletConstructionScaling(s*stepSize);
        // next update should advance the solution at the Dirichlet boundary by 1/numIncSteps of the Dirichlet BC
        assembler.setDirichletAssemblyScaling(stepSize);
        if (!computeUpdate())
            goto abort;

        // current solution satisfies a fraction (s+1)/numIncSteps of the Dirichlet BC
        assembler.setDirichletConstructionScaling((s+1)*stepSize);
        // further updates should be 0 at the Dirichelt boundary
        assembler.setDirichletAssemblyScaling(0.);

        index_t maxIter = (s == m_options.getInt("NumIncSteps") - 1 ?
                               m_options.getInt("MaxIters") : m_options.getInt("MaxItersInter"));
        if (maxIter < 1)
            maxIter = m_options.getInt("MaxIters");
        while (m_status == newton_status::working)
        {
            if (numIterations >= maxIter)
            {
                m_status = newton_status::interrupted;
                break;
            }
            if (!computeUpdate())
                goto abort;

            if (residualNorm < m_options.getReal("AbsTol") ||
                updateNorm < m_options.getReal("AbsTol") ||
                residualNorm/initResidualNorm < m_options.getReal("RelTol") ||
                updateNorm/initUpdateNorm < m_options.getReal("RelTol"))
                m_status = newton_status::converged;
        }

        if (m_options.getInt("Verbosity") != newton_verbosity::none)
            gsInfo << status() << std::endl;
    }
    abort:;
}

template <class T>
void gsNewton<T>::solveWithGuess()
{
    // assuming that the initial guess satisfies Dirichlet BC,
    assembler.setDirichletConstructionScaling(1.);
    // further updates should be 0 at the Dirichelt boundary
    assembler.setDirichletAssemblyScaling(0.);

    while (m_status == newton_status::working)
    {
        if (!computeUpdate())
            goto abort;
        if (residualNorm < m_options.getReal("AbsTol") ||
            updateNorm < m_options.getReal("AbsTol") ||
            residualNorm/initResidualNorm < m_options.getReal("RelTol") ||
            updateNorm/initUpdateNorm < m_options.getReal("RelTol"))
            m_status = newton_status::converged;
        else if (numIterations == m_options.getInt("MaxIters"))
            m_status = newton_status::interrupted;
    }

    if (m_options.getInt("Verbosity") != newton_verbosity::none)
        gsInfo << status() << std::endl;
    abort:;
}

template <class T>
bool gsNewton<T>::computeUpdate()
{
    preProcessingFunction(solVector);

    if (!assembler.assemble(solVector))
    {
        m_status = newton_status::bad_solution;
        if (m_options.getInt("Verbosity") != newton_verbosity::none)
            gsInfo << status() << std::endl;
        return false;
    }

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

    if (numIterations == 0)
    {
        initUpdateNorm = updateNorm;
        initResidualNorm = residualNorm;
    }

    numIterations++;
    if (m_options.getInt("Verbosity") == newton_verbosity::all)
        gsInfo << status() << std::endl;

    postProcessingFunction(solVector);
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
    else if (m_status == newton_status::working)
        statusString = "It: " + util::to_string(numIterations) +
                 ", resAbs: " + util::to_string(residualNorm) +
                 ", resRel: " + util::to_string(residualNorm/initResidualNorm) +
                 ", updAbs: " + util::to_string(updateNorm) +
                 ", updRel: " + util::to_string(updateNorm/initUpdateNorm);
    else
        statusString = "Newton's method was interrupted due to an invalid solution.";
    return statusString;
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

} // namespace ends
