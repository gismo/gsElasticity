/** @file gsNewtonIETI.h

    @brief A version of the Newton iterator taylored for nonlinear elasticity.

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s):
        O. Weeger    (2012 - 2015, TU Kaiserslautern),
        A.Shamanskiy (2016 - ...., TU Kaiserslautern)
*/

#pragma once

#include <gsPde/gsNewtonIterator.h>
#include <gsTrilinos/gsTrilinos.h>

namespace gismo
{

/** @brief Provides Newton's method taylored for nonlinear elasticity
 *         a bit better than the parent class.
*/
template <class T>
class gsElasticityNewton : public gsNewtonIterator<T>
{
public:

    typedef gsNewtonIterator<T> Base;

    gsElasticityNewton(gsAssembler<T> & assembler)
        : gsNewtonIterator<T>(assembler),
          m_verbose(true),
          m_doFirst(true),
          m_initResidueNorm(-1.),
          m_initUpdateNorm(-1.),
          m_firstDone(false)
    { }

    gsElasticityNewton(gsAssembler<T> & assembler, gsMatrix<T> & solVector)
        : gsNewtonIterator<T>(assembler),
          m_verbose(true),
          m_doFirst(false),
          m_initResidueNorm(-1.),
          m_firstDone(false)

    {
        m_assembler.constructSolution(solVector,m_curSolution);
        m_updnorm = m_initUpdateNorm = solVector.norm();
    }


    /// \brief Applies Newton method and performs Newton iterations
    /// until convergence or maximum iterations.
    void solve();

    /// \brief Solves linear system obtained using linear elasticity
    /// in first step and computes residual
    void firstIteration();

    /// \brief Solves linear system in each iteration based on last
    /// solution and computes residual
    void nextIteration();

    /// \brief Returns the solution after the first iteration of Newton's method.
    const gsMultiPatch<T> & linearSolution() const
    {
        GISMO_ASSERT(m_firstDone,"Initial solution is not constructed\n");
        return m_linSolution;
    }

    /// \brief Sets verbosity of the solver
    void setVerbosity(bool verbose) { m_verbose = verbose; }

    /// brief Prints current status of the solver to the console
    void printStatus() const;

protected:
    // main
    using Base::m_assembler;
    /// \brief If true prints useful information to the console. Default: true
    bool m_verbose;
    /// \brief If true, make firstIteration()
    bool m_doFirst;

    // solution variables
    using Base::m_curSolution;
    gsVector<T> m_updateVector;
    /// \brief Solution vector of the linear problem, a.k.a. first iteration of Newton's method.
    /// Can be used for comparison.
    gsMultiPatch<T> m_linSolution;

    // stopping criterion variables
    using Base::m_maxIterations;
    using Base::m_tolerance;
    using Base::m_residue;
    using Base::m_updnorm;
    /// \brief Initial value for the residual norm. Used by the stopping criterion
    T m_initResidueNorm;
    /// \brief Initial value for the update norm. Used by the stopping criterion
    T m_initUpdateNorm;

    //status variables
    using Base::m_numIterations;
    using Base::m_converged;
    /// \brief Sets to true after firstIteration() is done
    bool m_firstDone;
};

template <class T>
void gsElasticityNewton<T>::solve()
{
    if (m_doFirst && !m_firstDone)
        firstIteration();

    // ----- Iterations start -----
    while (m_numIterations < m_maxIterations && !m_converged)
    {
        nextIteration();

        if (m_updnorm/m_initUpdateNorm < m_tolerance ||
            m_residue/m_initResidueNorm < m_tolerance ||
            m_updnorm < m_tolerance || m_residue < m_tolerance)
            m_converged = true;
    }

    if (m_verbose)
    {
        if (m_converged)
            gsInfo << "Newton's method converged after " << m_numIterations << " iterations.\n";
        else
            gsInfo << "Newton's method didn't converged after exceeding a maximum number of iterations: "
                   << m_maxIterations << ".\n";
    }

}

template <class T>
void gsElasticityNewton<T>::firstIteration()
{
    m_converged = false;
    m_numIterations = 0;

    m_assembler.assemble();
    Base::m_solver.compute(m_assembler.matrix());
    m_updateVector = Base::m_solver.solve(m_assembler.rhs());

    m_assembler.constructSolution(m_updateVector,m_curSolution);
    m_assembler.constructSolution(m_updateVector,m_linSolution);

    // Homogenize Dirichlet dofs (values are now copied in m_curSolution)
    m_assembler.homogenizeFixedDofs(-1);

    m_residue = m_initResidueNorm = m_assembler.rhs().norm();
    m_updnorm = m_initUpdateNorm = m_updateVector.norm();

    if (m_verbose)
        printStatus();

    m_firstDone = true;
    m_numIterations++;
}

template <class T>
void gsElasticityNewton<T>::nextIteration()
{
    m_assembler.assemble(m_curSolution);
    Base::m_solver.compute(m_assembler.matrix());
    m_updateVector = Base::m_solver.solve(m_assembler.rhs());
    m_assembler.updateSolution(m_updateVector, m_curSolution);

    m_residue = m_assembler.rhs().norm();
    if (m_initResidueNorm < 0)
        m_initResidueNorm = m_residue;
    m_updnorm = m_updateVector.norm();
    if (m_initUpdateNorm < 0)
        m_initUpdateNorm = m_updnorm;

    if (m_verbose)
        printStatus();

    m_numIterations++;
}

template <class T>
void gsElasticityNewton<T>::printStatus() const
{
    gsInfo << "Iteration: " << m_numIterations
           << ", residue: " << m_residue
           << ", update norm: " << m_updnorm <<"\n";
}


} // namespace gismo
