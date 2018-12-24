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
        : gsNewtonIterator<T>(assembler)
    { }

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
    const gsMultiPatch<T> & linearSolution() { return m_linSolution; }

protected:

    using Base::m_curSolution;
    using Base::m_residue;
    using Base::m_updnorm;
    using Base::m_numIterations;
    using Base::m_maxIterations;
    using Base::m_tolerance;
    using Base::m_converged;
    using Base::m_updateVector;
    using Base::m_assembler;

    /// \brief Solution vector of the linear problem, a.k.a. first iteration of Newton's method.
    /// Can be used for comparison.
    gsMultiPatch<T> m_linSolution;

};

template <class T>
void gsElasticityNewton<T>::solve()
{
    firstIteration();

    const T initResidue = m_residue;
    const T initUpdate = m_updnorm;

    // ----- Iterations start -----
    for (m_numIterations = 1; m_numIterations < m_maxIterations; ++m_numIterations)
    {
        nextIteration();

        if (abs(m_updnorm/initUpdate)  < m_tolerance ||
            abs(m_residue/initResidue) < m_tolerance)
        {
            m_converged = true;
            break;
        }
    }
}

template <class T>
void gsElasticityNewton<T>::firstIteration()
{
    m_converged = false;

    m_assembler.assemble();
     Base::m_solver.compute(m_assembler.matrix());
    m_updateVector = Base::m_solver.solve(m_assembler.rhs());

    index_t numUnk = m_assembler.system().numUnknowns();
    gsVector<index_t> unknowns(numUnk);
    for (index_t d = 0; d < numUnk; ++d)
        unknowns.at(d) = d;
    m_assembler.constructSolution(m_updateVector,m_curSolution,unknowns);
    m_assembler.constructSolution(m_updateVector,m_linSolution,unknowns);

    // Homogenize Dirichlet dofs (values are now copied in m_curSolution)
    m_assembler.homogenizeFixedDofs(-1);

    m_residue = m_assembler.rhs().norm();
    m_updnorm = m_updateVector.norm();

    gsInfo << "Iteration: " << 0
           << ", residue: " << m_residue
           << ", update norm: " << m_updnorm <<"\n";
}

template <class T>
void gsElasticityNewton<T>::nextIteration()
{
    m_assembler.assemble(m_curSolution);
    Base::m_solver.compute(m_assembler.matrix());
    m_updateVector = Base::m_solver.solve(m_assembler.rhs());
    m_assembler.updateSolution(m_updateVector, m_curSolution);

    m_residue = m_assembler.rhs().norm();
    m_updnorm = m_updateVector.norm();

    gsInfo << "Iteration: " << m_numIterations
           << ", residue: " << m_residue
           << ", update norm: " << m_updnorm <<"\n";
}


} // namespace gismo
