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

#include <gsElasticity/gsElasticityAssembler.h>

namespace gismo
{

/** @brief Provides Newton's method taylored for nonlinear elasticity
 *         a bit better than the parent class.
*/
template <class T>
class gsElasticityNewton
{
public:

    gsElasticityNewton(gsElasticityAssembler<T> & assembler)
        : m_assembler(assembler),
          m_verbose(true),
          m_doFirst(true),
          m_initResidueNorm(-1.),
          m_initUpdateNorm(-1.),
          m_firstDone(false)
    { }

    gsElasticityNewton(gsAssembler<T> & assembler, const gsMatrix<T> & solVector)
        : m_assembler(assembler),
          m_verbose(true),
          m_doFirst(false),
          m_initResidueNorm(-1.),
          m_firstDone(false)

    {
        m_assembler.constructSolution(solVector,m_curSolution);
        // bijectivity check for neo-Hookean material law
        if (static_cast<gsElasticityAssembler<T> &>(m_assembler).options().getInt("MaterialLaw") == 1)
        {
            index_t corruptedPatch = static_cast<gsElasticityAssembler<T> &>(m_assembler).checkSolution(m_curSolution);
            T J = (corruptedPatch == -1 ? 1 : -1);
            GISMO_ENSURE(J > 0,"Initial guess is not bijecive");
            (void)corruptedPatch;
        }
        m_updateNorm = m_initUpdateNorm = solVector.norm();
    }

    gsElasticityNewton(gsAssembler<T> & assembler, const gsMultiPatch<T> & solution)
        : m_assembler(assembler),
          m_verbose(true),
          m_doFirst(false),
          m_initResidueNorm(-1.),
          m_initUpdateNorm(-1.),
          m_firstDone(false)
    {
        m_curSolution.clear();

        for (index_t p = 0; p < solution.nPatches(); ++p)
            m_curSolution.addPatch(solution.patch(p).clone());
        // bijectivity check for neo-Hookean material law
        if (static_cast<gsElasticityAssembler<T> &>(m_assembler).options().getInt("MaterialLaw") == 1)
        {
            index_t corruptedPatch = static_cast<gsElasticityAssembler<T> &>(m_assembler).checkSolution(m_curSolution);
            T J = (corruptedPatch == -1 ? 1 : -1);
            GISMO_ENSURE(J > 0,"Initial guess is not bijecive");
            (void)corruptedPatch;
        }
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
        GISMO_ENSURE(m_firstDone,"Initial solution is not constructed\n");
        return m_linSolution;
    }

    /// \brief Returns the latest configuration
    const gsMultiPatch<T> & solution() const { return m_curSolution; }

    /// \brief Tells whether the Newton method converged
    bool converged() const {return m_converged;}

    /// \brief Returns the number of Newton iterations performed
    index_t numIterations() const { return m_numIterations; }

    /// \brief Sets verbosity of the solver
    void setVerbosity(bool verbose) { m_verbose = verbose; }

    /// brief Checks status of the solver; prints it if *m_verbose* is TRUE;
    /// return true if the solution is bijective
    bool status() const;

    /// \brief Set the maximum number of Newton iterations allowed
    void setMaxIterations(index_t nIter) {m_maxIterations = nIter;}

    /// \brief Set the tolerance for convergence
    void setTolerance(T tol) { m_tolerance = tol; }

protected:
    /// assembler object that generates the linear system for every Newton's step
    gsElasticityAssembler<T> & m_assembler;
    /// \brief If true prints useful information to the console. Default: true
    bool m_verbose;
    /// \brief If true, make firstIteration()
    bool m_doFirst;

    /// current displacement field
    gsMultiPatch<T> m_curSolution;

    gsVector<T> m_updateVector;
    /// \brief Solution vector of the linear problem, a.k.a. first iteration of Newton's method.
    /// Can be used for comparison.
    gsMultiPatch<T> m_linSolution;

    /// stopping criterion variables
    index_t m_maxIterations;
    T m_tolerance;
    T m_residueNorm;
    T m_updateNorm;
    /// \brief Initial value for the residual norm. Used by the stopping criterion
    T m_initResidueNorm;
    /// \brief Initial value for the update norm. Used by the stopping criterion
    T m_initUpdateNorm;

    /// status variables
    index_t m_numIterations;
    bool m_converged;
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

        if (m_updateNorm/m_initUpdateNorm < m_tolerance ||
            m_residueNorm/m_initResidueNorm < m_tolerance ||
            m_updateNorm < m_tolerance || m_residueNorm < m_tolerance)
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
    gsSparseSolver<>::LU solver(m_assembler.matrix());
    m_updateVector = solver.solve(m_assembler.rhs());

    m_assembler.constructSolution(m_updateVector,m_curSolution);
    m_assembler.constructSolution(m_updateVector,m_linSolution);

    // Homogenize Dirichlet dofs (values are now copied in m_curSolution)
    m_assembler.homogenizeFixedDofs(-1);

    m_residueNorm = m_initResidueNorm = m_assembler.rhs().norm();
    m_updateNorm = m_initUpdateNorm = m_updateVector.norm();

    bool bijective = status();
    (void)bijective;
    // bijectivity check for neo-Hookean material law
    if (static_cast<gsElasticityAssembler<T> &>(m_assembler).options().getInt("MaterialLaw") == 1)
        GISMO_ENSURE(bijective,"Failed to compute a valid initial guess. Reduce load or use incremental loading");


    m_firstDone = true;
    m_numIterations++;
}

template <class T>
void gsElasticityNewton<T>::nextIteration()
{   
    m_assembler.assemble(m_curSolution);
    gsSparseSolver<>::LU solver(m_assembler.matrix());
    m_updateVector = solver.solve(m_assembler.rhs());
    m_assembler.updateSolution(m_updateVector, m_curSolution);

    m_residueNorm = m_assembler.rhs().norm();
    if (m_initResidueNorm < 0)
        m_initResidueNorm = m_residueNorm;
    m_updateNorm = m_updateVector.norm();
    if (m_initUpdateNorm < 0)
        m_initUpdateNorm = m_updateNorm;

    bool bijective = status();

    // enforcing bijectivity for neo-Hookean material law
    if (static_cast<gsElasticityAssembler<T> &>(m_assembler).options().getInt("MaterialLaw") == 1)
        while (!bijective)
        {
            m_assembler.updateSolution(-1*m_updateVector, m_curSolution);
            m_updateVector *= 0.5;
            m_assembler.updateSolution(m_updateVector, m_curSolution);
            index_t corruptedPatch = static_cast<gsElasticityAssembler<T> &>(m_assembler).checkSolution(m_curSolution);
            bijective = (corruptedPatch == -1 ? true : false);
            gsInfo << "          Displacement update halved"
                   << ", J" << (corruptedPatch == -1 ? " > 0" : " < 0") << "\n";
        }

    m_numIterations++;
}

template <class T>
bool gsElasticityNewton<T>::status() const
{
    index_t corruptedPatch = static_cast<gsElasticityAssembler<T> &>(m_assembler).checkSolution(m_curSolution);

    if (m_verbose)
        gsInfo << "Iteration: " << m_numIterations
               << ", J" << (corruptedPatch == -1 ? " > 0" : " < 0")
               << ", residue: " << m_residueNorm
               << ", update norm: " << m_updateNorm <<"\n";

    return (corruptedPatch == -1 ? true : false);
}


} // namespace gismo
