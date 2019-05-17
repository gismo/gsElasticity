/** @file gsElasticityNewtonDeLuxe.hpp

    @brief A class providing Newton's method for nonlinear elasticity in mixed formulation

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s):
        O. Weeger    (2012 - 2015, TU Kaiserslautern),
        A.Shamanskiy (2016 - ...., TU Kaiserslautern)
*/

#pragma once

#include <gsElasticity/gsElasticityNewtonDeLuxe.h>
#include <gsElasticity/gsElasticityAssembler.h>
#include <gsCore/gsField.h>

namespace gismo
{

template <class T>
gsElasticityNewtonDeLuxe<T>::gsElasticityNewtonDeLuxe(gsElasticityAssembler<T> & elasticityAssembler)
    : assembler(elasticityAssembler),
      m_options(defaultOptions()) { }

template <class T>
gsOptionList gsElasticityNewtonDeLuxe<T>::defaultOptions()
{
    gsOptionList opt;
    /// stopping creteria
    opt.addInt("MaxIter","Maximum number of iterations per loading step",50);
    opt.addReal("AbsTol","Absolute tolerance for the convergence cretiria",1e-12);
    opt.addReal("RelTol","Relative tolerance for the stopping criteria",1e-6);
    /// incremental loading
    opt.addInt("NumIncStep","Number of incremental loading steps",1);
    opt.addInt("MaxIterNotLast","Maximum number of Newton's iterations per incremental loading step (not last!); "
                                "always full convergence at the last incremental loading step; "
                                "full convergence at every incremental loading step if < 1",0);
    opt.addSwitch("BijectivityCheck","Check bijectivity; stop the solution process if unable to preserve bijectivity",false);
    /// additional setting
    opt.addInt("Verbosity","Amount of information printed to the terminal: none, some, all",newtonVerbosity2::all);
    opt.addInt("Save","Amount of intermediate soluton information saved: "
                      "only the final solution, "
                      "the first and the last displacement fields at every incremental loading step, "
                      "all",newtonSave2::onlyFinal);
    return opt;
}

template <class T>
void gsElasticityNewtonDeLuxe<T>::solve()
{
    T absTol = m_options.getReal("AbsTol");
    T relTol = m_options.getReal("RelTol");
    bijective = true;
    index_t numIncSteps = m_options.getInt("NumIncStep");
    T stepSize = 1./numIncSteps;
    // save computed by the assembler Dirichet DoFs for further use
    std::vector<gsMatrix<T> > ddof(assembler.patches().dim());
    for (index_t d = 0; d < assembler.patches().dim(); ++d)
        ddof[d] = assembler.fixedDofs(d);

    for (index_t s = 0; s < numIncSteps; ++s)
    {
        // reset the status of Newton's method
        numIterations = 0;
        converged = false;

        if (m_options.getInt("Verbosity") != newtonVerbosity2::none)
            gsInfo << "Load: " << s*stepSize*100 << "% -> " << (s+1)*stepSize*100 << "%\n";
        // set load scaling for RHS and Neumann BC
        assembler.options().setReal("ForceScaling",(s+1)*stepSize);
        // set load scaling for Dirichlet BC;
        // a temporary variable is necessary because of the memory swap
        for (index_t d = 0; d < assembler.patches().dim(); ++d)
        {
            gsMatrix<T> tempDDof = stepSize * ddof[d];
            assembler.setFixedDofVector(tempDDof,d);
        }
        computeUpdate(true);
        if (!bijective)
        {
            if (m_options.getInt("Verbosity") != newtonVerbosity2::none)
                gsInfo << "Interrupted due to bijectivity violation!\n";
            goto abort;
        }
        // set Dirichlet BC to zero for further Newton's iterations at this ILS
        assembler.homogenizeFixedDofs(-1);
        // if this it not the last ILS, it can be interrupted earlier (to save computation)
        index_t maxNumIter = (s != numIncSteps-1) ? m_options.getInt("MaxIterNotLast") : m_options.getInt("MaxIter");
        while (!converged && numIterations < maxNumIter)
        {
            computeUpdate(false);
            if (!bijective)
            {
                if (m_options.getInt("Verbosity") != newtonVerbosity2::none)
                    gsInfo << "Interrupted due to bijectivity violation!\n";
                goto abort;
            }

            if (residualNorm < absTol || updateNorm < absTol ||
                residualNorm/initResidualNorm < relTol || updateNorm/initUpdateNorm < relTol)
                converged = true;
        }

        if (m_options.getInt("Verbosity") != newtonVerbosity2::none)
        {
            if (converged)
                gsInfo << "Newton's method converged after " << numIterations << " iterations\n";
            else
                gsInfo << "Newton's method interrupted after " << maxNumIter << " iterations\n";
        }
    }
    abort:;
}

template <class T>
void gsElasticityNewtonDeLuxe<T>::computeUpdate(bool initUpdate)
{
    if (displacements.empty()) // no previous displacement field
        assembler.assemble();
    else // use previous displacement field to assemble the problem
        assembler.assemble(displacements.back(),pressures.back());

    gsSparseSolver<>::SimplicialLDLT solver(assembler.matrix());
    gsVector<T> solVector = solver.solve(assembler.rhs());
    gsMultiPatch<T> incDisplacement, incPressure;
    assembler.constructSolution(solVector,incDisplacement,incPressure);


    if (m_options.getSwitch("BijectivityCheck") || assembler.options().getInt("MaterialLaw") == material_law::neo_hooke_ln)
        bijectivityCheck(incDisplacement);

    saveSolution(incDisplacement,incPressure);

    // compute norm for the stopping criteria
    updateNorm = solVector.norm();
    residualNorm = assembler.rhs().norm();
    // for the first Newton's iteraion of this ILS, save the initial norm for relative error
    if (initUpdate)
    {
        initUpdateNorm = updateNorm;
        initResidualNorm = residualNorm;
    }
    printStatus();
    numIterations++;
}

template <class T>
void gsElasticityNewtonDeLuxe<T>::printStatus()
{
    if (m_options.getInt("Verbosity") == newtonVerbosity2::all)
        gsInfo << "Iteration: " << numIterations
               << ", resAbs: " << residualNorm
               << ", resRel: " << residualNorm/initResidualNorm
               << ", updAbs: " << updateNorm
               << ", updRel: " << updateNorm/initUpdateNorm << std::endl;
}

template <class T>
void gsElasticityNewtonDeLuxe<T>::saveSolution(const gsMultiPatch<T> & incDisplacement,
                                         const gsMultiPatch<T> & incPressure)
{
    if (displacements.empty())
    {
        displacements.push_back(incDisplacement);
        pressures.push_back(incPressure);
    }
    else
        if (m_options.getInt("Save") == newtonSave2::onlyFinal)
            for (size_t p = 0; p < displacements.back().nPatches(); ++p)
            {
                displacements.back().patch(p).coefs() += incDisplacement.patch(p).coefs();
                pressures.back().patch(p).coefs() += incPressure.patch(p).coefs();
            }
        else if (m_options.getInt("Save") == newtonSave2::firstAndLastPerIncStep)
            if (numIterations == 0 || numIterations == 1)
            {
                displacements.push_back(incDisplacement);
                pressures.push_back(incPressure);
                for (size_t p = 0; p < displacements.back().nPatches(); ++p)
                {
                    displacements.back().patch(p).coefs() += displacements[displacements.size()-2].patch(p).coefs();
                    pressures.back().patch(p).coefs() += pressures[pressures.size()-2].patch(p).coefs();
                }
            }
            else
                for (size_t p = 0; p < displacements.back().nPatches(); ++p)
                {
                    displacements.back().patch(p).coefs() += incDisplacement.patch(p).coefs();
                    pressures.back().patch(p).coefs() += incPressure.patch(p).coefs();
                }
        else if (m_options.getInt("Save") == newtonSave2::all)
        {
            displacements.push_back(incDisplacement);
            pressures.push_back(incPressure);
            for (size_t p = 0; p < displacements.back().nPatches(); ++p)
            {
                displacements.back().patch(p).coefs() += displacements[displacements.size()-2].patch(p).coefs();
                pressures.back().patch(p).coefs() += pressures[pressures.size()-2].patch(p).coefs();
            }
        }
}

template <class T>
void gsElasticityNewtonDeLuxe<T>::bijectivityCheck(const gsMultiPatch<T> & incDisplacement)
{
    gsMultiPatch<T> tempDisplacement(incDisplacement);
    if (!displacements.empty())
        for (size_t p = 0; p < displacements.back().nPatches(); ++p)
            tempDisplacement.patch(p).coefs() += displacements.back().patch(p).coefs();
    index_t corruptedPatch = assembler.checkSolution(tempDisplacement);
    bijective = corruptedPatch == -1 ? true : false;
}

} // namespace ends
