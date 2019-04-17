/** @file gsElasticityNewton2.hpp

    @brief A class providing Newton's method for nonlinear elasticity.

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s):
        O. Weeger    (2012 - 2015, TU Kaiserslautern),
        A.Shamanskiy (2016 - ...., TU Kaiserslautern)
*/

#pragma once

#include <gsElasticity/gsElasticityNewton2.h>
#include <gsElasticity/gsElasticityAssembler.h>

namespace gismo
{

template <class T>
gsElasticityNewton2<T>::gsElasticityNewton2(gsElasticityAssembler<T> & elasticityAssembler)
    : assembler(elasticityAssembler),
      m_options(defaultOptions())
{

}

template <class T>
gsOptionList gsElasticityNewton2<T>::defaultOptions()
{
    gsOptionList opt;
    /// stopping creteria
    opt.addInt("MaxIter","Maximum number of iterations per loading step",50);
    opt.addReal("AbsTol","Absolute tolerance for the convergence cretiria",1e-12);
    opt.addReal("RelTol","Relative tolerance for the stopping criteria",1e-6);
    /// incremental loading
    opt.addInt("NumIncStep","Number of incremental loading steps",1);
    opt.addInt("NumIterPerStep","Maximum number of Newton's iterations per incremental loading step (not last!); "
                                "always full convergence at the last incremental loading step; "
                                "full convergence at every incremental loading step if < 1",0);
    /// step size control
    opt.addSwitch("Bijective","Check bijectivity; stop the solution process if unable to preserve bijectivity",false);
    opt.addInt("MaxHalving","Maximum number of increment halving in order to preserve bijectivity",0);
    opt.addReal("QualityRatio","Solution quality to preserve w.r.t. the last incremental step"
                               " (0 means to at least ensure bijectivity)",0.);
    opt.addSwitch("DampedNewton","Use damped Newton method for the step size control",false);
    /// additional setting
    opt.addInt("Verbosity","Amount of information printed to the terminal: none, some, all",newtonVerbosity::all);
    opt.addInt("Save","Amount of intermediate soluton information saved: "
                      "only the final solution, "
                      "the first and the last displacement fields at every incremental loading step, "
                      "all",newtonSave::onlyFinal);
    return opt;
}

template <class T>
void gsElasticityNewton2<T>::solve()
{
    index_t numIncSteps = m_options.getInt("NumIncStep");
    index_t maxNumIter = m_options.getInt("MaxIter");
    T absTol = m_options.getReal("AbsTol");
    T relTol = m_options.getReal("RelTol");

    T maxStepSize = 1./numIncSteps;
    T progress = 0.;

    std::vector<gsMatrix<T> > ddof(assembler.patches().dim());
    for (index_t d = 0; d < assembler.patches().dim(); ++d)
        ddof[d] = assembler.fixedDofs(d);

    // incremental loop
    while (abs(progress-1.) > 1e-10)
    {
        T stepSize = (1.-progress >= maxStepSize) ? maxStepSize : 1.-progress;
        if (m_options.getInt("Verbosity") != newtonVerbosity::none)
            gsInfo << "Load: " << progress*100 << "% -> " << (progress+stepSize)*100 << "%\n";
        assembler.options().setReal("ForceScaling",progress+stepSize);
        numIterations = 0;
        converged = false;

        for (index_t d = 0; d < assembler.patches().dim(); ++d)
        {
            gsMatrix<T> tempDDof = stepSize * ddof[d];
            assembler.setFixedDofVector(tempDDof,d);
        }
        computeDisplacementUpdate(true);

        assembler.homogenizeFixedDofs(-1);
        while (!converged && numIterations < maxNumIter)
        {
            computeDisplacementUpdate(false);

            if (residualNorm < absTol || updateNorm < absTol ||
                residualNorm/initResidualNorm < relTol || updateNorm/initUpdateNorm < relTol)
                converged = true;
        }

        progress += stepSize;
    }
}

template <class T>
void gsElasticityNewton2<T>::computeDisplacementUpdate(bool initUpdate)
{
    if (solutions.empty())
        assembler.assemble();
    else
        assembler.assemble(solutions.back());

    gsSparseSolver<>::LU solver(assembler.matrix());
    gsVector<T> solVector = solver.solve(assembler.rhs());

    if (solutions.empty())
    {
        solutions.push_back(gsMultiPatch<T>());
        assembler.constructSolution(solVector,solutions.back());
    }
    else
    {
        if (m_options.getInt("Save") == newtonSave::onlyFinal)
        {
            if (initUpdate)
            {
                gsMultiPatch<T> temp;
                assembler.constructSolution(solVector,temp);
                for (index_t p = 0; p < solutions.back().nPatches(); ++p)
                    solutions.back().patch(p).coefs() += temp.patch(p).coefs();
            }
            else
                assembler.updateSolution(solVector,solutions.back());
        }
        else if (m_options.getInt("Save") == newtonSave::firstAndLastPerIncStep)
        {
            if (numIterations == 0 || numIterations == 1)
            {
                solutions.push_back(gsMultiPatch<T>());
                assembler.constructSolution(solVector,solutions.back());
                for (index_t p = 0; p < solutions.back().nPatches(); ++p)
                    solutions.back().patch(p).coefs() += solutions[solutions.size()-2].patch(p).coefs();
            }
            else
                assembler.updateSolution(solVector,solutions.back());
        }
        else if (m_options.getInt("Save") == newtonSave::all)
        {
            solutions.push_back(gsMultiPatch<T>());
            assembler.constructSolution(solVector,solutions.back());
            for (index_t p = 0; p < solutions.back().nPatches(); ++p)
                solutions.back().patch(p).coefs() += solutions[solutions.size()-2].patch(p).coefs();
        }
    }

    updateNorm = solVector.norm();
    residualNorm = assembler.rhs().norm();

    if (initUpdate)
    {
        initUpdateNorm = updateNorm;
        initResidualNorm = residualNorm;
    }

    printStatus();
    numIterations++;
}

template <class T>
void gsElasticityNewton2<T>::printStatus()
{
    if (m_options.getInt("Verbosity") == newtonVerbosity::all)
        gsInfo << "Iteration: " << numIterations
               << ", resAbs: " << residualNorm
               << ", resRel: " << residualNorm/initResidualNorm
               << ", updAbs: " << updateNorm
               << ", updRel: " << updateNorm/initUpdateNorm << std::endl;
}

template <class T>
void gsElasticityNewton2<T>::plotDeformation(std::string fileName, index_t numSamplingPoints)
{

}

} // namespace ends
