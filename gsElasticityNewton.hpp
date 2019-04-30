/** @file gsElasticityNewton.hpp

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

#include <gsElasticity/gsElasticityNewton.h>
#include <gsElasticity/gsElasticityAssembler.h>
#include <gsCore/gsField.h>

namespace gismo
{

template <class T>
gsElasticityNewton<T>::gsElasticityNewton(gsElasticityAssembler<T> & elasticityAssembler)
    : assembler(elasticityAssembler),
      m_options(defaultOptions())
{

}

template <class T>
gsOptionList gsElasticityNewton<T>::defaultOptions()
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
    /// step size control
    opt.addSwitch("BijectivityCheck","Check bijectivity; stop the solution process if unable to preserve bijectivity",false);
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
void gsElasticityNewton<T>::solve()
{
    T absTol = m_options.getReal("AbsTol");
    T relTol = m_options.getReal("RelTol");
    bijective = true;
    // preferred number of ILSs; may increase to preserve bijectivity
    index_t numIncSteps = m_options.getInt("NumIncStep");
    T maxStepSize = 1./numIncSteps;
    // computed portion of total loading
    T progress = 0.;
    // save computed by the assembler Dirichet DoFs for further use
    std::vector<gsMatrix<T> > ddof(assembler.patches().dim());
    for (index_t d = 0; d < assembler.patches().dim(); ++d)
        ddof[d] = assembler.fixedDofs(d);

    // incremental loop
    while (abs(progress-1.) > 1e-10 && bijective)
    {
        // reset the status of Newton's method
        numIterations = 0;
        converged = false;
        // determine whether this is the last ILS
        T stepSize = 1.-progress > maxStepSize ? maxStepSize : 1. - progress;
        bool final = abs(1. - progress - stepSize) < 1e-10;

        if (m_options.getInt("Verbosity") != newtonVerbosity::none)
            gsInfo << "Load: " << progress*100 << "% -> " << (progress+stepSize)*100 << "%\n";
        // set load scaling for RHS and Neumann BC
        assembler.options().setReal("ForceScaling",progress+stepSize);
        // set load scaling for Dirichlet BC;
        // a temporary variable is necessary because of the memory swap
        for (index_t d = 0; d < assembler.patches().dim(); ++d)
        {
            gsMatrix<T> tempDDof = stepSize * ddof[d];
            assembler.setFixedDofVector(tempDDof,d);
        }
        computeDisplacementUpdate(true);
        if (!bijective)
        {
            if (m_options.getInt("Verbosity") != newtonVerbosity::none)
                gsInfo << "Interrupted due to bijectivity violation!\n";
            goto abort;
        }
        // if size of ILS was adaptively halved in an attempt to preserve bijectivity then the ILS is not the last one
        final = (numAdaptHalving == 0 && final) ? true : false;
        for (index_t i = 0; i < numAdaptHalving; ++i)
            stepSize *= 0.5;

        stepSize *= alpha;
        // set Dirichlet BC to zero for further Newton's iterations at this ILS
        assembler.homogenizeFixedDofs(-1);
        // if this it not the last ILS, it can be interrupted earlier (to save computation)
        index_t maxNumIter = (!final && m_options.getInt("MaxIterNotLast") > 0) ?
                             m_options.getInt("MaxIterNotLast") : m_options.getInt("MaxIter");
        if (m_options.getInt("Verbosity") != newtonVerbosity::none && final)
            gsInfo << "Final loading step!\n";
        // Newton's loop
        while (!converged && numIterations < maxNumIter)
        {
            computeDisplacementUpdate(false);
            if (!bijective)
            {
                if (m_options.getInt("Verbosity") != newtonVerbosity::none)
                    gsInfo << "Interrupted due to bijectivity violation!\n";
                goto abort;
            }

            if (residualNorm < absTol || updateNorm < absTol ||
                residualNorm/initResidualNorm < relTol || updateNorm/initUpdateNorm < relTol)
                converged = true;            
        }

        if (m_options.getInt("Verbosity") != newtonVerbosity::none)
        {
            if (converged)
                gsInfo << "Newton's method converged after " << numIterations << " iterations\n";
            else
                gsInfo << "Newton's method interrupted after " << maxNumIter << " iterations\n";
        }
        // advance the computed portion of total loading
        progress += stepSize;
    }
    abort:;
}

template <class T>
void gsElasticityNewton<T>::computeDisplacementUpdate(bool initUpdate)
{
    // reset the counter of adaptive halving
    numAdaptHalving = 0;
    alpha = 1.;

    if (solutions.empty()) // no previous displacement field
        assembler.assemble();
    else // use previous displacement field to assemble the problem
        assembler.assemble(solutions.back());

    gsSparseSolver<>::LU solver(assembler.matrix());
    gsVector<T> solVector = solver.solve(assembler.rhs());
    gsMultiPatch<T> incDisplacement;
    assembler.constructSolution(solVector,incDisplacement);


    if (m_options.getSwitch("BijectivityCheck") || assembler.options().getInt("MaterialLaw") == material_law::neo_hooke_ln)
        bijectivityCheck(incDisplacement);

    if (!bijective && m_options.getInt("MaxHalving") > 0)
        adaptiveHalving(incDisplacement);

    if (m_options.getSwitch("DampedNewton"))
        dampedNewton(incDisplacement);

    saveSolution(incDisplacement);

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
void gsElasticityNewton<T>::printStatus()
{
    if (m_options.getInt("Verbosity") == newtonVerbosity::all)
        gsInfo << "Iteration: " << numIterations
               << ", resAbs: " << residualNorm
               << ", resRel: " << residualNorm/initResidualNorm
               << ", updAbs: " << updateNorm
               << ", updRel: " << updateNorm/initUpdateNorm
               << (numAdaptHalving == 0 ? "" : ", adaptHalv: " + std::to_string(numAdaptHalving)) <<  std::endl;
}

template <class T>
void gsElasticityNewton<T>::plotDeformation(const gsMultiPatch<T> & initDomain,
                                             std::string fileName, index_t numSamplingPoints)
{
    if (m_options.getInt("Verbosity") != newtonVerbosity::none)
        gsInfo << "Plotting deformed configurations...\n";

    std::string fileNameOnly = fileName.substr(fileName.find_last_of("/\\")+1);
    gsParaviewCollection collectionMesh(fileName + "_mesh");
    gsParaviewCollection collectionJac(fileName + "_jac");
    index_t res;

    gsMultiPatch<T> configuration;
    for (index_t p = 0; p < initDomain.nPatches(); ++p)
        configuration.addPatch(initDomain.patch(p).clone());

    gsPiecewiseFunction<T> dets;
    for (index_t p = 0; p < configuration.nPatches(); ++p)
        dets.addPiecePointer(new gsDetFunction<T>(configuration,p));

    bool plotJac = true;
    if (numSamplingPoints == 0)
        plotJac = false;

    if (m_options.getInt("Verbosity") == newtonVerbosity::all)
        gsInfo << "Step: 0/" << solutions.size() << std::endl;

    gsField<T> detField(configuration,dets,true);
    std::map<std::string,const gsField<T> *> fields;
    fields["Jacobian"] = &detField;
    gsWriteParaviewMultiPhysics(fields,fileName+std::to_string(0),numSamplingPoints,true);

    for (index_t p = 0; p < configuration.nPatches(); ++p)
    {
        collectionMesh.addTimestep(fileNameOnly + std::to_string(0),p,0,"_mesh.vtp");
        if (plotJac)
            collectionJac.addTimestep(fileNameOnly + std::to_string(0),p,0,".vts");
        else
            res = system(("rm " + fileName + std::to_string(0) + std::to_string(p) + ".vts").c_str());
        GISMO_ENSURE(res == 0, "Problems with deleting files\n");
    }
    res = system(("rm " + fileName + std::to_string(0) + ".pvd").c_str());
    GISMO_ENSURE(res == 0, "Problems with deleting files\n");


    for (unsigned s = 0; s < solutions.size(); ++s)
    {
        if (m_options.getInt("Verbosity") == newtonVerbosity::all)
            gsInfo << "Step: " << s+1 << "/" << solutions.size() << std::endl;

        for (index_t p = 0; p < configuration.nPatches(); ++p)
        {
            configuration.patch(p).coefs() += solutions[s].patch(p).coefs();
            if (s > 0)
               configuration.patch(p).coefs() -= solutions[s-1].patch(p).coefs();
        }

        gsWriteParaviewMultiPhysics(fields,fileName+std::to_string(s+1),numSamplingPoints,true);
        for (index_t p = 0; p < configuration.nPatches(); ++p)
        {
            collectionMesh.addTimestep(fileNameOnly + std::to_string(s+1),p,s+1,"_mesh.vtp");

            if (plotJac)
                collectionJac.addTimestep(fileNameOnly + std::to_string(s+1),p,s+1,".vts");
            else
                res = system(("rm " + fileName + std::to_string(s+1) + std::to_string(p) + ".vts").c_str());

            GISMO_ENSURE(res == 0, "Problems with deleting files\n");
        }
        res = system(("rm " + fileName + std::to_string(s+1) + ".pvd").c_str());
        GISMO_ENSURE(res == 0, "Problems with deleting files\n");
    }

    collectionMesh.save();
    if (plotJac)
        collectionJac.save();

    (void)res;
}

template <class T>
void gsElasticityNewton<T>::saveSolution(const gsMultiPatch<T> & incDisplacement)
{
    if (solutions.empty())
        solutions.push_back(incDisplacement);
    else
        if (m_options.getInt("Save") == newtonSave::onlyFinal)
            for (index_t p = 0; p < solutions.back().nPatches(); ++p)
                solutions.back().patch(p).coefs() += incDisplacement.patch(p).coefs();
        else if (m_options.getInt("Save") == newtonSave::firstAndLastPerIncStep)
            if (numIterations == 0 || numIterations == 1)
            {
                solutions.push_back(incDisplacement);
                for (index_t p = 0; p < solutions.back().nPatches(); ++p)
                    solutions.back().patch(p).coefs() += solutions[solutions.size()-2].patch(p).coefs();
            }
            else
                for (index_t p = 0; p < solutions.back().nPatches(); ++p)
                    solutions.back().patch(p).coefs() += incDisplacement.patch(p).coefs();
        else if (m_options.getInt("Save") == newtonSave::all)
        {
            solutions.push_back(incDisplacement);
            for (index_t p = 0; p < solutions.back().nPatches(); ++p)
                solutions.back().patch(p).coefs() += solutions[solutions.size()-2].patch(p).coefs();
        }
}

template <class T>
void gsElasticityNewton<T>::bijectivityCheck(const gsMultiPatch<T> & incDisplacement)
{
    gsMultiPatch<T> tempDisplacement(incDisplacement);
    if (!solutions.empty())
        for (index_t p = 0; p < solutions.back().nPatches(); ++p)
            tempDisplacement.patch(p).coefs() += solutions.back().patch(p).coefs();
    index_t corruptedPatch = assembler.checkSolution(tempDisplacement);
    if (corruptedPatch == -1)
        bijective = true;
    else
        bijective = false;
}

template <class T>
void gsElasticityNewton<T>::adaptiveHalving(gsMultiPatch<T> & incDisplacement)
{
    T qualityRatio = m_options.getReal("QualityRatio");
    index_t maxAdaptHaling = m_options.getInt("MaxHalving");
    T ratio = -1.;
    T oldRatio = solutions.empty() ? 1. : assembler.solutionJacRatio(solutions.back());
    while (ratio < qualityRatio && numAdaptHalving < maxAdaptHaling)
    {
        for (index_t p = 0; p < incDisplacement.nPatches(); ++p)
            incDisplacement.patch(p).coefs() *= 0.5;
        gsMultiPatch<T> tempDisplacement(incDisplacement);
        if (!solutions.empty())
            for (index_t p = 0; p < solutions.back().nPatches(); ++p)
                tempDisplacement.patch(p).coefs() += solutions.back().patch(p).coefs();

        ratio = assembler.solutionJacRatio(tempDisplacement) / oldRatio;
        numAdaptHalving++;
    }

    if (ratio >= qualityRatio)
        bijective = true;
    else
        bijectivityCheck(incDisplacement);
}

template <class T>
void gsElasticityNewton<T>::dampedNewton(gsMultiPatch<T> & incDisplacement)
{
    gsInfo << "Damped Newton...\n";
    const index_t maxIter = 100;
    const T ratio = 0.8;

    const T g0 = (assembler.rhs().transpose()*assembler.rhs())(0,0);
    gsMatrix<T> G0 = assembler.rhs();
    gsInfo << "g0 " << g0 << std::endl;
    T alphaA, alphaB, gA, gB;
    alphaA = 0.;
    alphaB = 1.;
    gA = g0;

    gsMultiPatch<T> tempDisplacement(incDisplacement);
    if (!solutions.empty())
        for (index_t p = 0; p < solutions.back().nPatches(); ++p)
            tempDisplacement.patch(p).coefs() += solutions.back().patch(p).coefs();
    assembler.assemble(tempDisplacement,false);
    gB = (G0.transpose()*assembler.rhs())(0,0);
    gsInfo << "gB " << gB << std::endl;

    if (gA*gB < 0)
    {
        gsInfo << "GO...\n";
        for (int i = 0; i < maxIter; ++i)
        {
            alpha = (alphaA*gB-alphaB*gA)/(gB-gA);

            gsMultiPatch<T> tempDisplacement(incDisplacement);
            if (!solutions.empty())
                for (index_t p = 0; p < solutions.back().nPatches(); ++p)
                {
                    tempDisplacement.patch(p).coefs() *= alpha;
                    tempDisplacement.patch(p).coefs() += solutions.back().patch(p).coefs();
                }
            assembler.assemble(tempDisplacement,false);
            T g = (G0.transpose()*assembler.rhs())(0,0);
            gsInfo << "alpha " << alpha << " g " << g << std::endl;

                if (gA*g < 0)
                {
                    alphaB = alpha;
                    gB = g;
                }
                else
                {
                    alphaA = alpha;
                    gA = g;
                }
        }
        for (index_t p = 0; p < incDisplacement.nPieces(); ++p)
            incDisplacement.patch(p).coefs() *= alpha;
    }
}

} // namespace ends
