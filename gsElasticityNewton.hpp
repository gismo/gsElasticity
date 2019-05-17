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
gsElasticityNewton<T>::gsElasticityNewton(gsElasticityAssembler<T> & elasticityAssembler,
                                          elasticity_formulation form)
    : assembler(elasticityAssembler),
      formulation(form),
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
    opt.addReal("RelTol","Relative tolerance for the stopping criteria",1e-9);
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
    opt.addInt("Verbosity","Amount of information printed to the terminal: none, some, all",newton_verbosity::all);
    opt.addInt("Save","Amount of intermediate soluton information saved: "
                      "only the final solution, "
                      "the first and the last displacement fields at every incremental loading step, "
                      "all",newton_save::onlyFinal);
    return opt;
}

template <class T>
void gsElasticityNewton<T>::solve()
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

        if (m_options.getInt("Verbosity") != newton_verbosity::none)
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
            if (m_options.getInt("Verbosity") != newton_verbosity::none)
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
                if (m_options.getInt("Verbosity") != newton_verbosity::none)
                    gsInfo << "Interrupted due to bijectivity violation!\n";
                goto abort;
            }

            if (residualNorm < absTol || updateNorm < absTol ||
                residualNorm/initResidualNorm < relTol || updateNorm/initUpdateNorm < relTol)
                converged = true;
        }

        if (m_options.getInt("Verbosity") != newton_verbosity::none)
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
void gsElasticityNewton<T>::solveAdaptive()
{
    GISMO_ENSURE(formulation == elasticity_formulation::displacement,
                 "Adaptive algorithm currently only available for displacement formulation");

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

        if (m_options.getInt("Verbosity") != newton_verbosity::none)
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
        computeUpdate(true);
        if (!bijective)
        {
            if (m_options.getInt("Verbosity") != newton_verbosity::none)
                gsInfo << "Interrupted due to bijectivity violation!\n";
            goto abort;
        }
        // if size of ILS was adaptively halved in an attempt to preserve bijectivity then the ILS is not the last one
        final = (numAdaptHalving == 0 && final && alpha == 1.) ? true : false;
        for (index_t i = 0; i < numAdaptHalving; ++i)
            stepSize *= 0.5;
        stepSize *= alpha;
        // set Dirichlet BC to zero for further Newton's iterations at this ILS
        assembler.homogenizeFixedDofs(-1);
        // if this it not the last ILS, it can be interrupted earlier (to save computation)
        index_t maxNumIter = (!final && m_options.getInt("MaxIterNotLast") > 0) ?
                             m_options.getInt("MaxIterNotLast") : m_options.getInt("MaxIter");
        if (m_options.getInt("Verbosity") != newton_verbosity::none && final)
            gsInfo << "Final loading step!\n";
        // Newton's loop
        while (!converged && numIterations < maxNumIter)
        {
            computeUpdate(false);
            if (!bijective)
            {
                if (m_options.getInt("Verbosity") != newton_verbosity::none)
                    gsInfo << "Interrupted due to bijectivity violation!\n";
                goto abort;
            }

            if (residualNorm < absTol || updateNorm < absTol ||
                residualNorm/initResidualNorm < relTol || updateNorm/initUpdateNorm < relTol)
                converged = true;            
        }

        if (m_options.getInt("Verbosity") != newton_verbosity::none)
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
void gsElasticityNewton<T>::computeUpdate(bool initUpdate)
{
    // reset the counter of adaptive halving
    numAdaptHalving = 0;
    alpha = 1.;

    if (displacements.empty()) // no previous displacement field
        assembler.assemble();
    else if (formulation == elasticity_formulation::displacement)
        assembler.assemble(displacements.back());
    else // mixed dispalcement-pressure formulation
        assembler.assemble(displacements.back(),pressures.back());

    gsSparseSolver<>::SimplicialLDLT solver(assembler.matrix());
    gsVector<T> solVector = solver.solve(assembler.rhs());

    gsMultiPatch<T> incDisplacement, incPressure;
    if (formulation == elasticity_formulation::displacement)
        assembler.constructSolution(solVector,incDisplacement);
    else // mixed dispalcement-pressure formulation
        assembler.constructSolution(solVector,incDisplacement,incPressure);

    if (m_options.getSwitch("BijectivityCheck") || assembler.options().getInt("MaterialLaw") == material_law::neo_hooke_ln)
        bijectivityCheck(incDisplacement);

    if (!bijective && m_options.getInt("MaxHalving") > 0 && formulation == elasticity_formulation::displacement)
        adaptiveHalving(incDisplacement);

    if (m_options.getSwitch("DampedNewton") && formulation == elasticity_formulation::displacement)
        dampedNewton(incDisplacement);

    if (formulation == elasticity_formulation::displacement)
        saveSolution(incDisplacement);
    else // mixed dispalcement-pressure formulation
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
void gsElasticityNewton<T>::printStatus()
{
    if (m_options.getInt("Verbosity") == newton_verbosity::all)
        gsInfo << "Iteration: " << numIterations
               << ", resAbs: " << residualNorm
               << ", resRel: " << residualNorm/initResidualNorm
               << ", updAbs: " << updateNorm
               << ", updRel: " << updateNorm/initUpdateNorm
               << (numAdaptHalving == 0 ? "" : ", adaptHalv: " + std::to_string(numAdaptHalving))
               << (alpha == 1 ? "" : ", damping: " + std::to_string(alpha)) << std::endl;
}

template <class T>
void gsElasticityNewton<T>::plotDeformation(const gsMultiPatch<T> & initDomain,
                                             std::string fileName, index_t numSamplingPoints)
{
    if (m_options.getInt("Verbosity") != newton_verbosity::none)
        gsInfo << "Plotting deformed configurations...\n";

    std::string fileNameOnly = fileName.substr(fileName.find_last_of("/\\")+1);
    gsParaviewCollection collectionMesh(fileName + "_mesh");
    gsParaviewCollection collectionJac(fileName + "_jac");
    index_t res;

    gsMultiPatch<T> configuration;
    for (size_t p = 0; p < initDomain.nPatches(); ++p)
        configuration.addPatch(initDomain.patch(p).clone());

    gsPiecewiseFunction<T> dets;
    for (size_t p = 0; p < configuration.nPatches(); ++p)
        dets.addPiecePointer(new gsDetFunction<T>(configuration,p));

    bool plotJac = true;
    if (numSamplingPoints == 0)
        plotJac = false;

    if (m_options.getInt("Verbosity") == newton_verbosity::all)
        gsInfo << "Step: 0/" << displacements.size() << std::endl;

    gsField<T> detField(configuration,dets,true);
    std::map<std::string,const gsField<T> *> fields;
    fields["Jacobian"] = &detField;
    gsWriteParaviewMultiPhysics(fields,fileName+std::to_string(0),numSamplingPoints,true);

    for (size_t p = 0; p < configuration.nPatches(); ++p)
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


    for (unsigned s = 0; s < displacements.size(); ++s)
    {
        if (m_options.getInt("Verbosity") == newton_verbosity::all)
            gsInfo << "Step: " << s+1 << "/" << displacements.size() << std::endl;

        for (size_t p = 0; p < configuration.nPatches(); ++p)
        {
            configuration.patch(p).coefs() += displacements[s].patch(p).coefs();
            if (s > 0)
               configuration.patch(p).coefs() -= displacements[s-1].patch(p).coefs();
        }

        gsWriteParaviewMultiPhysics(fields,fileName+std::to_string(s+1),numSamplingPoints,true);
        for (size_t p = 0; p < configuration.nPatches(); ++p)
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
    if (displacements.empty())
        displacements.push_back(incDisplacement);
    else
        if (m_options.getInt("Save") == newton_save::onlyFinal)
            for (size_t p = 0; p < displacements.back().nPatches(); ++p)
                displacements.back().patch(p).coefs() += incDisplacement.patch(p).coefs();
        else if (m_options.getInt("Save") == newton_save::firstAndLastPerIncStep)
            if (numIterations == 0 || numIterations == 1)
            {
                displacements.push_back(incDisplacement);
                for (size_t p = 0; p < displacements.back().nPatches(); ++p)
                    displacements.back().patch(p).coefs() += displacements[displacements.size()-2].patch(p).coefs();
            }
            else
                for (size_t p = 0; p < displacements.back().nPatches(); ++p)
                    displacements.back().patch(p).coefs() += incDisplacement.patch(p).coefs();
        else if (m_options.getInt("Save") == newton_save::all)
        {
            displacements.push_back(incDisplacement);
            for (size_t p = 0; p < displacements.back().nPatches(); ++p)
                displacements.back().patch(p).coefs() += displacements[displacements.size()-2].patch(p).coefs();
        }
}

template <class T>
void gsElasticityNewton<T>::saveSolution(const gsMultiPatch<T> & incDisplacement,
                                         const gsMultiPatch<T> & incPressure)
{
    if (displacements.empty())
    {
        displacements.push_back(incDisplacement);
        pressures.push_back(incPressure);
    }
    else
        if (m_options.getInt("Save") == newton_save::onlyFinal)
            for (size_t p = 0; p < displacements.back().nPatches(); ++p)
            {
                displacements.back().patch(p).coefs() += incDisplacement.patch(p).coefs();
                pressures.back().patch(p).coefs() += incPressure.patch(p).coefs();
            }
        else if (m_options.getInt("Save") == newton_save::firstAndLastPerIncStep)
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
        else if (m_options.getInt("Save") == newton_save::all)
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
void gsElasticityNewton<T>::bijectivityCheck(const gsMultiPatch<T> & incDisplacement)
{
    gsMultiPatch<T> tempDisplacement(incDisplacement);
    if (!displacements.empty())
        for (size_t p = 0; p < displacements.back().nPatches(); ++p)
            tempDisplacement.patch(p).coefs() += displacements.back().patch(p).coefs();
    index_t corruptedPatch = assembler.checkSolution(tempDisplacement);
    bijective = corruptedPatch == -1 ? true : false;
}

template <class T>
void gsElasticityNewton<T>::adaptiveHalving(gsMultiPatch<T> & incDisplacement)
{
    T qualityRatio = m_options.getReal("QualityRatio");
    index_t maxAdaptHaling = m_options.getInt("MaxHalving");
    T ratio = -1.;
    T oldRatio = displacements.empty() ? 1. : assembler.solutionJacRatio(displacements.back());
    while (ratio < qualityRatio && numAdaptHalving < maxAdaptHaling)
    {
        for (size_t p = 0; p < incDisplacement.nPatches(); ++p)
            incDisplacement.patch(p).coefs() *= 0.5;
        gsMultiPatch<T> tempDisplacement(incDisplacement);
        if (!displacements.empty())
            for (size_t p = 0; p < displacements.back().nPatches(); ++p)
                tempDisplacement.patch(p).coefs() += displacements.back().patch(p).coefs();

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
    if (!bijective)
    {
        gsInfo << "Can't apply dampling due to bijectivity violation!\n";
        return;
    }

    const index_t maxIter = 20;
    index_t iter = 0;
    const T ratio = 0.1;

    const T g0 = (assembler.rhs().transpose()*assembler.rhs())(0,0); // initial energy value
    gsMatrix<T> G0 = assembler.rhs(); // initial residual
    T alphaA, alphaB, gA, gB, g;
    alphaA = 0.;
    alphaB = 1.;
    gA = g = g0;

    gsMultiPatch<T> tempDisplacement(incDisplacement);
    if (!displacements.empty())
        for (size_t p = 0; p < displacements.back().nPatches(); ++p)
            tempDisplacement.patch(p).coefs() += displacements.back().patch(p).coefs();
    assembler.assemble(tempDisplacement,false); // *false* to assemble only the rhs
    gB = (G0.transpose()*assembler.rhs())(0,0);

    if (gA*gB < 0)
    {
        while (g > ratio*g0 || iter < maxIter)
        {
            alpha = (alphaA*gB-alphaB*gA)/(gB-gA); // symmetric secant
            gsMultiPatch<T> tempDisplacement(incDisplacement);
            if (!displacements.empty())
                for (size_t p = 0; p < displacements.back().nPatches(); ++p)
                {
                    tempDisplacement.patch(p).coefs() *= alpha;
                    tempDisplacement.patch(p).coefs() += displacements.back().patch(p).coefs();
                }
            assembler.assemble(tempDisplacement,false);
            g = (G0.transpose()*assembler.rhs())(0,0);

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
            iter++;
        }
        for (index_t p = 0; p < incDisplacement.nPieces(); ++p)
            incDisplacement.patch(p).coefs() *= alpha;
    }
}

} // namespace ends
