/** @file gsPartitionedFSI.hpp

    @brief Implementation of gsPartitionedFSI.

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s):
        A.Shamanskiy (2016 - ...., TU Kaiserslautern)
*/

#pragma once

#include <gsElasticity/gsPartitionedFSI.h>

#include <gsElasticity/gsNsTimeIntegrator.h>
#include <gsElasticity/gsElTimeIntegrator.h>
#include <gsElasticity/gsALE.h>
#include <gsUtils/gsStopwatch.h>
#include <gsElasticity/gsGeoUtils.h>

namespace gismo
{

template <class T>
gsPartitionedFSI<T>::gsPartitionedFSI(gsNsTimeIntegrator<T> & nsSolver,
                                      gsMultiPatch<T> & velocity, gsMultiPatch<T> & pressure,
                                      gsElTimeIntegrator<T> & elSolver,
                                      gsMultiPatch<T> & displacement,
                                      gsALE<T> & aleSolver,
                                      gsMultiPatch<T> & aleDisplacement, gsMultiPatch<T> & aleVelocity) :
    m_nsSolver(nsSolver),
    m_velocity(velocity),
    m_pressure(pressure),
    m_elSolver(elSolver),
    m_displacement(displacement),
    m_aleSolver(aleSolver),
    m_ALEdisplacment(aleDisplacement),
    m_ALEvelocity(aleVelocity),
    m_options(defaultOptions())
{

}

template <class T>
gsOptionList gsPartitionedFSI<T>::defaultOptions()
{
    gsOptionList opt;
    opt.addInt("MaxIter","Maximum number of coupling iterations per time step",10);
    opt.addReal("AbsTol","Absolute tolerance for the convergence creterion",1e-10);
    opt.addReal("RelTol","Absolute tolerance for the convergence creterion",1e-6);
    opt.addInt("Verbosity","Amount of information printed to the terminal: none, some, all",solver_verbosity::none);
    return opt;
}

template <class T>
bool gsPartitionedFSI<T>::makeTimeStep(T timeStep)
{
    // save states of the component solvers at the beginning of the time step
    m_nsSolver.saveState();
    m_elSolver.saveState();
    m_aleSolver.saveState();

    // reset the solver
    numIter = 0;
    converged = false;
    omega = 1.;

    // reset time profiling info
    gsStopwatch clock;
    nsTime = elTime = aleTime = 0.;

    gsMultiPatch<> dispOldOld, dispOld, dispOldGuess;

    while (numIter < m_options.getInt("MaxIter") && !converged)
    {
        // ================== Structure section ================ //
        clock.restart();
        if (numIter > 0) // recover the solver state from the time step beginning
            m_elSolver.recoverState();

        m_elSolver.makeTimeStep(timeStep);

        if (numIter == 0) // save displacement i-2, no correction
        {
            m_elSolver.constructSolution(dispOldOld);
            m_elSolver.constructSolution(m_displacement);
        }
        else if (numIter == 1) // save displacement i-1 as a guess and a corrected solution
        {
            m_elSolver.constructSolution(dispOld);
            m_elSolver.constructSolution(dispOldGuess);
            m_elSolver.constructSolution(m_displacement);
            gsMatrix<> vecA, vecB;
            formVector(dispOldOld,vecA);
            formVector(m_displacement,vecB);
            absResNorm = initResNorm = (vecB-vecA).norm()/sqrt(vecB.rows());
        }
        else // save displacement as a current guess i and apply Aitken relaxation
        {
            m_elSolver.constructSolution(m_displacement);
            aitken(dispOldOld,dispOldGuess,dispOld,m_displacement);
        }

        if (numIter > 0 && m_options.getInt("Verbosity") == solver_verbosity::all)
            gsInfo << numIter << ": absRes " << absResNorm << ", relRes " << absResNorm/initResNorm << std::endl;


        elTime += clock.stop();
        // =================================================================== //


        // ============= Flow mesh/ALE section ===================== //
        clock.restart();
        // recover ALE at the start of timestep
        if (numIter > 0)
            m_aleSolver.recoverState();

        // undo last ALE deformation of the flow domain
        for (index_t p = 0; p < m_nsSolver.aleInterface().patches.size(); ++p)
        {
            index_t pFlow = m_nsSolver.aleInterface().patches[p].second;
            index_t pALE = m_nsSolver.aleInterface().patches[p].first;
            m_nsSolver.assembler().patches().patch(pFlow).coefs() -= m_ALEdisplacment.patch(pALE).coefs();
            m_nsSolver.mAssembler().patches().patch(pFlow).coefs() -= m_ALEdisplacment.patch(pALE).coefs();
        }

        // save ALE displacement at the beginning of the time step for ALE velocity computation
        m_aleSolver.constructSolution(m_ALEvelocity);
        // update ALE
        if (m_aleSolver.updateMesh() != -1)
            return false; // if the new ALE deformation is not bijective, stop the simulation
        // construct new ALE displacement
        m_aleSolver.constructSolution(m_ALEdisplacment);
        for (index_t p = 0; p < m_ALEvelocity.nPatches(); ++p)
            m_ALEvelocity.patch(p).coefs() = (m_ALEdisplacment.patch(p).coefs() - m_ALEvelocity.patch(p).coefs()) / timeStep;

        // apply new ALE deformation to the flow domain
        for (index_t p = 0; p < m_nsSolver.aleInterface().patches.size(); ++p)
        {
            index_t pFlow = m_nsSolver.aleInterface().patches[p].second;
            index_t pALE = m_nsSolver.aleInterface().patches[p].first;
            m_nsSolver.assembler().patches().patch(pFlow).coefs() += m_ALEdisplacment.patch(pALE).coefs();
            m_nsSolver.mAssembler().patches().patch(pFlow).coefs() += m_ALEdisplacment.patch(pALE).coefs();
        }

        aleTime += clock.stop();
        // =================================================================== //


        // ======================= Flow section ============================== //
        clock.restart();
        if (numIter > 0) // recover the solver state from the time step beginning
            m_nsSolver.recoverState();

        // set velocity boundary condition on the FSI interface; velocity comes from the ALE velocity;
        // FSI inteface info is contained in the Navier-Stokes solver
        for (index_t p = 0; p < m_nsSolver.aleInterface().sidesA.size(); ++p)
        {
            index_t pFlow = m_nsSolver.aleInterface().sidesB[p].patch;
            boxSide sFlow = m_nsSolver.aleInterface().sidesB[p].side();
            index_t pALE = m_nsSolver.aleInterface().sidesA[p].patch;
            boxSide sALE = m_nsSolver.aleInterface().sidesA[p].side();
            m_nsSolver.assembler().setFixedDofs(pFlow,sFlow,m_ALEvelocity.patch(pALE).boundary(sALE)->coefs());
        }

        m_nsSolver.makeTimeStep(timeStep);
        m_nsSolver.constructSolution(m_velocity,m_pressure);

        nsTime += clock.stop();
        // =================================================================== //


        ++numIter;
    }

    if (m_options.getInt("Verbosity") != solver_verbosity::none && numIter > 1)
    {
        if (converged)
            gsInfo << "Converged after " << numIter << " iters, absRes "
                   << absResNorm << ", relRes " << absResNorm/initResNorm << std::endl;
        else
            gsInfo << "Terminated after " << numIter << " iters, absRes "
                   << absResNorm << ", relRes " << absResNorm/initResNorm << std::endl;
    }

    return true;
}

template <class T>
void gsPartitionedFSI<T>::formVector(const gsMultiPatch<T> & disp, gsMatrix<T> & vector)
{
    index_t dim = disp.patch(0).parDim();

    index_t totalSize = 0;
    for (index_t i = 0; i < m_aleSolver.interface().sidesA.size(); ++i)
    {
        index_t patch = m_aleSolver.interface().sidesA[i].patch;
        boxSide side = m_aleSolver.interface().sidesA[i].side();
        totalSize += disp.patch(patch).boundary(side)->coefs().rows();
    }

    vector.setZero(totalSize*dim,1);
    index_t filledSize = 0;
    for (index_t i = 0; i < m_aleSolver.interface().sidesA.size(); ++i)
    {
        index_t patch = m_aleSolver.interface().sidesA[i].patch;
        boxSide side = m_aleSolver.interface().sidesA[i].side();
        index_t size = disp.patch(patch).boundary(side)->coefs().rows();
        for (index_t d = 0; d < dim;++d)
        {
            vector.middleRows(filledSize,size) = disp.patch(patch).boundary(side)->coefs().col(d);
            filledSize += size;
        }
    }
}

template <class T>
void gsPartitionedFSI<T>::aitken(gsMultiPatch<T> & dispOO, gsMultiPatch<T> & dispOG,
                                 gsMultiPatch<T> & dispO, gsMultiPatch<T> & dispN)
{
    gsMatrix<> vecOO,vecOG,vecO,vecN;
    formVector(dispOO,vecOO);
    formVector(dispOG,vecOG);
    formVector(dispO,vecO);
    formVector(dispN,vecN);

    gsMatrix<> vecTemp = vecN - vecO - vecOG + vecOO;
    omega = -1*omega * ((vecOG - vecOO).transpose()*vecTemp)(0,0) /
            (vecTemp.transpose()*vecTemp)(0,0);

    for (index_t p = 0; p < dispOO.nPatches(); ++p)
    {
        dispOO.patch(p).coefs() = dispO.patch(p).coefs();
        dispOG.patch(p).coefs() = dispN.patch(p).coefs();
        dispN.patch(p).coefs() = omega * dispN.patch(p).coefs() + (1-omega) * dispO.patch(p).coefs();
        dispO.patch(p).coefs() = dispN.patch(p).coefs();
    }

    formVector(dispN,vecN);
    absResNorm = ((vecN-vecO)*omega).norm()/sqrt(vecN.rows());
    if (absResNorm < m_options.getReal("AbsTol") || absResNorm/initResNorm < m_options.getReal("RelTol"))
        converged = true;
}

} // namespace ends
