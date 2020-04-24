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
    return opt;
}

template <class T>
bool gsPartitionedFSI<T>::makeTimeStep(T timeStep)
{
    m_nsSolver.saveState();
    m_elSolver.saveState();
    m_aleSolver.saveState();
    gsStopwatch clock;

    numIter = 0;
    nsTime = elTime = aleTime = 0.;
    converged = false;

    gsMultiPatch<> dispOldOld, dispNew,dispDiff, dispGuess;
    omega = 1.;
    absResNorm = 0.;
    initResNorm = 1.;

    while (numIter < m_options.getInt("MaxIter") && !converged)
    {
        // Beam
        clock.restart();
        if (numIter > 0)
            m_elSolver.recoverState();
        m_elSolver.constructSolution(dispDiff);
        m_elSolver.makeTimeStep(timeStep);
        if (numIter == 0)
        {
            m_elSolver.constructSolution(dispOldOld);
            m_elSolver.constructSolution(m_displacement);
            gsMatrix<> vec;
            formVector(m_displacement,vec);
        }
        else if (numIter == 1)
        {
            m_elSolver.constructSolution(m_displacement);
            m_elSolver.constructSolution(dispGuess);
            gsMatrix<> vecA, vecB;
            formVector(dispOldOld,vecA);
            formVector(m_displacement,vecB);
            initResNorm = (vecB-vecA).norm();
        }
        else
        {
            m_elSolver.constructSolution(dispNew);
            aitken(dispOldOld,m_displacement,dispGuess,dispNew);
        }
        for (index_t p = 0; p < dispDiff.nPatches(); ++p)
            dispDiff.patch(p).coefs() = m_displacement.patch(p).coefs() - dispDiff.patch(p).coefs();
        elTime += clock.stop();

        // ALE
        clock.restart();
        // undo last ALE
        for (index_t p = 0; p < m_nsSolver.alePatches().uniquePatches.size(); ++p)
        {
            index_t pFlow = m_nsSolver.alePatches().uniquePatches[p].second;
            index_t pALE = m_nsSolver.alePatches().uniquePatches[p].first;
            m_nsSolver.assembler().patches().patch(pFlow).coefs() -= m_ALEdisplacment.patch(pALE).coefs();
            m_nsSolver.mAssembler().patches().patch(pFlow).coefs() -= m_ALEdisplacment.patch(pALE).coefs();
        }
        // recover ALE at the start of timestep
        if (numIter > 0)
            m_aleSolver.recoverState();
        m_aleSolver.constructSolution(m_ALEdisplacment);
        // update
        index_t badPatch = m_aleSolver.updateMesh();
        // construct ALE velocity
        m_aleSolver.constructSolution(m_ALEvelocity);
        for (index_t p = 0; p < m_ALEvelocity.nPatches(); ++p)
        {
            m_ALEvelocity.patch(p).coefs() -= m_ALEdisplacment.patch(p).coefs();
            m_ALEvelocity.patch(p).coefs() /= timeStep;
        }
        // construct ALE displacement
        m_aleSolver.constructSolution(m_ALEdisplacment);
        // apply new ALE to the flow domain
        for (index_t p = 0; p < m_nsSolver.alePatches().uniquePatches.size(); ++p)
        {
            index_t pFlow = m_nsSolver.alePatches().uniquePatches[p].second;
            index_t pALE = m_nsSolver.alePatches().uniquePatches[p].first;
            m_nsSolver.assembler().patches().patch(pFlow).coefs() += m_ALEdisplacment.patch(pALE).coefs();
            m_nsSolver.mAssembler().patches().patch(pFlow).coefs() += m_ALEdisplacment.patch(pALE).coefs();
        }
        aleTime += clock.stop();

        // test if any patch is not bijective
        if (badPatch != -1)
            return false;

        // FLOW
        clock.restart();
        for (index_t p = 0; p < m_nsSolver.alePatches().sidesA.size(); ++p)
            m_nsSolver.assembler().setFixedDofs(m_nsSolver.alePatches().sidesB[p].patch,m_nsSolver.alePatches().sidesB[p].side(),
                                                m_ALEvelocity.patch(m_nsSolver.alePatches().sidesA[p].patch).boundary(
                                                    m_nsSolver.alePatches().sidesA[p].side())->coefs());
        if (numIter > 0)
            m_nsSolver.recoverState();
        m_nsSolver.makeTimeStep(timeStep);
        m_nsSolver.constructSolution(m_velocity,m_pressure);
        nsTime += clock.stop();

        ++numIter;
    }
    return true;
}

template <class T>
void gsPartitionedFSI<T>::formVector(const gsMultiPatch<T> & disp, gsMatrix<T> & vector)
{
    index_t dim = disp.patch(0).parDim();

    index_t totalSize = 0;
    for (index_t i = 0; i < m_aleSolver.interfaceStr2Mesh().sidesA.size(); ++i)
    {
        index_t patch = m_aleSolver.interfaceStr2Mesh().sidesA[i].patch;
        boxSide side = m_aleSolver.interfaceStr2Mesh().sidesA[i].side();
        totalSize += disp.patch(patch).boundary(side)->coefs().rows();
    }

    vector.setZero(totalSize*dim,1);
    index_t filledSize = 0;
    for (index_t i = 0; i < m_aleSolver.interfaceStr2Mesh().sidesA.size(); ++i)
    {
        index_t patch = m_aleSolver.interfaceStr2Mesh().sidesA[i].patch;
        boxSide side = m_aleSolver.interfaceStr2Mesh().sidesA[i].side();
        index_t size = disp.patch(patch).boundary(side)->coefs().rows();
        for (index_t d = 0; d < dim;++d)
        {
            vector.middleRows(filledSize,size) = disp.patch(patch).boundary(side)->coefs().col(d);
            filledSize += size;
        }
    }
}

template <class T>
void gsPartitionedFSI<T>::aitken(gsMultiPatch<T> & dispA, gsMultiPatch<T> & dispB,
                                 gsMultiPatch<T> & dispB2, gsMultiPatch<T> & dispC)
{
    gsMatrix<> vecA,vecB,vecB2,vecC;
    formVector(dispA,vecA);
    formVector(dispB,vecB);
    formVector(dispB2,vecB2);
    formVector(dispC,vecC);

    omega = -1*omega*((vecB2-vecA).transpose()*(vecC-vecB -vecB2+vecA))(0,0)/((vecC-vecB -vecB2+vecA).transpose()*(vecC-vecB -vecB2+vecA))(0,0);
    for (index_t p = 0; p < dispA.nPatches(); ++p)
    {
        dispA.patch(p).coefs() = dispB.patch(p).coefs();
        dispB.patch(p).coefs() += omega*(dispC.patch(p).coefs()-dispB.patch(p).coefs());
        dispB2.patch(p).coefs() = dispC.patch(p).coefs();
    }
    absResNorm = ((vecC-vecB)*omega).norm();
    gsInfo << numIter << " " << absResNorm << std::endl;
    if (absResNorm < m_options.getReal("AbsTol") || absResNorm/initResNorm < m_options.getReal("RelTol"))
        converged = true;
}

} // namespace ends
