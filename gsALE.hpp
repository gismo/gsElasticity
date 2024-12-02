/** @file gsALE.hpp

    @brief Implementation of gsALE.

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s):
        A.Shamanskiy (2016 - ...., TU Kaiserslautern)
*/

#pragma once

#include <gsElasticity/gsALE.h>

#include <gsElasticity/gsElasticityAssembler.h>
#include <gsElasticity/gsElPoissonAssembler.h>
#include <gsElasticity/gsBiharmonicAssembler.h>
#include <gsElasticity/gsIterative.h>
#include <gsElasticity/gsGeoUtils.h>
#include <gsCore/gsConstantFunction.h>

namespace gismo
{

template <class T>
gsALE<T>::gsALE(gsMultiPatch<T> & geometry, const gsMultiPatch<T> & displacement,
                const gsBoundaryInterface & interfaceS2M, ale_method::method method)
    : disp(displacement),
      m_interface(interfaceS2M),
      methodALE(method),
      m_options(defaultOptions()),
      initialized(false),
      hasSavedState(false)
{
    // create input for the assembler
    gsMultiBasis<T> basis(geometry);
    gsBoundaryConditions<T> bcInfo;
    for (gsMultiPatch<>::const_biterator it = geometry.bBegin(); it != geometry.bEnd(); ++it)
        for (index_t d = 0; d < geometry.parDim(); ++d)
            bcInfo.addCondition(it->patch,it->side(),condition_type::dirichlet,0,d);
    gsConstantFunction<T> rhs(gsVector<T>::Zero(geometry.parDim()),geometry.parDim());

    // define assembler according to the method
    if (methodALE == ale_method::TINE || methodALE == ale_method::TINE_StVK || methodALE == ale_method::ILE || methodALE == ale_method::LE)
    {
        assembler = typename gsBaseAssembler<T>::uPtr(new gsElasticityAssembler<T>(geometry,basis,bcInfo,rhs));
        if (methodALE == ale_method::TINE || methodALE == ale_method::TINE_StVK)
        {
            assembler->options().setSwitch("Check",false);
            solverNL = typename gsIterative<T>::uPtr(new gsIterative<T>(*assembler,
                                                                        gsMatrix<>::Zero(assembler->numDofs(),1),
                                                                        assembler->allFixedDofs()));
            solverNL->options().setInt("Verbosity",solver_verbosity::none);
            solverNL->options().setInt("Solver",linear_solver::LDLT);
        }
        if (methodALE == ale_method::TINE)
            assembler->options().setInt("MaterialLaw",material_law::neo_hooke_ln);
        if (methodALE == ale_method::TINE_StVK)
            assembler->options().setInt("MaterialLaw",material_law::saint_venant_kirchhoff);

        assembler->constructSolution(gsMatrix<T>::Zero(assembler->numDofs(),1),assembler->allFixedDofs(),ALEdisp);
    }
    else if (methodALE == ale_method::HE || methodALE == ale_method::IHE)
    {
        assembler = typename gsBaseAssembler<T>::uPtr(new gsElPoissonAssembler<T>(geometry,basis,bcInfo,rhs));
        assembler->constructSolution(gsMatrix<T>::Zero(assembler->numDofs(),geometry.parDim()),assembler->allFixedDofs(),ALEdisp);
    }
    else if (methodALE == ale_method::BHE || methodALE == ale_method::IBHE)
    {
        gsBoundaryConditions<T> bcInfoBHE;
        for (gsMultiPatch<>::const_biterator it = geometry.bBegin(); it != geometry.bEnd(); ++it)
            bcInfoBHE.addCondition(it->patch,it->side(),condition_type::dirichlet,0);
        assembler = typename gsBaseAssembler<T>::uPtr(new gsBiharmonicAssembler<T>(geometry,basis,bcInfoBHE,rhs));
        assembler->constructSolution(gsMatrix<T>::Zero(assembler->numDofs(),geometry.parDim()),assembler->allFixedDofs(),ALEdisp);
    }
    else
        for (size_t p = 0; p < geometry.nPatches(); ++p)
        {
            ALEdisp.addPatch(geometry.patch(p).clone());
            ALEdisp.patch(p).coefs() *= 0.;
        }
}

template <class T>
gsOptionList gsALE<T>::defaultOptions()
{
    gsOptionList opt;
    opt.addReal("PoissonsRatio","Poisson's ratio of the material (only for elasticity-based methods)",0.4);
    opt.addReal("LocalStiff","Stiffening degree for the Jacobian-based local stiffening",0.);
    opt.addSwitch("Check","Check bijectivity of the resulting ALE displacement field",true);
    opt.addInt("NumIter","Number of iterations for nonlinear methods",1);
    return opt;
}


template <class T>
void gsALE<T>::constructSolution(gsMultiPatch<T> & solution) const
{
    solution.clear();
    for (size_t p = 0; p < ALEdisp.nPatches(); ++p)
        solution.addPatch(ALEdisp.patch(p).clone());
}

template <class T>
void gsALE<T>::initialize()
{
    assembler->options().setReal("LocalStiff",m_options.getReal("LocalStiff"));
    if (methodALE == ale_method::LE || methodALE == ale_method::ILE || methodALE == ale_method::TINE || methodALE == ale_method::TINE_StVK)
        assembler->options().setReal("PoissonsRatio",m_options.getReal("PoissonsRatio"));
    if (methodALE == ale_method::LE || methodALE == ale_method::HE || methodALE == ale_method::BHE)
        assembler->assemble(true);
    if (methodALE == ale_method::TINE || methodALE == ale_method::TINE_StVK)
        solverNL->options().setInt("MaxIters",m_options.getInt("NumIter"));

    initialized = true;
}

template <class T>
index_t gsALE<T>::updateMesh()
{
    if (!initialized)
        initialize();

    switch (methodALE)
    {
    case ale_method::HE: return linearMethod(); break;
    case ale_method::IHE: return linearIncrementalMethod(); break;
    case ale_method::LE: return linearMethod(); break;
    case ale_method::ILE: return linearIncrementalMethod(); break;
    case ale_method::TINE: return nonlinearMethod(); break;
    case ale_method::TINE_StVK: return nonlinearMethod(); break;
    case ale_method::BHE: return linearMethod(); break;
    case ale_method::IBHE: return linearIncrementalMethod(); break;
    default: return -1;
    }
}

template <class T>
index_t gsALE<T>::linearMethod()
{

    for (size_t i = 0; i < m_interface.sidesA.size(); ++i)
        assembler->setFixedDofs(m_interface.sidesB[i].patch,
                                m_interface.sidesB[i].side(),
                                disp.patch(m_interface.sidesA[i].patch).boundary(m_interface.sidesA[i].side())->coefs(),
                                methodALE == ale_method::LE ? false : true);
    assembler->eliminateFixedDofs();

#ifdef GISMO_WITH_PARDISO
    gsSparseSolver<>::PardisoLDLT solver(assembler->matrix());
    gsMatrix<> solVector = solver.solve(assembler->rhs());
#else
    gsSparseSolver<>::SimplicialLDLT solver(assembler->matrix());
    gsMatrix<> solVector = solver.solve(assembler->rhs());
#endif

    assembler->constructSolution(solVector,assembler->allFixedDofs(),ALEdisp);
    if (m_options.getSwitch("Check"))
        return checkDisplacement(assembler->patches(),ALEdisp);
    else
        return -1;
}

template <class T>
index_t gsALE<T>::linearIncrementalMethod()
{
    for (size_t i = 0; i < m_interface.sidesA.size(); ++i)
        assembler->setFixedDofs(m_interface.sidesB[i].patch,
                                m_interface.sidesB[i].side(),
                                disp.patch(m_interface.sidesA[i].patch).boundary(m_interface.sidesA[i].side())->coefs() -
                                ALEdisp.patch(m_interface.sidesB[i].patch).boundary(m_interface.sidesB[i].side())->coefs(),
                                methodALE == ale_method::ILE ? false : true);
    assembler->assemble();

#ifdef GISMO_WITH_PARDISO
    gsSparseSolver<>::PardisoLDLT solver(assembler->matrix());
    gsMatrix<> solVector = solver.solve(assembler->rhs());
#else
    gsSparseSolver<>::SimplicialLDLT solver(assembler->matrix());
    gsMatrix<> solVector = solver.solve(assembler->rhs());
#endif

    gsMultiPatch<T> ALEupdate;
    assembler->constructSolution(solVector,assembler->allFixedDofs(),ALEupdate);
    for (size_t p = 0; p < ALEupdate.nPatches(); ++p)
    {
        ALEdisp.patch(p).coefs() += ALEupdate.patch(p).coefs();
        assembler->patches().patch(p).coefs() += ALEupdate.patch(p).coefs();
    }
    if (m_options.getSwitch("Check"))
        return checkGeometry(assembler->patches());
    else
        return -1;
}

template <class T>
index_t gsALE<T>::nonlinearMethod()
{
    for (size_t i = 0; i < m_interface.sidesA.size(); ++i)
        assembler->setFixedDofs(m_interface.sidesB[i].patch,
                                m_interface.sidesB[i].side(),
                                disp.patch(m_interface.sidesA[i].patch).boundary(m_interface.sidesA[i].side())->coefs() -
                                ALEdisp.patch(m_interface.sidesB[i].patch).boundary(m_interface.sidesB[i].side())->coefs());
    solverNL->reset();
    solverNL->solve();
    assembler->constructSolution(solverNL->solution(),solverNL->allFixedDofs(),ALEdisp);

    if (m_options.getSwitch("Check"))
        return checkDisplacement(assembler->patches(),ALEdisp);
    else
        return -1;
}

template <class T>
void gsALE<T>::saveState()
{
    if (methodALE == ale_method::TINE || methodALE == ale_method::TINE_StVK)
        solverNL->saveState();
    ALEdispSaved.clear();
    for (size_t p = 0; p < ALEdisp.nPatches(); ++p)
        ALEdispSaved.addPatch(ALEdisp.patch(p).clone());
    hasSavedState = true;
}

template <class T>
void gsALE<T>::recoverState()
{
    GISMO_ENSURE(hasSavedState,"No state saved!");
    if (methodALE == ale_method::TINE || methodALE == ale_method::TINE_StVK)
        solverNL->recoverState();
    if (methodALE == ale_method::IHE || methodALE == ale_method::ILE || methodALE == ale_method::IBHE)
        for (size_t p = 0; p < ALEdisp.nPatches(); ++p)
            assembler->patches().patch(p).coefs() += ALEdispSaved.patch(p).coefs() - ALEdisp.patch(p).coefs();
    for (size_t p = 0; p < ALEdisp.nPatches(); ++p)
        ALEdisp.patch(p).coefs() = ALEdispSaved.patch(p).coefs();
}

} // namespace ends
