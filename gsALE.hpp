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
#include <gsElasticity/gsIterative.h>
#include <gsElasticity/gsGeoUtils.h>

namespace gismo
{

template <class T>
gsALE<T>::gsALE(gsMultiPatch<T> & geometry, const gsMultiPatch<T> & displacement,
                const gsInterfaceFSI & interface, ale_method::method method)
    : disp(displacement),
      fsiInterface(interface),
      methodALE(method),
      m_options(defaultOptions()),
      hasSavedState(false)
{
    // create input for the assembler
    gsMultiBasis<T> basis(geometry);
    gsBoundaryConditions<T> bcInfo;
    for (auto it = geometry.bBegin(); it != geometry.bEnd(); ++it)
        for (index_t d = 0; d < geometry.parDim(); ++d)
            bcInfo.addCondition(it->patch,it->side(),condition_type::dirichlet,0,d);
    gsConstantFunction<T> rhs(gsVector<>::Zero(geometry.parDim()),geometry.parDim());

    // define assembler according to the method
    if (methodALE == ale_method::TINE || methodALE == ale_method::ILE || methodALE == ale_method::LE)
    {
        assembler = typename gsBaseAssembler<T>::uPtr(new gsElasticityAssembler<T>(geometry,basis,bcInfo,rhs));
        if (methodALE == ale_method::TINE)
        {
        assembler->options().setInt("MaterialLaw",material_law::neo_hooke_ln);
        assembler->options().setSwitch("Check",false);
        solverNL = typename gsIterative<T>::uPtr(new gsIterative<T>(*assembler,
                                                                    gsMatrix<>::Zero(assembler->numDofs(),1),
                                                                    assembler->allFixedDofs()));
        solverNL->options().setInt("Verbosity",solver_verbosity::none);
        solverNL->options().setInt("MaxIters",1);
        solverNL->options().setInt("Solver",linear_solver::LDLT);
        }
        assembler->constructSolution(gsMatrix<T>::Zero(assembler->numDofs(),1),assembler->allFixedDofs(),ALEdisp);
    }
    else if (methodALE == ale_method::HE || methodALE == ale_method::IHE)
    {
        assembler = typename gsBaseAssembler<T>::uPtr(new gsElPoissonAssembler<T>(geometry,basis,bcInfo,rhs));
        assembler->constructSolution(gsMatrix<T>::Zero(assembler->numDofs(),geometry.parDim()),assembler->allFixedDofs(),ALEdisp);
    }
    else
        for (index_t p = 0; p < geometry.nPatches(); ++p)
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
    return opt;
}


template <class T>
void gsALE<T>::constructSolution(gsMultiPatch<T> & solution) const
{
    solution.clear();
    for (index_t p = 0; p < ALEdisp.nPatches(); ++p)
        solution.addPatch(ALEdisp.patch(p).clone());
}

template <class T>
index_t gsALE<T>::updateMesh()
{
    switch (methodALE)
    {
        case ale_method::TINE: return TINE(); break;
        case ale_method::ILE: return ILE(); break;
        case ale_method::LE: return LE(); break;
        case ale_method::HE: return HE(); break;
        case ale_method::IHE: return IHE(); break;
        case ale_method::BHE: return BHE(); break;
        default: return -1;
    }
}

template <class T>
index_t gsALE<T>::TINE()
{
    assembler->options().setReal("PoissonsRatio",m_options.getReal("PoissonsRatio"));
    assembler->options().setReal("LocalStiff",m_options.getReal("LocalStiff"));
    for (index_t i = 0; i < fsiInterface.fluidSides.size(); ++i)
        assembler->setFixedDofs(fsiInterface.fluidSides[i].patch,
                                fsiInterface.fluidSides[i].side(),
                                disp.patch(fsiInterface.solidSides[i].patch).boundary(fsiInterface.solidSides[i].side())->coefs()-
                                ALEdisp.patch(fsiInterface.fluidSides[i].patch).boundary(fsiInterface.fluidSides[i].side())->coefs());
    solverNL->reset();
    solverNL->solve();
    assembler->constructSolution(solverNL->solution(),solverNL->allFixedDofs(),ALEdisp);
    if (m_options.getSwitch("Check"))
        return checkDisplacement(assembler->patches(),ALEdisp);
    else
        return -1;
}

template <class T>
index_t gsALE<T>::ILE()
{
    assembler->options().setReal("PoissonsRatio",m_options.getReal("PoissonsRatio"));
    assembler->options().setReal("LocalStiff",m_options.getReal("LocalStiff"));
    for (index_t i = 0; i < fsiInterface.fluidSides.size(); ++i)
        assembler->setFixedDofs(fsiInterface.fluidSides[i].patch,
                                fsiInterface.fluidSides[i].side(),
                                disp.patch(fsiInterface.solidSides[i].patch).boundary(fsiInterface.solidSides[i].side())->coefs()-
                                ALEdisp.patch(fsiInterface.fluidSides[i].patch).boundary(fsiInterface.fluidSides[i].side())->coefs());
    assembler->assemble();

#ifdef GISMO_WITH_PARDISO
    gsSparseSolver<>::PardisoLDLT solver(assembler->matrix());
    gsVector<> solVector = solver.solve(assembler->rhs());
#else
    gsSparseSolver<>::SimplicialLDLT solver(assembler->matrix());
    gsVector<> solVector = solver.solve(assembler->rhs());
#endif

    gsMultiPatch<T> ALEupdate;
    assembler->constructSolution(solVector,assembler->allFixedDofs(),ALEupdate);
    for (index_t p = 0; p < ALEupdate.nPatches(); ++p)
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
index_t gsALE<T>::LE()
{
    assembler->options().setReal("PoissonsRatio",m_options.getReal("PoissonsRatio"));
    assembler->options().setReal("LocalStiff",m_options.getReal("LocalStiff"));
    for (index_t i = 0; i < fsiInterface.fluidSides.size(); ++i)
        assembler->setFixedDofs(fsiInterface.fluidSides[i].patch,
                                fsiInterface.fluidSides[i].side(),
                                disp.patch(fsiInterface.solidSides[i].patch).boundary(fsiInterface.solidSides[i].side())->coefs());
    assembler->assemble();

#ifdef GISMO_WITH_PARDISO
    gsSparseSolver<>::PardisoLDLT solver(assembler->matrix());
    gsVector<> solVector = solver.solve(assembler->rhs());
#else
    gsSparseSolver<>::SimplicialLDLT solver(assembler->matrix());
    gsVector<> solVector = solver.solve(assembler->rhs());
#endif

    assembler->constructSolution(solVector,assembler->allFixedDofs(),ALEdisp);
    if (m_options.getSwitch("Check"))
        return checkDisplacement(assembler->patches(),ALEdisp);
    else
        return -1;
}

template <class T>
index_t gsALE<T>::HE()
{
    assembler->options().setReal("LocalStiff",m_options.getReal("LocalStiff"));
    for (index_t i = 0; i < fsiInterface.fluidSides.size(); ++i)
        assembler->setFixedDofs(fsiInterface.fluidSides[i].patch,
                                fsiInterface.fluidSides[i].side(),
                                disp.patch(fsiInterface.solidSides[i].patch).boundary(fsiInterface.solidSides[i].side())->coefs(),true);
    assembler->assemble();

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
index_t gsALE<T>::IHE()
{
    assembler->options().setReal("LocalStiff",m_options.getReal("LocalStiff"));
    for (index_t i = 0; i < fsiInterface.fluidSides.size(); ++i)
        assembler->setFixedDofs(fsiInterface.fluidSides[i].patch,
                                fsiInterface.fluidSides[i].side(),
                                disp.patch(fsiInterface.solidSides[i].patch).boundary(fsiInterface.solidSides[i].side())->coefs()-
                                ALEdisp.patch(fsiInterface.fluidSides[i].patch).boundary(fsiInterface.fluidSides[i].side())->coefs(),true);
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
    for (index_t p = 0; p < ALEupdate.nPatches(); ++p)
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
index_t gsALE<T>::BHE()
{
    GISMO_NO_IMPLEMENTATION;
    return checkDisplacement(assembler->patches(),ALEdisp);
}

template <class T>
void gsALE<T>::saveState()
{
    if (methodALE == ale_method::TINE)
        solverNL->saveState();
    ALEdispSaved.clear();
    for (index_t p = 0; p < ALEdisp.nPatches(); ++p)
        ALEdispSaved.addPatch(ALEdisp.patch(p).clone());
    hasSavedState = true;
}

template <class T>
void gsALE<T>::recoverState()
{
    GISMO_ENSURE(hasSavedState,"No state saved!");
    if (methodALE == ale_method::TINE)
        solverNL->recoverState();
    if (methodALE == ale_method::IHE || methodALE == ale_method::ILE)
        for (index_t p = 0; p < ALEdisp.nPatches(); ++p)
            assembler->patches().patch(p).coefs() += ALEdispSaved.patch(p).coefs() - ALEdisp.patch(p).coefs();
    for (index_t p = 0; p < ALEdisp.nPatches(); ++p)
        ALEdisp.patch(p).coefs() = ALEdispSaved.patch(p).coefs();
}

} // namespace ends
