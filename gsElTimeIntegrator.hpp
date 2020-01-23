/** @file gsElTimeIntegrator.hpp

    @brief A class providing time integration for dynamical elasticity.

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s):
        O. Weeger    (2012 - 2015, TU Kaiserslautern),
        A.Shamanskiy (2016 - ...., TU Kaiserslautern)
*/

#pragma once

#include <gsElasticity/gsElTimeIntegrator.h>

#include <gsElasticity/gsElasticityAssembler.h>
#include <gsElasticity/gsMassAssembler.h>
#include <gsElasticity/gsIterative.h>


namespace gismo
{

template <class T>
gsElTimeIntegrator<T>::gsElTimeIntegrator(gsElasticityAssembler<T> & stiffAssembler_,
                                          gsMassAssembler<T> & massAssembler_)
    : stiffAssembler(stiffAssembler_),
      massAssembler(massAssembler_)
{
    m_options = defaultOptions();
    m_ddof = stiffAssembler.allFixedDofs();
}

template <class T>
gsOptionList gsElTimeIntegrator<T>::defaultOptions()
{
    gsOptionList opt = Base::defaultOptions();
    opt.addInt("Scheme","Time integration scheme",time_integration::implicit_linear);
    opt.addReal("Beta","Parameter beta for the time integration scheme, see Wriggers, Nonlinear FEM, p.213 ",0.25);
    opt.addReal("Gamma","Parameter gamma for the time integration scheme, see Wriggers, Nonlinear FEM, p.213 ",0.5);
    opt.addInt("Verbosity","Amount of information printed to the terminal: none, some, all",solver_verbosity::none);
    return opt;
}

template <class T>
void gsElTimeIntegrator<T>::initialize()
{
    massAssembler.assemble();

    if (dispVector.rows() != stiffAssembler.numDofs())
        dispVector.setZero(stiffAssembler.numDofs(),1);
    if (velVector.rows() != stiffAssembler.numDofs())
        velVector.setZero(stiffAssembler.numDofs(),1);

    stiffAssembler.assemble(dispVector,m_ddof,false);

    gsSparseSolver<>::SimplicialLDLT solver(massAssembler.matrix());
    accVector = solver.solve(stiffAssembler.rhs());
}

template <class T>
void gsElTimeIntegrator<T>::makeTimeStep(T timeStep)
{
    tStep = timeStep;
    gsMatrix<T> newDispVector;
    if (m_options.getInt("Scheme") == time_integration::implicit_linear)
        newDispVector = implicitLinear();
    if (m_options.getInt("Scheme") == time_integration::implicit_nonlinear)
        newDispVector = implicitNonlinear();
    gsMatrix<T> tempVelVector = velVector;
    velVector = alpha4()*(newDispVector - dispVector) + alpha5()*tempVelVector + alpha6()*accVector;
    accVector = alpha1()*(newDispVector - dispVector) - alpha2()*tempVelVector - alpha3()*accVector;
    dispVector = newDispVector;
}

template <class T>
gsMatrix<T> gsElTimeIntegrator<T>::implicitLinear()
{
    m_matrix = alpha1()*massAssembler.matrix() + stiffAssembler.matrix();
    m_matrix.makeCompressed();
    m_rhs = massAssembler.matrix()*(alpha1()*dispVector + alpha2()*velVector + alpha3()*accVector) + stiffAssembler.rhs();
#ifdef GISMO_WITH_PARDISO
    gsSparseSolver<>::PardisoLDLT solver(m_matrix);
    return solver.solve(m_rhs);
#else
    gsSparseSolver<>::SimplicialLDLT solver(m_matrix);
    return solver.solve(m_rhs);
#endif
}

template <class T>
gsMatrix<T> gsElTimeIntegrator<T>::implicitNonlinear()
{
    gsIterative<T> solver(*this,dispVector);
    solver.options().setInt("Verbosity",m_options.getInt("Verbosity"));
    solver.options().setInt("Solver",linear_solver::LDLT);
    solver.solve();
    return solver.solution();
}

template <class T>
bool gsElTimeIntegrator<T>::assemble(const gsMatrix<T> & solutionVector,
                                     const std::vector<gsMatrix<T> > & fixedDoFs,
                                     bool assembleMatrix)
{
    stiffAssembler.assemble(solutionVector,fixedDoFs,assembleMatrix);
    Base::m_system.matrix() = alpha1()*massAssembler.matrix() + stiffAssembler.matrix();
    Base::m_system.matrix().makeCompressed();
    Base::m_system.rhs() = stiffAssembler.rhs() + massAssembler.matrix()*(alpha1()*(dispVector-solutionVector) + alpha2()*velVector + alpha3()*accVector);
    return true;
}

template <class T>
void gsElTimeIntegrator<T>::saveState()
{
    dispVecSaved = dispVector;
    velVecSaved = velVector;
    accVecSaved = accVector;
    ddofsSaved = m_ddof;
}

template <class T>
void gsElTimeIntegrator<T>::recoverState()
{
    GISMO_ASSERT(dispVecSaved.rows() == dispVector.rows(),"No state saved!");
    dispVector = dispVecSaved;
    velVector = velVecSaved;
    accVector = accVecSaved;
    m_ddof = ddofsSaved;
}


} // namespace ends
