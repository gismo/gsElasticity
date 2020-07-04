/** @file gsElasticityAssembler.hpp

    @brief Provides linear and nonlinear elasticity systems for 2D plain strain and 3D continua.

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

#include <gsUtils/gsPointGrid.h>
#include <gsElasticity/gsBaseUtils.h>
#include <gsElasticity/gsGeoUtils.h>
#include <gsElasticity/gsBasePde.h>

// Element visitors
#include <gsElasticity/gsVisitorLinearElasticity.h>
#include <gsElasticity/gsVisitorMixedLinearElasticity.h>
#include <gsElasticity/gsVisitorMixedNonLinearElasticity.h>
#include <gsElasticity/gsVisitorNonLinearElasticity.h>
#include <gsElasticity/gsVisitorElasticityNeumann.h>

namespace gismo
{

template<class T>
gsElasticityAssembler<T>::gsElasticityAssembler(const gsMultiPatch<T> & patches,
                                                const gsMultiBasis<T> & basis,
                                                const gsBoundaryConditions<T> & bconditions,
                                                const gsFunction<T> & body_force)
{
    // Originally concieved as a meaningful class, now gsPde is just a container for
    // the domain, boundary conditions and the right-hand side;
    // TUDO: change/remove gsPde from gsAssembler logic
    gsPiecewiseFunction<T> rightHandSides;
    rightHandSides.addPiece(body_force);
    typename gsPde<T>::Ptr pde( new gsBasePde<T>(patches,bconditions,rightHandSides) );
    // gsAssembler<>::initialize requires a vector of bases, one for each unknown;
    // different bases are used to compute Dirichlet DoFs;
    // but always the first basis is used for the assembly;
    // TODO: change gsAssembler logic
    m_dim = body_force.targetDim();
    for (short_t d = 0; d < m_dim; ++d)
        m_bases.push_back(basis);

    Base::initialize(pde, m_bases, defaultOptions());
}

template<class T>
gsElasticityAssembler<T>::gsElasticityAssembler(gsMultiPatch<T> const & patches,
                                                gsMultiBasis<T> const & basisDisp,
                                                gsMultiBasis<T> const & basisPres,
                                                gsBoundaryConditions<T> const & bconditions,
                                                const gsFunction<T> & body_force)
{
    // same as above
    gsPiecewiseFunction<T> rightHandSides;
    rightHandSides.addPiece(body_force);
    typename gsPde<T>::Ptr pde( new gsBasePde<T>(patches,bconditions,rightHandSides) );
    // same as above
    m_dim = body_force.targetDim();
    for (short_t d = 0; d < m_dim; ++d)
        m_bases.push_back(basisDisp);
    m_bases.push_back(basisPres);

    Base::initialize(pde, m_bases, defaultOptions());
    m_options.setInt("MaterialLaw",material_law::mixed_hooke);
}

template <class T>
gsOptionList gsElasticityAssembler<T>::defaultOptions()
{
    gsOptionList opt = Base::defaultOptions();
    opt.addReal("YoungsModulus","Youngs modulus of the material",200e9);
    opt.addReal("PoissonsRatio","Poisson's ratio of the material",0.33);
    opt.addReal("ForceScaling","Force scaling parameter",1.);
    opt.addInt("MaterialLaw","Material law: 0 for St. Venant-Kirchhof, 1 for Neo-Hooke",material_law::hooke);
    opt.addReal("LocalStiff","Stiffening degree for the Jacobian-based local stiffening",0.);
    opt.addSwitch("Check","Check bijectivity of the displacement field before matrix assebmly",false);
    return opt;
}

template <class T>
void gsElasticityAssembler<T>::reserve()
{
    // Pick up values from options
    const T bdA       = m_options.getReal("bdA");
    const index_t bdB = m_options.getInt("bdB");
    const T bdO       = m_options.getReal("bdO");

    index_t deg = 0;
    for (index_t d = 0; d < m_bases[0][0].dim(); ++d )
        if (m_bases[0][0].degree(d) > deg)
            deg = m_bases[0][0].degree(d);

    if (m_bases.size() == unsigned(m_dim)) // displacement formulation
    {
        // m_dim velocity*velocity blocks
        index_t numElPerColumn = pow((bdA*deg+bdB),m_dim)*m_dim;
        m_system.reserve(numElPerColumn*(1+bdO),1);
    }
    else // mixed formulation (displacement + pressure)
    {
        // m_dim velocity*velocity blocks + 1 pressure*velocity block (additioanal factor 2 for subgrid element)
        index_t numElPerColumn = pow((bdA*deg+bdB),m_dim)*m_dim + pow((2*bdA*deg+bdB),m_dim);
        m_system.reserve(numElPerColumn*(1+bdO),1);
    }
}

template <class T>
void gsElasticityAssembler<T>::refresh()
{
    GISMO_ENSURE(m_dim == m_pde_ptr->domain().parDim(), "The RHS dimension and the domain dimension don't match!");
    GISMO_ENSURE(m_dim == 2 || m_dim == 3, "Only two- and three-dimenstion domains are supported!");

    std::vector<gsDofMapper> m_dofMappers(m_bases.size());
    for (unsigned d = 0; d < m_bases.size(); d++)
        m_bases[d].getMapper((dirichlet::strategy)m_options.getInt("DirichletStrategy"),
                             iFace::glue,m_pde_ptr->bc(),m_dofMappers[d],d,true);

    m_system = gsSparseSystem<T>(m_dofMappers, gsVector<index_t>::Ones(m_bases.size()));
    reserve();

    for (unsigned d = 0; d < m_bases.size(); ++d)
        Base::computeDirichletDofs(d);
}

//--------------------- SYSTEM ASSEMBLY ----------------------------------//

template<class T>
void gsElasticityAssembler<T>::assemble(bool saveEliminationMatrix)
{
    m_system.matrix().setZero();
    reserve();
    m_system.rhs().setZero();

    // Compute volumetric integrals and write to the global linear system
    if (m_bases.size() == unsigned(m_dim)) // displacement formulation
    {
        GISMO_ENSURE(m_options.getInt("MaterialLaw") == material_law::hooke,
                     "Material law not specified OR not supported!");
        if (saveEliminationMatrix)
        {
            eliminationMatrix.resize(Base::numDofs(),Base::numFixedDofs());
            eliminationMatrix.setZero();
            eliminationMatrix.reservePerColumn(m_system.numColNz(m_bases[0],m_options));
        }

        gsVisitorLinearElasticity<T> visitor(*m_pde_ptr, saveEliminationMatrix ? &eliminationMatrix : nullptr);
        Base::template push<gsVisitorLinearElasticity<T> >(visitor);

        if (saveEliminationMatrix)
        {
            Base::rhsWithZeroDDofs = m_system.rhs();
            eliminationMatrix.makeCompressed();
        }

    }
    else // mixed formulation (displacement + pressure)
    {
        GISMO_ENSURE(m_options.getInt("MaterialLaw") == material_law::mixed_hooke,
                     "Material law not specified OR not supported!");
        gsVisitorMixedLinearElasticity<T> visitor(*m_pde_ptr);
        Base::template push<gsVisitorMixedLinearElasticity<T> >(visitor);
    }

    // Compute surface integrals and write to the global rhs vector
    Base::template push<gsVisitorElasticityNeumann<T> >(m_pde_ptr->bc().neumannSides());

    m_system.matrix().makeCompressed();
}

template <class T>
bool gsElasticityAssembler<T>::assemble(const gsMatrix<T> & solutionVector,
                                        const std::vector<gsMatrix<T> > & fixedDoFs)
{
    gsMultiPatch<T> displacement;
    constructSolution(solutionVector,fixedDoFs,displacement);
    if (m_options.getSwitch("Check"))
        if (checkDisplacement(m_pde_ptr->patches(),displacement) != -1)
            return false;

    if (m_bases.size() == unsigned(m_dim)) // displacement formulation 
        assemble(displacement);
    else // mixed formulation (displacement + pressure)
    {
        gsMultiPatch<T> pressure;
        constructPressure(solutionVector,fixedDoFs,pressure);
        assemble(displacement,pressure);
    }
    return true;
}

template<class T>
void gsElasticityAssembler<T>::assemble(const gsMultiPatch<T> & displacement)
{
    GISMO_ENSURE(m_options.getInt("MaterialLaw") == material_law::saint_venant_kirchhoff ||
                 m_options.getInt("MaterialLaw") == material_law::neo_hooke_ln ||
                 m_options.getInt("MaterialLaw") == material_law::neo_hooke_quad,
                 "Material law not specified OR not supported!");
    m_system.matrix().setZero();
    reserve();
    m_system.rhs().setZero();

    // Compute volumetric integrals and write to the global linear system
    gsVisitorNonLinearElasticity<T> visitor(*m_pde_ptr,displacement);
    Base::template push<gsVisitorNonLinearElasticity<T> >(visitor);
    // Compute surface integrals and write to the global rhs vector
    // change to reuse rhs from linear system
    Base::template push<gsVisitorElasticityNeumann<T> >(m_pde_ptr->bc().neumannSides());

    m_system.matrix().makeCompressed();
}

template<class T>
void gsElasticityAssembler<T>::assemble(const gsMultiPatch<T> & displacement,
                                        const gsMultiPatch<T> & pressure)
{
    GISMO_ENSURE(m_options.getInt("MaterialLaw") == material_law::mixed_neo_hooke_ln,
                 "Material law not specified OR not supported!");
    m_options.setInt("MaterialLaw",material_law::mixed_neo_hooke_ln);
    m_system.matrix().setZero();
    reserve();
    m_system.rhs().setZero();

    // Compute volumetric integrals and write to the global linear systemz
    gsVisitorMixedNonLinearElasticity<T> visitor(*m_pde_ptr,displacement,pressure);
    Base::template push<gsVisitorMixedNonLinearElasticity<T> >(visitor);
    // Compute surface integrals and write to the global rhs vector
    // change to reuse rhs from linear system
    Base::template push<gsVisitorElasticityNeumann<T> >(m_pde_ptr->bc().neumannSides());

    m_system.matrix().makeCompressed();
}

//--------------------- SOLUTION CONSTRUCTION ----------------------------------//

template <class T>
void gsElasticityAssembler<T>::constructSolution(const gsMatrix<T> & solVector,
                                                 const std::vector<gsMatrix<T> > & fixedDoFs,
                                                 gsMultiPatch<T> & result) const
{
    gsVector<index_t> unknowns(m_dim);
    for (short_t d = 0; d < m_dim; ++d)
        unknowns.at(d) = d;
    Base::constructSolution(solVector,fixedDoFs,result,unknowns);
}

template <class T>
void gsElasticityAssembler<T>::constructSolution(const gsMatrix<T>& solVector,
                                                 const std::vector<gsMatrix<T> > & fixedDoFs,
                                                 gsMultiPatch<T>  & displacement, gsMultiPatch<T> & pressure) const
{
    GISMO_ENSURE(m_bases.size() == unsigned(m_dim) + 1, "Not a mixed formulation: can't construct pressure.");
    // construct displacement
    constructSolution(solVector,fixedDoFs,displacement);
    // construct pressure
    constructPressure(solVector,fixedDoFs,pressure);
}

template <class T>
void gsElasticityAssembler<T>::constructPressure(const gsMatrix<T>& solVector,
                                                 const std::vector<gsMatrix<T> > & fixedDoFs,
                                                 gsMultiPatch<T>& pressure) const
{
    GISMO_ENSURE(m_bases.size() == unsigned(m_dim) + 1, "Not a mixed formulation: can't construct pressure.");
    gsVector<index_t> unknowns(1);
    unknowns.at(0) = m_dim;
    Base::constructSolution(solVector,fixedDoFs,pressure,unknowns);
}

//--------------------- SPECIALS ----------------------------------//

template <class T>
void gsElasticityAssembler<T>::constructCauchyStresses(const gsMultiPatch<T> & displacement,
                                                       gsPiecewiseFunction<T> & result,
                                                       stress_components::components comp) const
{
    if (comp == stress_components::all_2D_vector || comp == stress_components::all_2D_matrix)
        GISMO_ENSURE(m_dim == 2, "Invalid stress components for a 2D problem");
    if (comp == stress_components::normal_3D_vector || comp == stress_components::shear_3D_vector ||
        comp == stress_components::all_3D_matrix)
        GISMO_ENSURE(m_dim == 3, "Invalid stress type for a 3D problem");
    GISMO_ENSURE(m_options.getInt("MaterialLaw") == material_law::hooke ||
                 m_options.getInt("MaterialLaw") == material_law::neo_hooke_ln ||
                 m_options.getInt("MaterialLaw") == material_law::saint_venant_kirchhoff ||
                 m_options.getInt("MaterialLaw") == material_law::neo_hooke_quad,
                 "Pressure field not provided! Can't compute stresses with the chosen material law.");
    result.clear();

    for (size_t p = 0; p < m_pde_ptr->domain().nPatches(); ++p )
        result.addPiecePointer(new gsCauchyStressFunction<T>(p,comp,m_options,
                                                             &(m_pde_ptr->domain()),&displacement));
}

template <class T>
void gsElasticityAssembler<T>::constructCauchyStresses(const gsMultiPatch<T> & displacement,
                                                       const gsMultiPatch<T> & pressure,
                                                       gsPiecewiseFunction<T> & result,
                                                       stress_components::components comp) const
{
    if (comp == stress_components::all_2D_vector || comp == stress_components::all_2D_matrix)
        GISMO_ENSURE(m_dim == 2, "Invalid stress components for a 2D problem");
    if (comp == stress_components::normal_3D_vector || comp == stress_components::shear_3D_vector ||
        comp == stress_components::all_3D_matrix)
        GISMO_ENSURE(m_dim == 3, "Invalid stress type for a 3D problem");
    GISMO_ENSURE(m_options.getInt("MaterialLaw") == material_law::mixed_neo_hooke_ln ||
                 m_options.getInt("MaterialLaw") == material_law::mixed_hooke,
                 "Pressure field is not necessary to compute stresses with the chosen material law.");
    result.clear();

    for (size_t p = 0; p < m_pde_ptr->domain().nPatches(); ++p )
        result.addPiecePointer(new gsCauchyStressFunction<T>(p,comp,m_options,
                                                             &(m_pde_ptr->domain()), &displacement, &pressure));
}

}// namespace gismo ends
