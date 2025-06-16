/** @file gsSolidAssembler.hpp

    @brief Provides assemblers for elasticity problems

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s):
        H.M. Verhelst
*/

#pragma once

#include <gsElasticity/gsSolidAssembler.h>
#include <gsPde/gsBoundaryConditions.h>
#include <gsElasticity/material_expr.h>
#include <gsElasticity/voigt_expr.h>

namespace gismo
{

template <short_t DIM, class T, class Material>
gsSolidAssembler<DIM,T,Material>::gsSolidAssembler(const gsMultiPatch<T> & mp,
                                                   const gsMultiBasis<T> & mb,
                                                   const gsBoundaryConditions<T> & bcs,
                                                   const Material *         material,
                                                   const gsFunctionSet<T> * g)
:
m_patches(mp),
m_basis(mb),
m_spaceBasis(&mb),
m_bcs(bcs),
m_forcing(g),
m_material(material),
m_initialized(false)
{
    this->_defaultOptions();
};

template <short_t DIM, class T, class Material>
gsSolidAssembler<DIM,T,Material>& gsSolidAssembler<DIM,T,Material>::operator=( const gsSolidAssembler& other )
{
    if (this!=&other)
    {
        m_continuity=other.m_continuity;

        m_patches=other.m_patches;
        m_basis=other.m_basis;
        m_spaceBasis=other.m_spaceBasis;
        m_bcs=other.m_bcs;

        m_options=other.m_options;

        // To do: make copy constructor for the gsExprAssembler
        m_assembler.setIntegrationElements(m_basis);
        m_assembler.setOptions(m_options);
    }
    return *this;
}

template <short_t DIM, class T, class Material>
gsSolidAssembler<DIM,T,Material>& gsSolidAssembler<DIM,T,Material>::operator=( gsSolidAssembler&& other )
{
    m_continuity=give(other.m_continuity);

    m_patches=give(other.m_patches);
    m_basis=give(other.m_basis);
    m_spaceBasis=give(other.m_spaceBasis);
    m_bcs=give(other.m_bcs);

    m_options=give(other.m_options);

    // To do: make copy constructor for the gsExprAssembler
    m_assembler.setIntegrationElements(m_basis);
    m_assembler.setOptions(m_options);
    return *this;
}

template <short_t DIM, class T, class Material>
void gsSolidAssembler<DIM,T,Material>::_defaultOptions()
{
    m_options.addInt("Continuity","Continuity between patches: C^{-1} (-1) or C^0 (0, default)",0);
    m_options.addSwitch("SmallStrain","Small strain formulation: 0 (no) or 1 (yes, default)",1);
    // Assembler options
    gsOptionList assemblerOptions = m_assembler.defaultOptions().wrapIntoGroup("ExprAssembler");
    m_options.update(assemblerOptions,gsOptionList::addIfUnknown);
}

template <short_t DIM, class T, class Material>
void gsSolidAssembler<DIM,T,Material>::_getOptions()
{
    m_continuity = m_options.getInt("Continuity");
    m_smallStrain = m_options.getSwitch("SmallStrain");

    GISMO_ENSURE(m_options.hasGroup("ExprAssembler"),"The option list does not contain options with the label 'ExprAssembler'!");
    m_assembler.setOptions(m_options.getGroup("ExprAssembler"));
}

template <short_t DIM, class T, class Material>
void gsSolidAssembler<DIM,T,Material>::setOptions(gsOptionList & options)
{
    m_options.update(options,gsOptionList::ignoreIfUnknown);
}

template <short_t DIM, class T, class Material>
void gsSolidAssembler<DIM,T,Material>::initialize()
{
    this->_getOptions();

    // Elements used for numerical integration
    m_assembler.setIntegrationElements(m_basis);
    m_assembler.setOptions(m_options.getGroup("ExprAssembler"));

    GISMO_ASSERT(m_bcs.hasGeoMap(),"No geometry map was assigned to the boundary conditions. Use bc.setGeoMap to assign one!");

    geometryMap G = m_assembler.getMap(m_patches);
    space u = m_assembler.getSpace(*m_spaceBasis, DIM, 0); // last argument is the space ID

    // Setup the system
    u.setup(m_bcs,
            m_assembler.options().askInt("DirichletValues",dirichlet::interpolation),
            m_options.askInt("Continuity",-1));
    m_assembler.initSystem();
    // Compute sparsity patter: this is done automatically - but
    // is needed if assemble(.) is called twice
    auto delta_FF = ijac(u, G);
    auto delta_EE = 0.5*(delta_FF.cwisetr() + delta_FF); //  check if symmetrize is possible
    auto delta_EE_voigt = gismo::expr::voigt<DIM>(delta_EE);
    m_assembler.computePattern(delta_EE_voigt*delta_EE_voigt.tr());
    m_initialized = true;
}

template <short_t DIM, class T, class Material>
void gsSolidAssembler<DIM,T,Material>::assemble()
{
    GISMO_ENSURE(m_initialized,"The assembler has not been initialized yet. Call initialize() before assembling the system.");
    GISMO_ENSURE(m_material,"The material has not been set yet. Call setMaterial() before assembling the system.");
    m_assembler.clearMatrix(); // Resets to zero the values of the already allocated to matrix (LHS)
    m_assembler.clearRhs(); // Resets to zero the values of the already allocated to rhs (RHS)

    // Set the geometry map
    geometryMap G = m_assembler.getMap(m_patches);
    space       u = m_assembler.trialSpace(0);
    auto        f = m_assembler.getCoeff(*m_forcing);

    // Unpack the material parameters
    index_t numParams = m_material->numParameters();
    std::vector<expr::gsFeVariable<T>> materialParams;
    materialParams.reserve(numParams);
    for (size_t i = 0; i < m_material->numParameters(); ++i)
        materialParams.push_back(m_assembler.getCoeff(*m_material->getParameter(i)));

    // Define expressions
    auto delta_FF = ijac(u, G);
    auto delta_EE = 0.5*(delta_FF.cwisetr() + delta_FF); //  check if symmetrize is possible
    auto delta_EE_voigt = gismo::expr::voigt<DIM>(delta_EE);

    auto C = gismo::expr::material<Material, false, gsMaterialOutput::C>(G,materialParams);
    auto bilinear_form = delta_EE_voigt*C*delta_EE_voigt.tr();
    // Assemble the system
    if (m_forcing == nullptr)
        m_assembler.assemble(bilinear_form*meas(G));
    else
    {
        auto f = m_assembler.getCoeff(*m_forcing);
        m_assembler.assemble(bilinear_form*meas(G),u*f*meas(G));
    }
    auto g_N = m_assembler.getBdrFunction(G);
    m_assembler.assembleBdr(m_bcs.get("Neumann"), u * g_N * nv(G).norm() );
}

template <short_t DIM, class T, class Material>
void gsSolidAssembler<DIM,T,Material>::assemble(const gsMatrix<T> & uvec)
{
    GISMO_ENSURE(m_initialized,"The assembler has not been initialized yet. Call initialize() before assembling the system.");
    GISMO_ENSURE(m_material,"The material has not been set yet. Call setMaterial() before assembling the system.");
    GISMO_ASSERT(uvec.rows() == m_assembler.numDofs(),"The size of the input vector does not match the number of degrees of freedom in the assembler.");
    gsMatrix<T> uvec_copy = uvec; // Make a copy of the input vector to avoid modifying the original one
    m_assembler.clearMatrix(); // Resets to zero the values of the already allocated to matrix (LHS)
    m_assembler.clearRhs(); // Resets to zero the values of the already allocated to rhs (RHS)

    // Set the geometry map
    geometryMap G = m_assembler.getMap(m_patches);
    space       u = m_assembler.trialSpace(0);
    solution  sol = m_assembler.getSolution(u, uvec_copy);

    // Unpack the material parameters
    index_t numParams = m_material->numParameters();
    std::vector<expr::gsFeVariable<T>> materialParams;
    materialParams.reserve(numParams);
    for (size_t i = 0; i < m_material->numParameters(); ++i)
        materialParams.push_back(m_assembler.getCoeff(*m_material->getParameter(i)));

    // Define expressions
    auto delta_FF = ijac(u, G);
    auto delta_EE = 0.5*(delta_FF.cwisetr() + delta_FF); //  check if symmetrize is possible
    auto delta_EE_voigt = gismo::expr::voigt<DIM>(delta_EE);

    auto C = gismo::expr::material<Material, true, gsMaterialOutput::C>(G,sol,materialParams);
    auto bilinear_form = delta_EE_voigt*C*delta_EE_voigt.tr();

    gsVector<T> pt(2);
    pt.setConstant(0.5);
    gsExprEvaluator<T> ev(m_assembler);
    if (m_forcing == nullptr)
        m_assembler.assemble(bilinear_form*meas(G));
    else
    {
        auto f = m_assembler.getCoeff(*m_forcing);
        m_assembler.assemble(bilinear_form*meas(G),u*f*meas(G));
    }
    auto g_N = m_assembler.getBdrFunction(G);
    m_assembler.assembleBdr(m_bcs.get("Neumann"), u * g_N * nv(G).norm() );
}

template <short_t DIM, class T, class Material>
void gsSolidAssembler<DIM,T,Material>::assembleMass()
{
    GISMO_ENSURE(m_initialized,"The assembler has not been initialized yet. Call initialize() before assembling the system.");
    GISMO_ENSURE(m_material,"The material has not been set yet. Call setMaterial() before assembling the system.");
    GISMO_ASSERT(m_material->hasDensity(),"The material does not have a density function defined. Please set the density function using setDensity() method.");

    m_assembler.clearMatrix(); // Resets to zero the values of the already allocated to matrix (LHS)

    // Set the geometry map
    geometryMap G = m_assembler.getMap(m_patches);
    space       u = m_assembler.trialSpace(0);
    u.setup(m_bcs, dirichlet::homogeneous, m_options.askInt("Continuity",-1));

    // Define expressions
    auto rho = m_assembler.getCoeff(*m_material->getDensity());
    m_assembler.assemble(u*u.tr()*meas(G));
}

template <short_t DIM, class T, class Material>
void gsSolidAssembler<DIM,T,Material>::constructSolution(gsMatrix<T> & uvec,
                                                         gsMultiPatch<T> & displacement) const
{
    auto w  = m_assembler.trialSpace(0);
    auto u  = m_assembler.getSolution(w,  uvec);
    u.extract(displacement);
}

template <short_t DIM, class T, class Material>
void gsSolidAssembler<DIM,T,Material>::constructDeformed(gsMatrix<T> & uvec,
                                                         gsMultiPatch<T> & deformed) const
{
    this->constructSolution(uvec, deformed);
    GISMO_ASSERT(deformed.nPatches() == m_patches.nPatches(), "The number of patches in the deformed solution does not match the number of patches in the original geometry.");
    // Apply the deformation to the patches
    for (size_t p = 0; p < m_patches.nPatches(); ++p)
    {
        GISMO_ASSERT(m_patches.patch(p).size()==deformed.patch(p).size(),"The size of the patch in the deformed solution does not match the size of the original patch.");
        deformed.patch(p).coefs()+= m_patches.patch(p).coefs();
    }
}

template <short_t DIM, class T, class Material>
void gsSolidAssembler<DIM,T,Material>::constructSolution(const gsMultiPatch<T> & displacement,
                                                               gsMatrix<T>     & uvec) const
{
    GISMO_ASSERT(displacement.geoDim()==DIM,"The displacement must be a vector field of dimension DIM.");
    GISMO_ASSERT(displacement.nPatches()==m_basis.nBases(),"Number of patches in the displacement must be equal to the number of bases in the assembler");
    auto w  = m_assembler.trialSpace(0);

    uvec.setZero(this->numDofs(),1);
    for (size_t b=0; b!=m_basis.nBases(); b++)
    {
        for (index_t i = 0; i < m_basis.basis(b).size(); i++)
        {
            GISMO_ASSERT(displacement.basis(b).size()==m_basis.basis(b).size(),"Number of basis functions in the displacement must be equal to the number of basis functions in the assembler");
            for (short_t d = 0; d < DIM; ++d) // Loop over the spatial dimensions
            {
                if (w.mapper().is_free(i,b))
                    uvec(w.mapper().index(i,b)) = displacement.patch(b).coefs()(i,d);
            }
        }
    }
}

template <short_t DIM, class T, class Material>
void gsSolidAssembler<DIM,T,Material>::setSpaceBasis(const gsFunctionSet<T> & spaceBasis)
{
    m_spaceBasis = &spaceBasis;
}

}// namespace gismo
