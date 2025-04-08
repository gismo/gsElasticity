/** @file gsPhaseFieldAssembler.hpp

    @brief Provides assemblers for the Cahn-Hilliard equation

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s):
        H.M. Verhelst   (2019-..., TU Delft)
        m_assembler. Mantzaflaris (2019-..., Inria)
*/

#pragma once

#include <gsElasticity/gsPhaseFieldAssembler.h>
#include <gsMSplines/gsMappedBasis.h>
#include <gsMSplines/gsMappedSpline.h>
#include <gsPde/gsBoundaryConditions.h>
#include <gsCore/gsFunctionExpr.h>
#include <gsCore/gsBasis.h>
#include <gsCore/gsPiecewiseFunction.h>

namespace gismo
{

template <class T, enum PForder order, enum PFmode mode>
gsPhaseFieldAssembler<T,order,mode>::gsPhaseFieldAssembler(const gsMultiPatch<T> & mp,
                                                const gsMultiBasis<T> & mb,
                                                const gsBoundaryConditions<T> & bcs
                                                    )
:
m_patches(mp),
m_basis(mb),
m_spaceBasis(&mb),
m_bcs(bcs),
m_initialized(false)
{
    this->_defaultOptions();
};

template <class T, enum PForder order, enum PFmode mode>
gsPhaseFieldAssembler<T,order,mode>& gsPhaseFieldAssembler<T,order,mode>::operator=( const gsPhaseFieldAssembler& other )
{
    if (this!=&other)
    {
        m_l0=other.m_l0;
        m_Gc=other.m_Gc;
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

template <class T, enum PForder order, enum PFmode mode>
gsPhaseFieldAssembler<T,order,mode>& gsPhaseFieldAssembler<T,order,mode>::operator=( gsPhaseFieldAssembler&& other )
{
    m_l0=give(other.m_l0);
    m_Gc=give(other.m_Gc);
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

template <class T, enum PForder order, enum PFmode mode>
void gsPhaseFieldAssembler<T,order,mode>::_defaultOptions()
{
    m_options.addReal("l0","l0 parameter",1e-4);
    m_options.addReal("cw","cw parameter",1e0);
    m_options.addReal("Gc","Gc parameter",1e0);
    m_options.addInt("Continuity","Continuity between patches: C^{-1} (-1) or C^0 (0, default)",0);
    // m_options.addReal("Penalty","Penalty parameter for Nitsche boundary conditions (default: 1e4)",1e4);
    // m_options.addSwitch("AssembleWeakBCs","Assemble Nitsche boundary conditions in every iteration",false);

    /* UNUSED:
    m_options.addInt("Mobility","Mobility function: 0 for constant, 1 for double well",0);
    m_options.addReal("M0","M0 parameter",1e0);
     */
    // Assembler options
    gsOptionList assemblerOptions = m_assembler.defaultOptions().wrapIntoGroup("ExprAssembler");
    m_options.update(assemblerOptions,gsOptionList::addIfUnknown);
}

template <class T, enum PForder order, enum PFmode mode>
void gsPhaseFieldAssembler<T,order,mode>::_getOptions()
{
    m_l0 = m_options.getReal("l0");
    m_cw = m_options.getReal("cw");
    m_Gc = m_options.getReal("Gc");
    m_continuity = m_options.getInt("Continuity");

    GISMO_ENSURE(m_options.hasGroup("ExprAssembler"),"The option list does not contain options with the label 'ExprAssembler'!");
    m_assembler.setOptions(m_options.getGroup("ExprAssembler"));
}

template <class T, enum PForder order, enum PFmode mode>
void gsPhaseFieldAssembler<T,order,mode>::setOptions(gsOptionList & options)
{
    m_options.update(options,gsOptionList::ignoreIfUnknown);
}

template <class T, enum PForder order, enum PFmode mode>
void gsPhaseFieldAssembler<T,order,mode>::initialize()
{
    this->_getOptions();

    // Elements used for numerical integration
    m_assembler.setIntegrationElements(m_basis);
    m_assembler.setOptions(m_options.getGroup("ExprAssembler"));

    GISMO_ASSERT(m_bcs.hasGeoMap(),"No geometry map was assigned to the boundary conditions. Use bc.setGeoMap to assign one!");

    typename gsExprAssembler<T>::space u = m_assembler.getSpace(*m_spaceBasis, 1, 0); // last argument is the space ID

    // Setup the system
    u.setup(m_bcs,
            m_assembler.options().askInt("DirichletValues",dirichlet::l2Projection),
            m_options.askInt("Continuity",-1));
    m_assembler.initSystem();

    // Compute sparsity patter: this is done automatically - but
    // is needed if assemble(.) is called twice
    m_assembler.computePattern( igrad(u) * igrad(u).tr() );

    m_initialized = true;
}

template <class T, enum PForder order, enum PFmode mode>
void gsPhaseFieldAssembler<T,order,mode>::assembleResidual(const gsFunctionSet<T> & C, const gsFunctionSet<T> & DC)
{
    assembleResidual_impl<order,mode>(C,DC);
}

template <class T, enum PForder order, enum PFmode mode>
template <enum PForder _order, enum PFmode _mode>
typename std::enable_if<(_order==PForder::Second && _mode==PFmode::AT1), void>::type
gsPhaseFieldAssembler<T,order,mode>::assembleResidual_impl(const gsFunctionSet<T> & C, const gsFunctionSet<T> & DC)
{
    GISMO_NO_IMPLEMENTATION;
}

template <class T, enum PForder order, enum PFmode mode>
template <enum PForder _order, enum PFmode _mode>
typename std::enable_if<(_order==PForder::Fourth && _mode==PFmode::AT1), void>::type
gsPhaseFieldAssembler<T,order,mode>::assembleResidual_impl(const gsFunctionSet<T> & C, const gsFunctionSet<T> & DC)
{
    GISMO_NO_IMPLEMENTATION;
}

template <class T, enum PForder order, enum PFmode mode>
template <enum PForder _order, enum PFmode _mode>
typename std::enable_if<(_order==PForder::Second && _mode==PFmode::AT2), void>::type
gsPhaseFieldAssembler<T,order,mode>::assembleResidual_impl(const gsFunctionSet<T> & C, const gsFunctionSet<T> & DC)
{
    GISMO_UNUSED(DC);

    GISMO_ENSURE(m_initialized,"The assembler has not been initialized yet. Call initialize() before assembling the system.");
    m_assembler.clearRhs(); // Resets to zero the values of the already allocated to matrix (LHS)
    m_assembler.initSystem();
    // Set the geometry map
    geometryMap G = m_assembler.getMap(m_patches);

    // Get the solution and its derivative
    auto c  = m_assembler.getCoeff(C);
    // auto dc = m_assembler.getCoeff(DC);

    // Set the discretization space
    auto w = m_assembler.trialSpace(0);

    m_assembler.assemble(
                m_Gc / m_l0 *
                (
                w*c.tr() +
                math::pow(m_l0,2)*igrad(w,G) * igrad(c,G).tr()
                ) * meas(G)
                );
}

template <class T, enum PForder order, enum PFmode mode>
template <enum PForder _order, enum PFmode _mode>
typename std::enable_if<(_order==PForder::Fourth && _mode==PFmode::AT2), void>::type
gsPhaseFieldAssembler<T,order,mode>::assembleResidual_impl(const gsFunctionSet<T> & C, const gsFunctionSet<T> & DC)
{
    GISMO_UNUSED(DC);

    GISMO_ENSURE(m_initialized,"The assembler has not been initialized yet. Call initialize() before assembling the system.");
    m_assembler.clearRhs(); // Resets to zero the values of the already allocated to matrix (LHS)
    m_assembler.initSystem();

    // Set the geometry map
    geometryMap G = m_assembler.getMap(m_patches);

    // Get the solution and its derivative
    auto c  = m_assembler.getCoeff(C);
    // auto dc = m_assembler.getCoeff(DC);

    // Set the discretization space
    auto w = m_assembler.trialSpace(0);

    m_assembler.assemble(
                m_Gc / m_l0 *
                (
                w*w.tr() +
                1./2. * math::pow(m_l0,2) *igrad(w,G) * igrad(c,G).tr() +
                1./16.* math::pow(m_l0,4) * ilapl(w,G) * ilapl(c,G).tr()
                ) * meas(G)
                );
}

template <class T, enum PForder order, enum PFmode mode>
void gsPhaseFieldAssembler<T,order,mode>::assemblePsiVector(const gsFunctionSet<T> & PSI)
{
    GISMO_ENSURE(m_initialized,"The assembler has not been initialized yet. Call initialize() before assembling the system.");
    m_assembler.clearRhs(); // Resets to zero the values of the already allocated to matrix (LHS)
    m_assembler.initSystem();

    // Set the geometry map
    geometryMap G = m_assembler.getMap(m_patches);

    auto psi= m_assembler.getCoeff(PSI);

    // Set the discretization space
    auto w = m_assembler.trialSpace(0);

    m_assembler.assemble(
                2.0 * psi.val() * w * meas(G)
                );
}


template <class T, enum PForder order, enum PFmode mode>
void gsPhaseFieldAssembler<T,order,mode>::assemblePsiMatrix(const gsFunctionSet<T> & PSI)
{
    GISMO_ENSURE(m_initialized,"The assembler has not been initialized yet. Call initialize() before assembling the system.");
    m_assembler.clearMatrix(); // Resets to zero the values of the already allocated to matrix (LHS)
    m_assembler.initSystem();
    // Set the geometry map
    geometryMap G = m_assembler.getMap(m_patches);

    // auto dc = m_assembler.getCoeff(DC);
    auto psi= m_assembler.getCoeff(PSI);

    // Set the discretization space
    auto w = m_assembler.trialSpace(0);

    m_assembler.assemble(
                2.0 * psi.val() * w * w.tr() * meas(G)
                );
}

// template <class T, enum PForder order, enum PFmode mode>
// void gsPhaseFieldAssembler<T,order,mode>::assembleMassMatrix()
// {
//     GISMO_ENSURE(m_initialized,"The assembler has not been initialized yet. Call initialize() before assembling the system.");
//     m_assembler.clearMatrix(); // Resets to zero the values of the already allocated to matrix (LHS)

//     // Set the geometry map
//     geometryMap G = m_assembler.getMap(m_patches);

//     // Set the discretization space
//     auto w = m_assembler.trialSpace(0);

//     // Initialize the system
//     m_assembler.clearMatrix(); // Resets to zero the values of the already allocated to matrix (LHS)
//     m_assembler.assemble(w*w.tr()*meas(G));// K_m
// }

template <class T, enum PForder order, enum PFmode mode>
void gsPhaseFieldAssembler<T,order,mode>::assembleMatrix()
{
    _assembleMatrix_impl<order,mode>();
}

template <class T, enum PForder order, enum PFmode mode>
template <enum PForder _order, enum PFmode _mode>
typename std::enable_if<(_order==PForder::Second && _mode==PFmode::AT1), void>::type
gsPhaseFieldAssembler<T,order,mode>::_assembleMatrix_impl()
{
    GISMO_ENSURE(m_initialized,"The assembler has not been initialized yet. Call initialize() before assembling the system.");
    m_assembler.clearMatrix(); // Resets to zero the values of the already allocated to matrix (LHS)
    m_assembler.initSystem();

    // Set the geometry map
    geometryMap G = m_assembler.getMap(m_patches);

    // Set the discretization space
    auto w = m_assembler.trialSpace(0);

    m_assembler.assemble( 3./4. * m_Gc * ( igrad(w,G) * igrad(w,G).tr() * m_l0 ) * meas(G) );
}

template <class T, enum PForder order, enum PFmode mode>
template <enum PForder _order, enum PFmode _mode>
typename std::enable_if<(_order==PForder::Fourth && _mode==PFmode::AT1), void>::type
gsPhaseFieldAssembler<T,order,mode>::_assembleMatrix_impl()
{
    GISMO_ENSURE(m_initialized,"The assembler has not been initialized yet. Call initialize() before assembling the system.");
    m_assembler.clearMatrix(); // Resets to zero the values of the already allocated to matrix (LHS)
    m_assembler.initSystem();

    // Set the geometry map
    geometryMap G = m_assembler.getMap(m_patches);

    // Set the discretization space
    auto w = m_assembler.trialSpace(0);

    m_assembler.assemble(
                (m_Gc / m_cw) *
                (
                2. *igrad(w,G) * igrad(w,G).tr() * m_l0 +
                2. * ilapl(w,G) * ilapl(w,G).tr()* math::pow(m_l0,3)
                ) * meas(G));
}

template <class T, enum PForder order, enum PFmode mode>
template <enum PForder _order, enum PFmode _mode>
typename std::enable_if<(_order==PForder::Second && _mode==PFmode::AT2), void>::type
gsPhaseFieldAssembler<T,order,mode>::_assembleMatrix_impl()
{
    GISMO_ENSURE(m_initialized,"The assembler has not been initialized yet. Call initialize() before assembling the system.");
    m_assembler.clearMatrix(); // Resets to zero the values of the already allocated to matrix (LHS)
    m_assembler.initSystem();

    // Set the geometry map
    geometryMap G = m_assembler.getMap(m_patches);

    // Set the discretization space
    auto w = m_assembler.trialSpace(0);

    m_assembler.assemble(
                m_Gc *
                (
                w*w.tr() / m_l0 +
                igrad(w,G) * igrad(w,G).tr() * m_l0
                ) * meas(G)
                );
}

template <class T, enum PForder order, enum PFmode mode>
template <enum PForder _order, enum PFmode _mode>
typename std::enable_if<(_order==PForder::Fourth && _mode==PFmode::AT2), void>::type
gsPhaseFieldAssembler<T,order,mode>::_assembleMatrix_impl()
{
    GISMO_ENSURE(m_initialized,"The assembler has not been initialized yet. Call initialize() before assembling the system.");
    m_assembler.clearMatrix(); // Resets to zero the values of the already allocated to matrix (LHS)
    m_assembler.initSystem();

    // Set the geometry map
    geometryMap G = m_assembler.getMap(m_patches);

    // Set the discretization space
    auto w = m_assembler.trialSpace(0);

    m_assembler.clearMatrix(); // Resets to zero the values of the already allocated to matrix (LHS)

    m_assembler.assemble(
                m_Gc *
                (
                w*w.tr() / m_l0 +
                1./2. *igrad(w,G) * igrad(w,G).tr() * m_l0 +
                1./16.* ilapl(w,G) * ilapl(w,G).tr()* math::pow(m_l0,3)
                ) * meas(G));
}

template <class T, enum PForder order, enum PFmode mode>
void gsPhaseFieldAssembler<T,order,mode>::assembleVector()
{
    _assembleVector_impl<order,mode>();
}

template <class T, enum PForder order, enum PFmode mode>
template <enum PForder _order, enum PFmode _mode>
typename std::enable_if<(_order==PForder::Second && _mode==PFmode::AT1), void>::type
gsPhaseFieldAssembler<T,order,mode>::_assembleVector_impl()
{
    GISMO_ENSURE(m_initialized,"The assembler has not been initialized yet. Call initialize() before assembling the system.");
    m_assembler.clearRhs(); // Resets to zero the values of the already allocated to matrix (LHS)
    m_assembler.initSystem();

    // Set the geometry map
    geometryMap G = m_assembler.getMap(m_patches);

    // Set the discretization space
    auto w = m_assembler.trialSpace(0);

    m_assembler.assemble( 3./8. * m_Gc  * w / m_l0  * meas(G) );
}

template <class T, enum PForder order, enum PFmode mode>
template <enum PForder _order, enum PFmode _mode>
typename std::enable_if<(_order==PForder::Fourth && _mode==PFmode::AT1), void>::type
gsPhaseFieldAssembler<T,order,mode>::_assembleVector_impl()
{
    GISMO_ENSURE(m_initialized,"The assembler has not been initialized yet. Call initialize() before assembling the system.");
    m_assembler.clearRhs(); // Resets to zero the values of the already allocated to matrix (LHS)
    m_assembler.initSystem();

    // Set the geometry map
    geometryMap G = m_assembler.getMap(m_patches);

    // Set the discretization space
    auto w = m_assembler.trialSpace(0);

    m_assembler.assemble( m_Gc / m_cw  * w / m_l0  * meas(G) );
}

template <class T, enum PForder order, enum PFmode mode>
template <enum PForder _order, enum PFmode _mode>
typename std::enable_if<(_order==PForder::Second && _mode==PFmode::AT2), void>::type
gsPhaseFieldAssembler<T,order,mode>::_assembleVector_impl()
{
    GISMO_ENSURE(m_initialized,"The assembler has not been initialized yet. Call initialize() before assembling the system.");
    m_assembler.clearRhs(); // Resets to zero the values of the already allocated to matrix (LHS)
    m_assembler.initVector();
}

template <class T, enum PForder order, enum PFmode mode>
template <enum PForder _order, enum PFmode _mode>
typename std::enable_if<(_order==PForder::Fourth && _mode==PFmode::AT2), void>::type
gsPhaseFieldAssembler<T,order,mode>::_assembleVector_impl()
{
    GISMO_ENSURE(m_initialized,"The assembler has not been initialized yet. Call initialize() before assembling the system.");
    m_assembler.clearRhs(); // Resets to zero the values of the already allocated to matrix (LHS)
    m_assembler.initVector();
}


template <class T, enum PForder order, enum PFmode mode>
void gsPhaseFieldAssembler<T,order,mode>::constructSolution(gsMatrix<T>     & Cvec,
                                                   gsMultiPatch<T> & C) const
{
    auto w  = m_assembler.trialSpace(0);
    auto c  = m_assembler.getSolution(w,  Cvec);
    c.extract(C);
}

template <class T, enum PForder order, enum PFmode mode>
void gsPhaseFieldAssembler<T,order,mode>::constructSolution(gsMatrix<T>         & Cvec,
                                                   gsMappedSpline<2,T> & C) const
{
    auto w  = m_assembler.trialSpace(0);
    auto c  = m_assembler.getSolution(w,  Cvec);
    c.extract(C);
}

template <class T, enum PForder order, enum PFmode mode>
void gsPhaseFieldAssembler<T,order,mode>::constructSolution(const gsMultiPatch<T> & C,
                                                         gsMatrix<T>     & Cvec) const
{
    GISMO_ASSERT(C.geoDim()==1,"C must be a scalar function");
    GISMO_ASSERT(C.nPatches()==m_basis.nBases(),"Number of patches in C must be equal to the number of bases in the assembler");
    auto w  = m_assembler.trialSpace(0);

    Cvec.setZero(this->numDofs(),1);
    for (size_t b=0; b!=m_basis.nBases(); b++)
    {
        for (index_t i = 0; i < m_basis.basis(b).size(); i++)
        {
            GISMO_ASSERT(C.basis(b).size()==m_basis.basis(b).size(),"Number of basis functions in C must be equal to the number of basis functions in the assembler");
            if (w.mapper().is_free(i,b))
                Cvec(w.mapper().index(i,b)) = C.patch(b).coefs()(i,0);
        }
    }
}

template <class T, enum PForder order, enum PFmode mode>
void gsPhaseFieldAssembler<T,order,mode>::setSpaceBasis(const gsFunctionSet<T> & spaceBasis)
{
    m_spaceBasis = &spaceBasis;
}

}// namespace gismo
