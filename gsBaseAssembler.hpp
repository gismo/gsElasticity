/** @file gsBaseAssembler.hpp

    @brief Provides implementation for gsBaseAssembler.

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s):
        A.Shamanskiy (2016 - ...., TU Kaiserslautern)
*/

#pragma once

#include <gsElasticity/gsBaseAssembler.h>

namespace gismo
{

template <class T>
void gsBaseAssembler<T>::constructSolution(const gsMatrix<T> & solVector,
                                           const std::vector<gsMatrix<T> > & fixedDoFs,
                                           gsMultiPatch<T> & result,
                                           const gsVector<index_t> & unknowns) const
{
    result.clear();
    GISMO_ENSURE(unknowns.rows() > 0, "No unknowns provided!");
    index_t nRhs = m_system.rhs().cols();
    unsigned idx;
    for (size_t p = 0; p < m_pde_ptr->domain().nPatches(); ++p)
    {
        index_t size  = m_bases[unknowns[0]][p].size();
        // check that all unknowns have the same basis size
        for (index_t i = 1; i < unknowns.rows(); ++i)
        {
            GISMO_ENSURE(m_bases[unknowns[i]][p].size() == size,"Unknowns do not have the same size!");
        }
        gsMatrix<T> coeffs(size,unknowns.rows()*nRhs);
        for (index_t unk = 0; unk < unknowns.rows(); ++unk)
            for (index_t i = 0; i < size; ++i)
                if (m_system.colMapper(unknowns[unk]).is_free(i,p)) // DoF is in the solVector
                {
                    m_system.mapToGlobalColIndex(i,p,idx,unknowns[unk]);
                    coeffs.block(i,unk*nRhs,1,nRhs) = solVector.block(idx,0,1,nRhs);
                }
                else // DoF is eliminated; it is the Dirichlet Dofs
                {
                    idx = m_system.colMapper(unknowns[unk]).bindex(i, p);
                    coeffs.block(i,unk*nRhs,1,nRhs) = fixedDoFs[unknowns[unk]].block(idx,0,1,nRhs);
                }
        result.addPatch(m_bases[unknowns[0]][p].makeGeometry(give(coeffs)));
    }
}

//--------------------- DIRICHLET BC SHENANIGANS ----------------------------------//


template <class T>
void gsBaseAssembler<T>::setFixedDofs(size_t patch, boxSide side, const gsMatrix<T> & ddofs, bool oneUnk)
{
    bool dirBcExists = false;
    typename gsBoundaryConditions<T>::const_iterator it = m_pde_ptr->bc().dirichletBegin();
    while (!dirBcExists && it != m_pde_ptr->bc().dirichletEnd())
    {
        if (it->patch() == patch && it->side() == side)
            dirBcExists = true;
        ++it;
    }
    GISMO_ENSURE(dirBcExists,"Side " + util::to_string(side) + " of patch " + util::to_string(patch)
                             + " does not belong to the Dirichlet boundary.");

    gsMatrix<unsigned> localBIndices = m_bases[0][patch].boundary(side);
    GISMO_ENSURE(localBIndices.rows() == ddofs.rows(),
                 "Wrong size of a given matrix with Dirichlet DoFs: " + util::to_string(ddofs.rows()) +
                 ". Must be:" + util::to_string(localBIndices.rows()));

    if (oneUnk)
    {
        gsMatrix<unsigned> globalIndices;
        m_system.mapColIndices(localBIndices, patch, globalIndices, 0);
        for (index_t i = 0; i < globalIndices.rows(); ++i)
            m_ddof[0].row(m_system.colMapper(0).global_to_bindex(globalIndices(i,0))) = ddofs.row(i);
    }
    else
        for (short_t d = 0; d < ddofs.cols(); ++d )
        {
            gsMatrix<unsigned> globalIndices;
            m_system.mapColIndices(localBIndices, patch, globalIndices, d);
            for (index_t i = 0; i < globalIndices.rows(); ++i)
                m_ddof[d](m_system.colMapper(d).global_to_bindex(globalIndices(i,0)),0) = ddofs(i,d);
        }
}

template <class T>
void gsBaseAssembler<T>::setFixedDofs(const std::vector<gsMatrix<T> > & ddofs)
{
    GISMO_ENSURE(ddofs.size() >= m_ddof.size(), "Wrong size of the container with fixed DoFs: " + util::to_string(ddofs.size()) +
                 ". Must be at least: " + util::to_string(m_ddof.size()));

    for (short_t d = 0; d < index_t(m_ddof.size()); ++d)
    {
        GISMO_ENSURE(m_ddof[d].rows() == ddofs[d].rows(),"Wrong number of fixed DoFs for " + util::to_string(d) + "component: " +
                     util::to_string(ddofs[d].rows()) + ". Must be: " + util::to_string(m_ddof[d].rows()));
        m_ddof[d] = ddofs[d];
    }
}

template <class T>
void gsBaseAssembler<T>::getFixedDofs(size_t patch, boxSide side, gsMatrix<T> & ddofs)
{
    bool dirBcExists = false;
    typename gsBoundaryConditions<T>::const_iterator it = m_pde_ptr->bc().dirichletBegin();
    while (!dirBcExists && it != m_pde_ptr->bc().dirichletEnd())
    {
        if (it->patch() == patch && it->side() == side)
            dirBcExists = true;
        ++it;
    }
    GISMO_ENSURE(dirBcExists,"Side " + util::to_string(side) + " of patch " + util::to_string(patch)
                             + " does not belong to the Dirichlet boundary.");

    short_t m_dim = m_pde_ptr->domain().targetDim();
    gsMatrix<unsigned> localBIndices = m_bases[0][patch].boundary(side);
    ddofs.clear();
    ddofs.resize(localBIndices.rows(),m_dim);
    for (short_t d = 0; d < m_dim; ++d )
    {
        gsMatrix<unsigned> globalIndices;
        m_system.mapColIndices(localBIndices, patch, globalIndices, d);

        for (index_t i = 0; i < globalIndices.rows(); ++i)
            ddofs(i,d) = m_ddof[d](m_system.colMapper(d).global_to_bindex(globalIndices(i,0)),0);
    }

}

}// namespace gismo ends
