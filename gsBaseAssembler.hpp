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

    const index_t dim = unknowns.rows();
    std::vector<gsDofMapper> mappers(dim);
    for (index_t unk = 0; unk < dim; ++unk)
        mappers[unk] = m_system.colMapper(unknowns[unk]);

    gsVector<index_t> basisIndices(dim);
    for (index_t unk = 0; unk < dim; ++unk)
        basisIndices[unk] = m_system.colBasis(unknowns[unk]);

    for (size_t p = 0; p < m_pde_ptr->domain().nPatches(); ++p)
    {
        const int size  = m_bases[basisIndices[0]][p].size();
        gsMatrix<T> coeffs(size,dim);
        for (index_t unk = 0; unk < dim; ++unk)
            for (index_t i = 0; i < size; ++i)
                if (mappers[unk].is_free(i, p) ) // DoF value is in the solVector
                {
                    unsigned idx;
                    m_system.mapToGlobalColIndex(i,p,idx,unknowns[unk]);
                    coeffs(i,unk) = solVector(idx,0);
                }
                else // eliminated DoF: fill with Dirichlet data
                    coeffs(i,unk) = fixedDoFs[unknowns[unk]](mappers[unk].bindex(i, p),0);

        result.addPatch(m_bases[basisIndices[0]][p].makeGeometry(give(coeffs)));
    }
}

//--------------------- DIRICHLET BC SHENANIGANS ----------------------------------//


template <class T>
void gsBaseAssembler<T>::setFixedDofs(size_t patch, boxSide side, const gsMatrix<T> & ddofs)
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
    GISMO_ENSURE(localBIndices.rows() == ddofs.rows() && m_dim == ddofs.cols(),
                 "Wrong size of a given matrix with Dirichlet DoFs: " + util::to_string(ddofs.rows()) +
                 " x " + util::to_string(ddofs.cols()) + ". Must be:" + util::to_string(localBIndices.rows()) +
                 " x " + util::to_string(m_dim));

    for (short_t d = 0; d < m_dim; ++d )
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

}// namespace gismo ends
