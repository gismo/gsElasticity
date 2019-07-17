/** @file gsBaseAssembler.h

    @brief Base class for assemblers of gsElasticity.

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s):
        O. Weeger    (2012 - 2015, TU Kaiserslautern),
        A.Shamanskiy (2016 - ...., TU Kaiserslautern)
*/

#pragma once

#include <gsAssembler/gsAssembler.h>

namespace gismo
{

/** @brief Extends the gsAssembler class by adding functionality necessary for a general nonlinear solver.
 * Potentially, can be merged back into gsAssembler.
 */
template <class T>
class gsBaseAssembler : public gsAssembler<T>
{
public:
    /// assembles the linear system given the current solution vector.
    /// checks if the current solution is valid (Newton's solver can exit safely if invalid).
    /// assembles only the RHS vector if \a assembleMatrix is set to false
    virtual bool assemble(const gsMatrix<T> & solutionVector, bool assembleMatrix = true) = 0;
    /// sets scaling of Dirichlet BC used for linear system assembly
    virtual void setDirichletAssemblyScaling(T factor) = 0;
    /// sets scaling of Dirichlet BC used for construction of the solution as a gsMultiPatch object
    virtual void setDirichletConstructionScaling(T factor) = 0;
    /// set scaling of the force loading (volume and surface loading)
    virtual void setForceScaling(T factor) = 0;

    virtual int numDofs() const { return gsAssembler<T>::numDofs(); }

    /// constructs solution as a gsMultiPatch object from the solution vector.
    /// exactly the same as the gsAssembler<>::constructSolution, but allows to control the scaling of Dirichlet BC
    virtual void constructSolution(const gsMatrix<T>& solVector,
                                   gsMultiPatch<T>& result,
                                   const gsVector<index_t> & unknowns) const
    {
        GISMO_ASSERT(solVector.cols()==1, "Vector valued output only works for single rhs");
        unsigned idx;

        const index_t dim = unknowns.rows();

        std::vector<gsDofMapper> mappers(dim);
        for(index_t unk = 0; unk<dim;++unk)
            mappers[unk] = gsAssembler<T>::m_system.colMapper(unknowns[unk]);

        result.clear(); // result is cleared first

        gsMatrix<T> coeffs;
        gsVector<index_t> basisIndices(dim);
        for(index_t unk = 0; unk<dim;++unk)
            basisIndices[unk] = gsAssembler<T>::m_system.colBasis(unknowns[unk]);

        for (size_t p = 0; p < gsAssembler<T>::m_pde_ptr->domain().nPatches(); ++p)
        {
            const int sz  = gsAssembler<T>::m_bases[basisIndices[0]][p].size(); //must be equal for all unk
            coeffs.resize(sz, dim);

            for(index_t unk = 0; unk<dim;++unk)
                for (index_t i = 0; i < sz; ++i)
                    if ( mappers[unk].is_free(i, p) ) // DoF value is in the solVector
                    {
                        gsAssembler<T>::m_system.mapToGlobalColIndex(i,p,idx,unknowns[unk]);
                        coeffs(i,unk) = solVector(idx,0);
                    }
                    else // eliminated DoF: fill with Dirichlet data
                        coeffs(i,unk) = gsAssembler<T>::m_ddof[unknowns[unk]](mappers[unk].bindex(i, p),0) *
                                        gsAssembler<T>::m_options.getReal("DirichletConstruction");

            result.addPatch( gsAssembler<T>::m_bases[basisIndices[0]][p].makeGeometry( give(coeffs) ) );
        }
    }
};

} // namespace ends
