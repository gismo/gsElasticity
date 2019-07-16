/** @file gsVisitorElMass.h

    @brief Visitor class for the mass matrix assembly for elasticity problems.

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s):
        O. Weeger    (2012 - 2015, TU Kaiserslautern),
        A.Shamanskiy (2016 - ...., TU Kaiserslautern)
*/

#pragma once

#include <gsAssembler/gsQuadrature.h>
#include <gsCore/gsFuncData.h>

namespace gismo
{

template <class T>
class gsVisitorElMass
{
public:
    gsVisitorElMass(const gsPde<T> & pde_) {}

    void initialize(const gsBasisRefs<T> & basisRefs,
                    const index_t patchIndex,
                    const gsOptionList & options,
                    gsQuadRule<T> & rule)
    {
        // parametric dimension of the first displacement component
        dim = basisRefs.front().dim();
        // a quadrature rule is defined by the basis for the first displacement component.
        rule = gsQuadrature::get(basisRefs.front(), options);
        // saving necessary info
        density = options.getReal("Density");
    }

    inline void evaluate(const gsBasisRefs<T> & basisRefs,
                         const gsGeometry<T> & geo,
                         const gsMatrix<T> & quNodes)
    {
        // store quadrature points of the element for geometry evaluation
        md.points = quNodes;
        // NEED_MEASURE to get the Jacobian determinant values for integration
        md.flags = NEED_MEASURE;
        // Compute the geometry mapping at the quadrature points
        geo.computeMap(md);
        // find local indices of the displacement basis functions active on the element
        basisRefs.front().active_into(quNodes.col(0),localIndicesDisp);
        N_D = localIndicesDisp.rows();
        // Evaluate displacement basis functions on the element
        basisRefs.front().eval_into(quNodes,basisValuesDisp);
    }

    inline void assemble(gsDomainIterator<T> & element,
                         const gsVector<T> & quWeights)
    {
        // initialize local matrix and rhs
        localMat.setZero(dim*N_D,dim*N_D);
        for (index_t q = 0; q < quWeights.rows(); ++q)
        {
            // Multiply quadrature weight by the geometry measure
            const T weight = density * quWeights[q] * md.measure(q);

            gsMatrix<T> block = weight * basisValuesDisp.col(q) * basisValuesDisp.col(q).transpose();
            for (short_t d = 0; d < dim; ++d)
                localMat.block(d*N_D,d*N_D,N_D,N_D) += block.block(0,0,N_D,N_D);
        }
    }

    inline void localToGlobal(const int patchIndex,
                              const std::vector<gsMatrix<T> > & eliminatedDofs,
                              gsSparseSystem<T> & system)
    {
        // number of unknowns: dim of displacement
        std::vector< gsMatrix<unsigned> > globalIndices(dim);
        gsVector<size_t> blockNumbers(dim);
        // computes global indices for displacement components
        for (short_t d = 0; d < dim; ++d)
        {
            system.mapColIndices(localIndicesDisp, patchIndex, globalIndices[d], d);
            blockNumbers.at(d) = d;
        }
        // push to global system
        system.pushToMatrix(localMat,globalIndices,eliminatedDofs,blockNumbers,blockNumbers);
    }

protected:
    // problem info
    short_t dim;
    //density
    T density;
    // geometry mapping
    gsMapData<T> md;
    // local components of the global linear system
    gsMatrix<T> localMat;
    // local indices (at the current patch) of the displacement basis functions active at the current element
    gsMatrix<unsigned> localIndicesDisp;
    // number of displacement basis functions active at the current element
    index_t N_D;
    // values of displacement basis functions at quadrature points at the current element stored as a N_D x numQuadPoints matrix;
    gsMatrix<T> basisValuesDisp;
};

} // namespace gismo

