/** @file gsVisitorElasticityNeumann.h

    @brief Visitor class for the surface load integration.

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
class gsVisitorElasticityNeumann
{
public:

    gsVisitorElasticityNeumann(const gsPde<T> & pde_,
                               const boundary_condition<T> & s)
        : neumannFunction_ptr(s.function().get()),
          patchSide(s.side()) {}

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
        forceScaling = options.getReal("ForceScaling");
        // resize containers for global indices
        globalIndices.resize(dim);
        blockNumbers.resize(dim);
    }

    inline void evaluate(const gsBasisRefs<T> & basisRefs,
                         const gsGeometry<T> & geo,
                         const gsMatrix<T> & quNodes)
    {
        // store quadrature points of the element for geometry evaluation
        md.points = quNodes;
        // NEED_VALUE to get points in the physical domain for evaluation of the load
        // NEED_MEASURE to get the Jacobian determinant values for integration
        md.flags = NEED_VALUE | NEED_MEASURE;
        // Compute image of the quadrature points plus gradient, jacobian and other necessary data
        geo.computeMap(md);
        // Evaluate the Neumann functon on the images of the quadrature points
        neumannFunction_ptr->eval_into(md.values[0], neumannValues);
        // find local indices of the displacement basis functions active on the element
        basisRefs.front().active_into(quNodes.col(0),localIndicesDisp);
        N_D = localIndicesDisp.rows();
        // Evaluate basis functions on element
        basisRefs.front().eval_into(quNodes,basisValuesDisp);
    }

    inline void assemble(gsDomainIterator<T> & element,
                         const gsVector<T> & quWeights)
    {
        // Initialize local matrix/rhs
        localRhs.setZero(dim*N_D,1);
        // loop over the quadrature nodes
        for (index_t q = 0; q < quWeights.rows(); ++q)
        {
			// Compute the outer normal vector on the side
            // normal length equals to the local area measure
            outerNormal(md, q, patchSide, unormal);
            // Collect the factors here: quadrature weight and geometry measure
            const T weight = quWeights[q] * unormal.norm() * forceScaling;

            for (short_t d = 0; d < dim; ++d)
                localRhs.middleRows(d*N_D,N_D).noalias() += weight * neumannValues(d,q) * basisValuesDisp.col(q);
        }
    }

    inline void localToGlobal(const int patchIndex,
                                  const std::vector<gsMatrix<T> > & eliminatedDofs,
                                  gsSparseSystem<T> & system)
        {
            // computes global indices for displacement components
            for (short_t d = 0; d < dim; ++d)
            {
                system.mapColIndices(localIndicesDisp, patchIndex, globalIndices[d], d);
                blockNumbers.at(d) = d;
            }
            // push to global system
            system.pushToRhs(localRhs,globalIndices,blockNumbers);
        }

protected:   
    // problem info
    short_t dim;
    const gsFunction<T> * neumannFunction_ptr;
    T forceScaling;
    boxSide patchSide;
    // geometry mapping
    gsMapData<T> md;
    // local components of the global linear system
    gsMatrix<T> localRhs;
    // local indices (at the current patch) of the displacement basis functions active at the current element
    gsMatrix<unsigned> localIndicesDisp;
    // number of displacement basis functions active at the current element
    index_t N_D;
    // values of displacement basis functions at quadrature points at the current element stored as a N_D x numQuadPoints matrix;
    gsMatrix<T> basisValuesDisp;
    // values of the boundary loading function stored as a dim x numQuadPoints matrix;
    gsMatrix<T> neumannValues;

    // all temporary matrices defined here for efficiency
    gsVector<T> unormal;
    // containers for global indices
    std::vector< gsMatrix<unsigned> > globalIndices;
    gsVector<size_t> blockNumbers;
};

} // namespace gismo
