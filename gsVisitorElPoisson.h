/** @file gsVisitorElPoisson.h

    @brief Visitor class for Poisson's equation.

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s):
        A.Shamanskiy (2016 - ...., TU Kaiserslautern)
*/

#pragma once

#include <gsAssembler/gsQuadrature.h>
#include <gsCore/gsFuncData.h>

namespace gismo
{

template <class T>
class gsVisitorElPoisson
{
public:

    gsVisitorElPoisson(const gsPde<T> & pde_, gsSparseMatrix<T> * elimMatrix = nullptr)
        : pde_ptr(static_cast<const gsPoissonPde<T>*>(&pde_)),
          elimMat(elimMatrix)
    {}

    void initialize(const gsBasisRefs<T> & basisRefs,
                    const index_t patchIndex,
                    const gsOptionList & options,
                    gsQuadRule<T> & rule)
    {
        // a quadrature rule is defined by the basis for the first displacement component.
        rule = gsQuadrature::get(basisRefs.front(), options);
        // saving necessary info
        localStiffening = options.getReal("LocalStiff");
        // resize containers for global indices
        globalIndices.resize(1);
        blockNumbers.resize(1);
    }

    inline void evaluate(const gsBasisRefs<T> & basisRefs,
                         const gsGeometry<T> & geo,
                         const gsMatrix<T> & quNodes)
    {
        // store quadrature points of the element for geometry evaluation
        md.points = quNodes;
        // NEED_MEASURE to get the Jacobian determinant values for integration
        md.flags = NEED_MEASURE | NEED_GRAD_TRANSFORM | NEED_VALUE;
        // Compute the geometry mapping at the quadrature points
        geo.computeMap(md);
        // Evaluate displacement basis functions on the element
        basisRefs.front().evalAllDers_into(quNodes,1,basisValues);
        // find local indices of the displacement basis functions active on the element
        basisRefs.front().active_into(quNodes.col(0),localIndices);
        N = localIndices.rows();
        pde_ptr->rhs()->eval_into(md.values[0],forceValues);
    }

    inline void assemble(gsDomainIterator<T> & element,
                         const gsVector<T> & quWeights)
    {
        // initialize local matrix and rhs
        localMat.setZero(N,N);
        localRhs.setZero(N,pde_ptr->numRhs());
        for (index_t q = 0; q < quWeights.rows(); ++q)
        {
            // Multiply quadrature weight by the geometry measure
            const T weightMatrix = quWeights[q] * pow(md.measure(q),1-localStiffening);
            const T weightRHS = quWeights[q] * md.measure(q);
            transformGradients(md,q,basisValues[1],physGrad);
            localMat.noalias() += weightMatrix * (physGrad.transpose() * physGrad);
            localRhs.noalias() += weightRHS * basisValues[0].col(q) * forceValues.col(q).transpose();
        }
    }

    inline void localToGlobal(const int patchIndex,
                              const std::vector<gsMatrix<T> > & eliminatedDofs,
                              gsSparseSystem<T> & system)
    {
        // computes global indices for displacement components
        system.mapColIndices(localIndices, patchIndex, globalIndices[0], 0);
        blockNumbers.at(0) = 0;
        // push to global system
        system.pushToMatrix(localMat,globalIndices,eliminatedDofs,blockNumbers,blockNumbers);
        system.pushToRhs(localRhs,globalIndices,blockNumbers);

        // push to the elimination matrix
        if (elimMat != nullptr)
        {
            index_t globalI, globalElimJ;
            for (index_t i = 0; i < N; ++i)
                if (system.colMapper(0).is_free_index(globalIndices[0].at(i)))
                {
                    system.mapToGlobalRowIndex(localIndices.at(i),patchIndex,globalI,0);
                    for (index_t j = 0; j < N; ++j)
                        if (!system.colMapper(0).is_free_index(globalIndices[0].at(j)))
                        {
                            globalElimJ = system.colMapper(0).global_to_bindex(globalIndices[0].at(j));
                            elimMat->coeffRef(globalI,globalElimJ) += localMat(i,j);
                        }
                }
        }
    }

protected:
    // geometry mapping
    gsMapData<T> md;
    const gsPoissonPde<T> * pde_ptr;
    // local components of the global linear system
    gsMatrix<T> localMat;
    gsMatrix<T> localRhs;
    // local indices (at the current patch) of the displacement basis functions active at the current element
    gsMatrix<index_t> localIndices;
    // number of displacement basis functions active at the current element
    index_t N;
    // values of displacement basis functions at quadrature points at the current element stored as a N_D x numQuadPoints matrix;
    std::vector<gsMatrix<T> >basisValues;
    // RHS values at quadrature points at the current element; stored as a dim x numQuadPoints matrix
    gsMatrix<T> forceValues;
    // elimination matrix to efficiently change Dirichlet degrees of freedom
    gsSparseMatrix<T> * elimMat;

    // all temporary matrices defined here for efficiency
    gsMatrix<T> physGrad;
    real_t localStiffening;
    // containers for global indices
    std::vector< gsMatrix<index_t> > globalIndices;
    gsVector<index_t> blockNumbers;
};

} // namespace gismo

