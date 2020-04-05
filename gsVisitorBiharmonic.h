/** @file gsVisitorBiharmonic.h

    @brief Visitor class for mixed formulation for the biharmonic equation

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
class gsVisitorBiharmonic
{
public:

    gsVisitorBiharmonic(const gsPde<T> & pde_)
        : pde_ptr(static_cast<const gsPoissonPde<T>*>(&pde_))
    {}

    void initialize(const gsBasisRefs<T> & basisRefs,
                    const index_t patchIndex,
                    const gsOptionList & options,
                    gsQuadRule<T> & rule)
    {
        // a quadrature rule is defined by the basis for the auxiliary variable.
        // the same rule is used for the main variable
        rule = gsQuadrature::get(basisRefs.back(), options);
        // saving necessary info
        localStiffening = options.getReal("LocalStiff");
    }

    inline void evaluate(const gsBasisRefs<T> & basisRefs,
                         const gsGeometry<T> & geo,
                         const gsMatrix<T> & quNodes)
    {
        // store quadrature points of the element for geometry evaluation
        md.points = quNodes;
        // NEED_VALUE to get points in the physical domain for evaluation of the RHS
        // NEED_MEASURE to get the Jacobian determinant values for integration
        // NEED_GRAD_TRANSFORM to get the Jacobian matrix to transform gradient from the parametric to physical domain
        md.flags = NEED_VALUE | NEED_MEASURE | NEED_GRAD_TRANSFORM;
        // Compute image of the quadrature points plus gradient, jacobian and other necessary data
        geo.computeMap(md);
        // find local indices of the main and auxiliary basis functions active on the element
        basisRefs.front().active_into(quNodes.col(0),localIndicesMain);
        N_M = localIndicesMain.rows();
        basisRefs.back().active_into(quNodes.col(0), localIndicesAux);
        N_A = localIndicesAux.rows();
        // Evaluate basis functions and their derivatives on the element
        basisRefs.front().evalAllDers_into(quNodes,1,basisValuesMain);
        basisRefs.back().evalAllDers_into(quNodes,1,basisValuesAux);
        // Evaluate right-hand side at the image of the quadrature points
        pde_ptr->rhs()->eval_into(md.values[0],forceValues);
    }

    inline void assemble(gsDomainIterator<T> & element,
                         const gsVector<T> & quWeights)
    {
        // Initialize local matrix/rhs                  // 0 | B^T = L
        localMat.setZero(N_M + N_A, N_M + N_A);         // --|--    matrix structure
        localRhs.setZero(N_M + N_A,pde_ptr->numRhs());  // B | A   = 0

        // Loop over the quadrature nodes
        for (index_t q = 0; q < quWeights.rows(); ++q)
        {
            // Multiply quadrature weight by the geometry measure
            const T weightMatrix = quWeights[q] * pow(md.measure(q),1-localStiffening);
            const T weightRHS = quWeights[q] * md.measure(q);
            // Compute physical gradients of the basis functions at q as a 1 x numActiveFunction matrix
            transformGradients(md, q, basisValuesAux[1], physGradAux);
            transformGradients(md, q, basisValuesMain[1], physGradMain);
            // matrix A
            block = weightMatrix * basisValuesAux[0].col(q) * basisValuesAux[0].col(q).transpose();
            localMat.block(N_M,N_M,N_A,N_A) += block.block(0,0,N_A,N_A);
            // matrix B
            block = weightMatrix * physGradAux.transpose()*physGradMain; // N_A x N_M
            localMat.block(0,N_M,N_M,N_A) += block.block(0,0,N_A,N_M).transpose();
            localMat.block(N_M,0,N_A,N_M) += block.block(0,0,N_A,N_M);
            // rhs contribution
            localRhs.middleRows(0,N_M).noalias() += weightRHS * basisValuesMain[0].col(q) * forceValues.col(q).transpose();
        }
    }

    inline void localToGlobal(const int patchIndex,
                              const std::vector<gsMatrix<T> > & eliminatedDofs,
                              gsSparseSystem<T> & system)
    {
        // number of unknowns: 2
        std::vector< gsMatrix<unsigned> > globalIndices(2);
        gsVector<size_t> blockNumbers(2);
        // compute global indices
        system.mapColIndices(localIndicesMain,patchIndex,globalIndices[0],0);
        system.mapColIndices(localIndicesAux,patchIndex,globalIndices[1],1);
        blockNumbers.at(0) = 0;
        blockNumbers.at(1) = 1;
        // push to global system
        system.pushToRhs(localRhs,globalIndices,blockNumbers);
        system.pushToMatrix(localMat,globalIndices,eliminatedDofs,blockNumbers,blockNumbers);
    }

protected:
    // problem info
    const gsPoissonPde<T> * pde_ptr;
    // geometry mapping
    gsMapData<T> md;
    // local components of the global linear system
    gsMatrix<T> localMat;
    gsMatrix<T> localRhs;
    // local indices (at the current patch) of basis functions active at the current element
    gsMatrix<unsigned> localIndicesMain;
    gsMatrix<unsigned> localIndicesAux;
    // number of main and auxiliary basis functions active at the current element
    index_t N_M, N_A;
    // values and derivatives of main basis functions at quadrature points at the current element
    // values are stored as a N_M x numQuadPoints matrix; not sure about derivatives, must be smth like N_M x numQuadPoints
    // same for the auxiliary basis functions
    std::vector<gsMatrix<T> > basisValuesMain;
    std::vector<gsMatrix<T> > basisValuesAux;
    // RHS values at quadrature points at the current element; stored as a 1 x numQuadPoints matrix
    gsMatrix<T> forceValues;
    // all temporary matrices defined here for efficiency
    gsMatrix<T> block, physGradMain, physGradAux;
    real_t localStiffening;
};

} // namespace gismo
