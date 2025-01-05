/** @file gsVisitorMixedLinearElasticity.h

    @brief Visitor class for volumetric integration of the mixed linear elasticity system.

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s):
        O. Weeger    (2012 - 2015, TU Kaiserslautern),
        A.Shamanskiy (2016 - ...., TU Kaiserslautern)
*/

#pragma once

#include <gsElasticity/gsVisitorElUtils.h>
#include <gsElasticity/gsBasePde.h>

#include <gsAssembler/gsQuadrature.h>
#include <gsCore/gsFuncData.h>

namespace gismo
{

template <class T>
class gsVisitorMixedLinearElasticity
{
public:
    gsVisitorMixedLinearElasticity(const gsPde<T> & pde_)
    : dim(0), N_D(0), pde_ptr(static_cast<const gsBasePde<T>*>(&pde_)) {}

    void initialize(const gsBasisRefs<T> & basisRefs,
                    const index_t patchIndex,
                    const gsOptionList & options,
                    gsQuadRule<T> & rule)
    {
        GISMO_UNUSED(patchIndex);
        // parametric dimension of the first displacement component
        dim = basisRefs.front().dim();
        // a quadrature rule is defined by the basis for the first velocity component.
        // the same rule is used for the presure
        rule = gsQuadrature::get(basisRefs.front(), options);
        // saving necessary info
        T E = options.getReal("YoungsModulus");
        T pr = options.getReal("PoissonsRatio");
        lambda_inv = ( 1. + pr ) * ( 1. - 2. * pr ) / E / pr ;
        mu     = E / ( 2. * ( 1. + pr ) );
        forceScaling = options.getReal("ForceScaling");
        I = gsMatrix<T>::Identity(dim,dim);
        // resize containers for global indices
        globalIndices.resize(dim+1);
        blockNumbers.resize(dim+1);
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
        // find local indices of the displacement and pressure basis functions active on the element
        basisRefs.front().active_into(quNodes.col(0),localIndicesDisp);
        N_D = localIndicesDisp.rows();
        basisRefs.back().active_into(quNodes.col(0), localIndicesPres);
        N_P = localIndicesPres.rows();
        // Evaluate displacement basis functions and their derivatives on the element
        basisRefs.front().evalAllDers_into(quNodes,1,basisValuesDisp);
        // Evaluate pressure basis functions on the element
        basisRefs.back().eval_into(quNodes,basisValuesPres);
        // Evaluate right-hand side at the image of the quadrature points
        pde_ptr->rhs()->eval_into(md.values[0],forceValues);
    }

    inline void assemble(gsDomainIterator<T> & element,
                         const gsVector<T> & quWeights)
    {
        GISMO_UNUSED(element);
        // Initialize local matrix/rhs                      // A | B^T
        localMat.setZero(dim*N_D + N_P, dim*N_D + N_P);     // --|--    matrix structure
        localRhs.setZero(dim*N_D + N_P,1);                  // B | C
        // elasticity tensor
        symmetricIdentityTensor<T>(C,I);
        C *= mu;
        // Loop over the quadrature nodes
        for (index_t q = 0; q < quWeights.rows(); ++q)
        {
            // Multiply quadrature weight by the geometry measure
            const T weight = quWeights[q] * md.measure(q);

            // Compute physical gradients of basis functions at q as a dim x numActiveFunction matrix
            transformGradients(md, q, basisValuesDisp[1], physGradDisp);
            // A-matrix: Loop over displacement basis functions
            for (index_t i = 0; i < N_D; i++)
            {
                setB<T>(B_i,I,physGradDisp.col(i));
                tempK = B_i.transpose() * C;
                // Loop over displacement basis functions
                for (index_t j = 0; j < N_D; j++)
                {
                    setB<T>(B_j,I,physGradDisp.col(j));
                    K = tempK * B_j;

                    for (short_t di = 0; di < dim; ++di)
                        for (short_t dj = 0; dj < dim; ++dj)
                            localMat(di*N_D+i, dj*N_D+j) += weight * K(di,dj);
                }
            }
            // B-matrix
            for (short_t d = 0; d < dim; ++d)
            {
                block = weight*basisValuesPres.col(q)*physGradDisp.row(d);
                localMat.block(dim*N_D,d*N_D,N_P,N_D) += block.block(0,0,N_P,N_D);
                localMat.block(d*N_D,dim*N_D,N_D,N_P) += block.transpose().block(0,0,N_D,N_P);
            }

            // C-matrix
            if (abs(lambda_inv) > 0)
                localMat.block(dim*N_D,dim*N_D,N_P,N_P) -=
                    (weight*lambda_inv*basisValuesPres.col(q)*basisValuesPres.col(q).transpose()).block(0,0,N_P,N_P);

            // rhs contribution
            for (short_t d = 0; d < dim; ++d)
                localRhs.middleRows(d*N_D,N_D).noalias() += weight * forceScaling * forceValues(d,q) * basisValuesDisp[0].col(q) ;
        }
    }

    inline void localToGlobal(const int patchIndex,
                              const std::vector<gsMatrix<T> > & eliminatedDofs,
                              gsSparseSystem<T> & system)
    {
        // computes global indices for displacement components
        for (short_t d = 0; d < dim; ++d)
        {
            system.mapColIndices(localIndicesDisp,patchIndex,globalIndices[d],d);
            blockNumbers.at(d) = d;
        }
        // computes global indices for pressure
        system.mapColIndices(localIndicesPres, patchIndex, globalIndices[dim], dim);
        blockNumbers.at(dim) = dim;
        // push to global system
        system.pushToRhs(localRhs,globalIndices,blockNumbers);
        system.pushToMatrix(localMat,globalIndices,eliminatedDofs,blockNumbers,blockNumbers);
    }

protected:
    // problem info
    short_t dim;
    const gsBasePde<T> * pde_ptr;
    // Lame coefficients and force scaling factor
    T lambda_inv, mu, forceScaling;
    // geometry mapping
    gsMapData<T> md;
    // local components of the global linear system
    gsMatrix<T> localMat;
    gsMatrix<T> localRhs;
    // local indices (at the current patch) of basis functions active at the current element
    gsMatrix<index_t> localIndicesDisp;
    gsMatrix<index_t> localIndicesPres;
    // number of displacement and pressure basis functions active at the current element
    index_t N_D, N_P;
    // values and derivatives of displacement basis functions at quadrature points at the current element
    // values are stored as a N_D x numQuadPoints matrix; not sure about derivatives, must be smth like N_D*dim x numQuadPoints
    std::vector<gsMatrix<T> > basisValuesDisp;
    // values of pressure basis functions active at the current element;
    // stores as a N_P x numQuadPoints matrix
    gsMatrix<T> basisValuesPres;
    // RHS values at quadrature points at the current element; stored as a dim x numQuadPoints matrix
    gsMatrix<T> forceValues;

    // all temporary matrices defined here for efficiency
    gsMatrix<T> C, physGradDisp, B_i, tempK, B_j, K, block, I;
    // containers for global indices
    std::vector< gsMatrix<index_t> > globalIndices;
    gsVector<index_t> blockNumbers;
};

} // namespace gismo
