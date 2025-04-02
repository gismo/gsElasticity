/** @file gsVisitorLinearElasticity.h

    @brief Visitor class for volumetric integration of the linear elasticity system.

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
class gsVisitorLinearElasticity
{
public:

    gsVisitorLinearElasticity(const gsPde<T> & pde_, gsSparseMatrix<T> * elimMatrix = nullptr)
    :
    dim(0),
    pde_ptr(static_cast<const gsBasePde<T>*>(&pde_)),
    N_D(0),
    elimMat(elimMatrix)
    {}

    void initialize(const gsBasisRefs<T> & basisRefs,
                    const index_t patchIndex,
                    const gsOptionList & options,
                    gsQuadRule<T> & rule)
    {
        GISMO_UNUSED(patchIndex);
        // parametric dimension of the first displacement component
        dim = basisRefs.front().dim();
        // a quadrature rule is defined by the basis for the first displacement component.
        rule = gsQuadrature::get(basisRefs.front(), options);
        // saving necessary info
        T E = options.getReal("YoungsModulus");
        T pr = options.getReal("PoissonsRatio");
        lambda = E * pr / ( ( 1. + pr ) * ( 1. - 2. * pr ) );
        mu     = E / ( 2. * ( 1. + pr ) );
        forceScaling = options.getReal("ForceScaling");
        localStiffening = options.getReal("LocalStiff");
        // linear elasticity tensor
        I = gsMatrix<T>::Identity(dim,dim);
        matrixTraceTensor<T>(C,I,I);
        C *= lambda;
        symmetricIdentityTensor<T>(Ctemp,I);
        C += mu*Ctemp;
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
        // NEED_VALUE to get points in the physical domain for evaluation of the RHS
        // NEED_MEASURE to get the Jacobian determinant values for integration
        // NEED_GRAD_TRANSFORM to get the Jacobian matrix to transform gradient from the parametric to physical domain
        md.flags = NEED_VALUE | NEED_MEASURE | NEED_GRAD_TRANSFORM;
        // Compute image of the quadrature points plus gradient, jacobian and other necessary data
        geo.computeMap(md);
        // find local indices of the displacement basis functions active on the element
        basisRefs.front().active_into(quNodes.col(0),localIndicesDisp);
        N_D = localIndicesDisp.rows();
        // Evaluate displacement basis functions and their derivatives on the element
        basisRefs.front().evalAllDers_into(quNodes,1,basisValuesDisp);
        // Evaluate right-hand side at the image of the quadrature points
        pde_ptr->rhs()->eval_into(md.values[0],forceValues);
    }

    inline void assemble(gsDomainIteratorWrapper<T> & element,
                         const gsVector<T> & quWeights)
    {
        GISMO_UNUSED(element);
        // initialize local matrix and rhs
        localMat.setZero(dim*N_D,dim*N_D);
        localRhs.setZero(dim*N_D,1);
        // Loop over the quadrature nodes
        for (index_t q = 0; q < quWeights.rows(); ++q)
        {
            // Multiply quadrature weight by the geometry measure
            const T weightForce = quWeights[q] * md.measure(q);
            const T weightBody = quWeights[q] * pow(md.measure(q),1-localStiffening);
            // Compute physical gradients of basis functions at q as a dim x numActiveFunction matrix
            transformGradients(md,q,basisValuesDisp[1],physGrad);
            // loop over active basis functions (v_j)
            for (index_t i = 0; i < N_D; i++)
            {
                // stiffness matrix K = B_i^T * C * B_j;
                setB<T>(B_i,I,physGrad.col(i));
                tempK = B_i.transpose() * C; //.reshapeCol(q,dim,dim)
                // loop over active basis functions (v_j)
                for (index_t j = 0; j < N_D; j++)
                {
                    setB<T>(B_j,I,physGrad.col(j));
                    K = tempK * B_j;
                    for (short_t di = 0; di < dim; ++di)
                        for (short_t dj = 0; dj < dim; ++dj)
                            localMat(di*N_D+i,dj*N_D+j) += weightBody * K(di,dj);
                }
            }
            // rhs contribution
            for (short_t d = 0; d < dim; ++d)
                localRhs.middleRows(d*N_D,N_D).noalias() += weightForce * forceScaling * forceValues(d,q) * basisValuesDisp[0].col(q) ;
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
        system.pushToMatrix(localMat,globalIndices,eliminatedDofs,blockNumbers,blockNumbers);
        // push to the elimination system
        if (elimMat != nullptr)
        {
            index_t globalI,globalElimJ;
            index_t elimSize = 0;
            for (short_t dJ = 0; dJ < dim; ++dJ)
            {
                for (short_t dI = 0; dI < dim; ++dI)
                    for (index_t i = 0; i < N_D; ++i)
                        if (system.colMapper(dI).is_free_index(globalIndices[dI].at(i)))
                        {
                            system.mapToGlobalRowIndex(localIndicesDisp.at(i),patchIndex,globalI,dI);
                            for (index_t j = 0; j < N_D; ++j)
                                if (!system.colMapper(dJ).is_free_index(globalIndices[dJ].at(j)))
                                {
                                    globalElimJ = system.colMapper(dJ).global_to_bindex(globalIndices[dJ].at(j));
                                    elimMat->coeffRef(globalI,elimSize+globalElimJ) += localMat(N_D*dI+i,N_D*dJ+j);
                                }
                        }
                elimSize += eliminatedDofs[dJ].rows();
            }
        }
    }

protected:
    // problem info
    short_t dim;
    // number of displacement basis functions active at the current element
    index_t N_D;
    const gsBasePde<T> * pde_ptr;
    // Lame coefficients and force scaling factor
    T lambda, mu, forceScaling, localStiffening;
    // geometry mapping
    gsMapData<T> md;
    // local components of the global linear system
    gsMatrix<T> localMat;
    gsMatrix<T> localRhs;
    // local indices (at the current patch) of the displacement basis functions active at the current element
    gsMatrix<index_t> localIndicesDisp;
    // values and derivatives of displacement basis functions at quadrature points at the current element
    // values are stored as a N_D x numQuadPoints matrix; not sure about derivatives, must be smth like N_D*dim x numQuadPoints
    std::vector<gsMatrix<T> > basisValuesDisp;
    // RHS values at quadrature points at the current element; stored as a dim x numQuadPoints matrix
    gsMatrix<T> forceValues;
    // elimination matrix to efficiently change Dirichlet degrees of freedom
    gsSparseMatrix<T> * elimMat;
    // all temporary matrices defined here for efficiency
    gsMatrix<T> C, Ctemp,physGrad, B_i, tempK, B_j, K, I;
    // containers for global indices
    std::vector< gsMatrix<index_t> > globalIndices;
    gsVector<index_t> blockNumbers;
};

} // namespace gismo
