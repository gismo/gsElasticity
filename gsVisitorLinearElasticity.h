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

#include <gsElasticity/gsVisitorBaseElasticity.h>

namespace gismo
{

template <class T>
class gsVisitorLinearElasticity : public gsVisitorBaseElasticity<T>
{
public:
    typedef gsVisitorBaseElasticity<T> Base;

    gsVisitorLinearElasticity(const gsPde<T> & pde_, bool assembleMatrix_)
        : Base(pde_,assembleMatrix_) {}

    void initialize(const gsBasisRefs<T> & basisRefs,
                    const index_t patchIndex,
                    const gsOptionList & options,
                    gsQuadRule<T> & rule)
    {
        Base::initialize(basisRefs,patchIndex,options,rule);

        T E = options.getReal("YoungsModulus");
        T pr = options.getReal("PoissonsRatio");
        lambda = E * pr / ( ( 1. + pr ) * ( 1. - 2. * pr ) );
        mu     = E / ( 2. * ( 1. + pr ) );
        forceScaling = options.getReal("ForceScaling");
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
        localIndices.resize(1);
        basisRefs.front().active_into(quNodes.col(0),localIndices[0]);
        N_D = localIndices[0].rows();
        // Evaluate displacement basis functions and their derivatives on the element
        basisRefs.front().evalAllDers_into(quNodes,1,basisValuesDisp);
        // Evaluate right-hand side at the image of the quadrature points
        pde_ptr->rhs()->eval_into(md.values[0],forceValues);
    }

    inline void assemble(gsDomainIterator<T> & element,
                         const gsVector<T> & quWeights)
    {
        // initialize local matrix and rhs
        if (assembleMatrix)
            localMat.setZero(dim*N_D,dim*N_D);
        localRhs.setZero(dim*N_D,1);
        // linear elasticity tensor
        gsMatrix<T> C;
        if (assembleMatrix)
            Base::setC(C,gsMatrix<T>::Identity(dim,dim),lambda,mu);
        // Loop over the quadrature nodes
        for (index_t q = 0; q < quWeights.rows(); ++q)
        {
            // Multiply quadrature weight by the geometry measure
            const T weight = quWeights[q] * md.measure(q);
            if (assembleMatrix)
            {
                // Compute physical gradients of basis functions at q as a dim x numActiveFunction matrix
                gsMatrix<T> physGrad;
                transformGradients(md,q,basisValuesDisp[1],physGrad);
                // loop over active basis functions (v_j)
                for (index_t i = 0; i < N_D; i++)
                {
                    // stiffness matrix K = B_i^T * C * B_j;
                    gsMatrix<T> B_i;
                    Base::setB(B_i,gsMatrix<T>::Identity(dim,dim),physGrad.col(i));
                    gsMatrix<T> tempK = B_i.transpose() * C;
                    // loop over active basis functions (v_j)
                    for (index_t j = 0; j < N_D; j++)
                    {
                        gsMatrix<T> B_j;
                        Base::setB(B_j,gsMatrix<T>::Identity(dim,dim),physGrad.col(j));
                        gsMatrix<T> K = tempK * B_j;

                        for (short_t di = 0; di < dim; ++di)
                            for (short_t dj = 0; dj < dim; ++dj)
                                localMat(di*N_D+i,dj*N_D+j) += weight * K(di,dj);
                    }
                }
            }
            // rhs contribution
            for (short_t d = 0; d < dim; ++d)
                localRhs.middleRows(d*N_D,N_D).noalias() += weight * forceScaling * forceValues(d,q) * basisValuesDisp[0].col(q) ;
        }
    }

protected:
    //------ inherited ------//
    using Base::dim;
    using Base::pde_ptr;
    using Base::assembleMatrix;
    using Base::md;
    using Base::localMat;
    using Base::localRhs;
    using Base::localIndices;
    using Base::N_D;
    using Base::basisValuesDisp;
    using Base::forceValues;

    //------ class specific ----//
    // Lame coefficients, density and force scaling factor
    T lambda, mu, forceScaling;
};

} // namespace gismo
