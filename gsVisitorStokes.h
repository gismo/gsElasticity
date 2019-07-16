/** @file gsVisitorStokes.h

    @brief Visitor class for volumetric integration of the Stokes system.

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s):
        A.Shamanskiy (2016 - ...., TU Kaiserslautern)
*/

#pragma once

#include <gsElasticity/gsVisitorBaseElasticity.h>

namespace gismo
{

template <class T>
class gsVisitorStokes : public gsVisitorBaseElasticity<T>
{
public:
    typedef gsVisitorBaseElasticity<T> Base;

    gsVisitorStokes(const gsPde<T> & pde_)
        : Base(pde_,true) {}

    void initialize(const gsBasisRefs<T> & basisRefs,
                    const index_t patchIndex,
                    const gsOptionList & options,
                    gsQuadRule<T> & rule)
    {
        Base::initialize(basisRefs,patchIndex,options,rule);

        viscosity = options.getReal("Viscosity");
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
        localIndices.resize(2);
        basisRefs.front().active_into(quNodes.col(0),localIndices[0]);
        N_D = localIndices[0].rows();
        basisRefs.back().active_into(quNodes.col(0), localIndices[1]);
        N_P = localIndices[1].rows();
        // Evaluate displacement basis functions and their derivatives on the element
        basisRefs.front().evalAllDers_into(quNodes,1,basisValuesDisp);
        // Evaluate pressure basis functions on the element
        basisRefs.back().eval_into(quNodes,basisValuesPres);
        // Evaluate right-hand side at the image of the quadrature points
        pde_ptr->rhs()->eval_into(md.values[0],forceValues);
    }

    inline void assemble(gsDomainIterator<T>    & element,
                         gsVector<T> const      & quWeights)
    {
        // Initialize local matrix/rhs                          // A | B^T
        localMat.setZero(dim*N_D + N_P, dim*N_D + N_P);         // --|--    matrix structure
        localRhs.setZero(dim*N_D + N_P,1);                      // B | 0

        // Loop over the quadrature nodes
        for (index_t q = 0; q < quWeights.rows(); ++q)
        {
            // Multiply quadrature weight by the geometry measure
            const T weight = quWeights[q] * md.measure(q);

            // Compute physical gradients of basis functions at q as a dim x numActiveFunction matrix
            gsMatrix<T> physGradDisp;
            transformGradients(md, q, basisValuesDisp[1], physGradDisp);

            // matrix A
            for (short_t d = 0; d < dim; ++d)
                localMat.block(d*N_D,d*N_D,N_D,N_D) += (weight*viscosity*(physGradDisp.transpose() * physGradDisp)).block(0,0,N_D,N_D);

            // matrix B
            for (short_t d = 0; d < dim; ++d)
            {
                gsMatrix<> block = weight*basisValuesPres.col(q)*physGradDisp.row(d);
                localMat.block(dim*N_D,d*N_D,N_P,N_D) += block.block(0,0,N_P,N_D);
                localMat.block(d*N_D,dim*N_D,N_D,N_P) += block.transpose().block(0,0,N_D,N_P);
            }

            // rhs contribution
            for (short_t d = 0; d < dim; ++d)
                localRhs.middleRows(d*N_D,N_D).noalias() += weight * forceValues(d,q) * basisValuesDisp[0].col(q) ;
        }
    }

protected:
    //------ inherited ------//
    using Base::dim;
    using Base::pde_ptr;
    using Base::md;
    using Base::localMat;
    using Base::localRhs;
    using Base::localIndices;
    using Base::N_D;
    using Base::basisValuesDisp;
    using Base::forceValues;

    //------ class specific ----//
    T viscosity;

    // number of pressure basis functions active at the current element
    index_t N_P;
    // values of pressure basis functions active at the current element;
    // stores as a N_P x numQuadPoints matrix
    gsMatrix<T> basisValuesPres;
};

} // namespace gismo
