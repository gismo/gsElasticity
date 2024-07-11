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

#include <gsAssembler/gsQuadrature.h>
#include <gsCore/gsFuncData.h>
#include <gsElasticity/gsBasePde.h>

namespace gismo
{

template <class T>
class gsVisitorStokes
{
public:

    gsVisitorStokes(const gsPde<T> & pde_)
    : dim(), pde_ptr(static_cast<const gsBasePde<T>*>(&pde_)),
      viscosity(), density(), N_V(), N_P()
    {}

    void initialize(const gsBasisRefs<T> & basisRefs,
                    const index_t patchIndex,
                    const gsOptionList & options,
                    gsQuadRule<T> & rule)
    {
        // parametric dimension of the first displacement component
        dim = basisRefs.front().dim();
        // a quadrature rule is defined by the basis for the first velocity component.
        // the same rule is used for the presure
        rule = gsQuadrature::get(basisRefs.front(), options);
        // saving necessary info
        viscosity = options.getReal("Viscosity");
        density = options.getReal("Density");
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
        // find local indices of the velocity and pressure basis functions active on the element
        basisRefs.front().active_into(quNodes.col(0),localIndicesVel);
        N_V = localIndicesVel.rows();
        basisRefs.back().active_into(quNodes.col(0), localIndicesPres);
        N_P = localIndicesPres.rows();
        // Evaluate velocity basis functions and their derivatives on the element
        basisRefs.front().evalAllDers_into(quNodes,1,basisValuesVel);
        // Evaluate pressure basis functions on the element
        basisRefs.back().eval_into(quNodes,basisValuesPres);
        // Evaluate right-hand side at the image of the quadrature points
        pde_ptr->rhs()->eval_into(md.values[0],forceValues);
    }

    inline void assemble(gsDomainIterator<T> & element,
                         const gsVector<T> & quWeights)
    {
        // Initialize local matrix/rhs                          // A | B^T
        localMat.setZero(dim*N_V + N_P, dim*N_V + N_P);         // --|--    matrix structure
        localRhs.setZero(dim*N_V + N_P,1);                      // B | 0

        // Loop over the quadrature nodes
        for (index_t q = 0; q < quWeights.rows(); ++q)
        {
            // Multiply quadrature weight by the geometry measure
            const T weight = quWeights[q] * md.measure(q);
            // Compute physical gradients of the velocity basis functions at q as a dim x numActiveFunction matrix
            transformGradients(md, q, basisValuesVel[1], physGradVel);
            // matrix A
            block = weight*viscosity*density * physGradVel.transpose()*physGradVel;
            for (short_t d = 0; d < dim; ++d)
                localMat.block(d*N_V,d*N_V,N_V,N_V) += block.block(0,0,N_V,N_V);
            // matrix B
            for (short_t d = 0; d < dim; ++d)
            {
                block = weight*basisValuesPres.col(q)*physGradVel.row(d);
                localMat.block(dim*N_V,d*N_V,N_P,N_V) -= block.block(0,0,N_P,N_V);
                localMat.block(d*N_V,dim*N_V,N_V,N_P) -= block.transpose().block(0,0,N_V,N_P);
            }
            // rhs contribution
            for (short_t d = 0; d < dim; ++d)
                localRhs.middleRows(d*N_V,N_V).noalias() += weight *density * forceValues(d,q) * basisValuesVel[0].col(q) ;
        }
    }

    inline void localToGlobal(const int patchIndex,
                              const std::vector<gsMatrix<T> > & eliminatedDofs,
                              gsSparseSystem<T> & system)
    {
        // computes global indices for velocity components
        for (short_t d = 0; d < dim; ++d)
        {
            system.mapColIndices(localIndicesVel, patchIndex, globalIndices[d], d);
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
    T viscosity, density;
    // geometry mapping
    gsMapData<T> md;
    // local components of the global linear system
    gsMatrix<T> localMat;
    gsMatrix<T> localRhs;
    // local indices (at the current patch) of basis functions active at the current element
    gsMatrix<index_t> localIndicesVel;
    gsMatrix<index_t> localIndicesPres;
    // number of velocity and pressure basis functions active at the current element
    index_t N_V, N_P;
    // values and derivatives of velocity basis functions at quadrature points at the current element
    // values are stored as a N_V x numQuadPoints matrix; not sure about derivatives, must be smth like N_V*dim x numQuadPoints
    std::vector<gsMatrix<T> > basisValuesVel;
    // values of pressure basis functions active at the current element;
    // stores as a N_P x numQuadPoints matrix
    gsMatrix<T> basisValuesPres;
    // RHS values at quadrature points at the current element; stored as a dim x numQuadPoints matrix
    gsMatrix<T> forceValues;

    // all temporary matrices defined here for efficiency
    gsMatrix<T> block, physGradVel;
    // containers for global indices
    std::vector< gsMatrix<index_t> > globalIndices;
    gsVector<index_t> blockNumbers;
};

} // namespace gismo
