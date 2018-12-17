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

#include <gsAssembler/gsQuadrature.h>

namespace gismo
{

template <class T>
class gsVisitorLinearElasticity
{
public:

    gsVisitorLinearElasticity(const gsPde<T> & pde_)
    {
        pde_ptr = static_cast<const gsPoissonPde<T>*>(&pde_);
    }

    void initialize(const gsBasis<T> & basis,
                    const index_t patchIndex,
                    const gsOptionList & options,
                    gsQuadRule<T> & rule)
    {
        rule = gsQuadrature::get(basis, options); // harmless slicing occurs here (?)
        md.flags = NEED_VALUE | NEED_MEASURE | NEED_GRAD_TRANSFORM;

        dim = basis.dim();
        rho = options.getReal("Density");
        timefactor = options.getReal("TimeFactor");

        T E = options.getReal("YoungsModulus");
        T pr = options.getReal("PoissonsRatio");
        lambda = E * pr / ( ( 1. + pr ) * ( 1. - 2. * pr ) );
        mu     = E / ( 2. * ( 1. + pr ) );
    }

    inline void evaluate(gsBasis<T> const       & basis, // to do: more unknowns
                         gsGeometry<T> const & geo,
                         gsMatrix<T> const      & quNodes)
    {
        // store quadrature points of the element for geometry evaluation
        md.points = quNodes;
        // find local indices which are active on the element
        basis.active_into(quNodes.col(0), localIndices);
        numActiveFunctions = localIndices.rows();
        // Evaluate basis functions on element
        basis.evalAllDers_into( quNodes, 1, basisValues);
        // Compute image of the quadrature points plus gradient, jacobian and other necessary data
        geo.computeMap(md);
        // Evaluate right-hand side at the image of the quadrature points
        pde_ptr->rhs()->eval_into( md.values[0], forceValues );
        // Initialize local matrix/rhs
        localMat.setZero(dim*numActiveFunctions, dim*numActiveFunctions);
        localRhs.setZero(dim*numActiveFunctions, 1);
    }

    inline void assemble(gsDomainIterator<T>    & element,
                         gsVector<T> const      & quWeights)
    {
        gsMatrix<T> & bVals  = basisValues[0];
        gsMatrix<T> & bGrads = basisValues[1];

        T lam2mu =  lambda + 2 * mu;
        index_t N = numActiveFunctions;

        // Loop over the quadrature nodes
        for (index_t q = 0; q < quWeights.rows(); ++q)
        {
            // Multiply quadrature weight by the geometry measure
            const T weight = quWeights[q] * md.measure(q);
            // Compute physical gradients at q as a dim x numActiveFunction matrix
            gsMatrix<T> physGrad;
            transformGradients(md, q, bGrads, physGrad);

            if (dim == 2)
                for (index_t i = 0; i < N; i++)
                    for (index_t j = 0; j < N; j++)
					{
                        localMat(0*N+i, 0*N+j) += weight *
                            ( lam2mu*physGrad(0,i)*physGrad(0,j) + mu*physGrad(1,i)*physGrad(1,j) );
                        localMat(0*N+i, 1*N+j) += weight *
                            ( lambda*physGrad(0,i)*physGrad(1,j) + mu*physGrad(1,i)*physGrad(0,j) );
                        localMat(1*N+i, 0*N+j) += weight *
                            ( lambda*physGrad(1,i)*physGrad(0,j) + mu*physGrad(0,i)*physGrad(1,j) );
                        localMat(1*N+i, 1*N+j) += weight *
                            ( lam2mu*physGrad(1,i)*physGrad(1,j) + mu*physGrad(0,i)*physGrad(0,j) );
                    }

            else if (dim == 3)
                for (index_t i = 0; i < N; i++)
                    for (index_t j = 0; j < N; j++)
					{
                        localMat(0*N+i, 0*N+j) += weight *
                            ( lam2mu*physGrad(0,i)*physGrad(0,j) + mu*(physGrad(1,i)*physGrad(1,j)+physGrad(2,i)*physGrad(2,j)) );
                        localMat(0*N+i, 1*N+j) += weight *
                            ( lambda*physGrad(0,i)*physGrad(1,j) + mu*physGrad(1,i)*physGrad(0,j) );
                        localMat(0*N+i, 2*N+j) += weight *
                            ( lambda*physGrad(0,i)*physGrad(2,j) + mu*physGrad(2,i)*physGrad(0,j) );
                        localMat(1*N+i, 0*N+j) += weight *
                            ( lambda*physGrad(1,i)*physGrad(0,j) + mu*physGrad(0,i)*physGrad(1,j) );
                        localMat(1*N+i, 1*N+j) += weight *
                            ( lam2mu*physGrad(1,i)*physGrad(1,j) + mu*(physGrad(0,i)*physGrad(0,j)+physGrad(2,i)*physGrad(2,j)) );
                        localMat(1*N+i, 2*N+j) += weight *
                            ( lambda*physGrad(1,i)*physGrad(2,j) + mu*physGrad(2,i)*physGrad(1,j) );
                        localMat(2*N+i, 0*N+j) += weight *
                            ( lambda*physGrad(2,i)*physGrad(0,j) + mu*physGrad(0,i)*physGrad(2,j) );
                        localMat(2*N+i, 1*N+j) += weight *
                            ( lambda*physGrad(2,i)*physGrad(1,j) + mu*physGrad(1,i)*physGrad(2,j) );
                        localMat(2*N+i, 2*N+j) += weight *
                            ( lam2mu*physGrad(2,i)*physGrad(2,j) + mu*(physGrad(1,i)*physGrad(1,j)+physGrad(0,i)*physGrad(0,j)) );
                    }

            for (index_t d = 0; d < dim; ++d)
                localRhs.middleRows(d*N,N).noalias() += weight * rho * forceValues(d,q) * timefactor * bVals.col(q) ;
        }
    }

    inline void localToGlobal(const int patchIndex,
                              const std::vector<gsMatrix<T> > & eliminatedDofs,
                              gsSparseSystem<T>     & system)
    {
        std::vector< gsMatrix<unsigned> > globalIndices(dim,localIndices);
        for (index_t d = 0; d < dim; ++d)
            system.mapColIndices(localIndices, patchIndex, globalIndices[d], d);

        system.push(localMat, localRhs, globalIndices,eliminatedDofs);
    }

protected:
    // general problem info
    index_t dim;
    const gsPoissonPde<T> * pde_ptr;
    // geometry mapping
    gsMapData<T> md;

    // Lame coefficients, density and time factor
    T lambda, mu, rho, timefactor;

    // local components of the global linear system
    gsMatrix<T> localMat;
    gsMatrix<T> localRhs;

    // Stored evaluations
    index_t numActiveFunctions;
    gsMatrix<unsigned> localIndices;
    std::vector<gsMatrix<T> > basisValues;
    gsMatrix<T> forceValues;
};

} // namespace gismo
