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

        T E = options.getReal("YoungsModulus");
        T pr = options.getReal("PoissonsRatio");
        lambda = E * pr / ( ( 1. + pr ) * ( 1. - 2. * pr ) );
        mu     = E / ( 2. * ( 1. + pr ) );
        forceScaling = options.getReal("ForceScaling");
    }

    inline void evaluate(const gsBasis<T> & basis,
                         const gsGeometry<T> & geo,
                         const gsMatrix<T> & quNodes)
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
        // values of basis functions, 1 x numActiveFunctions
        gsMatrix<T> & basisVals  = basisValues[0];
        // gradients of basis functions, dim x numActiveFunctions
        gsMatrix<T> & basisGrads = basisValues[1];

        index_t N = numActiveFunctions;

        // elasticity tensor
        gsMatrix<T> C;
        setC(C,gsMatrix<T>::Identity(dim,dim),lambda,2*mu);
        // Loop over the quadrature nodes
        for (index_t q = 0; q < quWeights.rows(); ++q)
        {
            // Multiply quadrature weight by the geometry measure
            const T weight = quWeights[q] * md.measure(q);
            // Compute physical gradients of basis functions at q as a dim x numActiveFunction matrix
            gsMatrix<T> physGrad;
            transformGradients(md, q, basisGrads, physGrad);
            // loop over active basis functions (v_j)
            for (index_t i = 0; i < N; i++)
            {
                // stiffness matrix K = B_i^T * C * B_j;
                gsMatrix<T> B_i;
                setB(B_i,gsMatrix<T>::Identity(dim,dim),physGrad.col(i));
                gsMatrix<T> tempK = B_i.transpose() * C;
                // loop over active basis functions (v_j)
                for (index_t j = 0; j < N; j++)
                {
                    gsMatrix<T> B_j;
                    setB(B_j,gsMatrix<T>::Identity(dim,dim),physGrad.col(j));
                    gsMatrix<T> K = tempK * B_j;

                    for (short_t di = 0; di < dim; ++di)
                        for (short_t dj = 0; dj < dim; ++dj)
                            localMat(di*N+i, dj*N+j) += weight * K(di,dj);
                }
            }

            for (short_t d = 0; d < dim; ++d)
                localRhs.middleRows(d*N,N).noalias() += weight * forceScaling * forceValues(d,q) * basisVals.col(q) ;
        }
    }

    inline void localToGlobal(const int patchIndex,
                              const std::vector<gsMatrix<T> > & eliminatedDofs,
                              gsSparseSystem<T>     & system)
    {
        std::vector< gsMatrix<unsigned> > globalIndices(dim,localIndices);
        gsVector<size_t> blockNumbers(dim);
        for (short_t d = 0; d < dim; ++d)
        {
           system.mapColIndices(localIndices, patchIndex, globalIndices[d], d);
           blockNumbers.at(d) = d;
        }
        system.pushToRhs(localRhs,globalIndices,blockNumbers);
        system.pushToMatrix(localMat,globalIndices,eliminatedDofs,blockNumbers,blockNumbers);
    }

protected:

    // construct an elasticity tensor C in Voigt notation assuming that
    // C = a R x R + b R (*) R in tensor notation,
    // where R is a symmetric second order tensoR (letter T is reserved).
    // (see Bernal, Calo, Collier, et. at., "PetIGA", ICCS 2013, p. 1608)
    inline void setC(gsMatrix<T> & C, const gsMatrix<T> & R, T a, T b)
    {
        index_t dimTensor = dim*(dim+1)/2;
        C.resize(dimTensor,dimTensor);

        // voigt indices
        gsMatrix<unsigned> v(dimTensor,2);
        if (dim == 2)
            v << 0,0,
                 1,1,
                 0,1;
        if (dim == 3)
            v << 0,0,
                 1,1,
                 2,2,
                 0,1,
                 1,2,
                 0,2;

        for (index_t i = 0; i < dimTensor; ++i)
            for (index_t j = 0; j < dimTensor; ++j)
                C(i,j) = a*R(v(i,0),v(i,1))*R(v(j,0),v(j,1)) +
                         b*0.5*(R(v(i,0),v(j,0))*R(v(i,1),v(j,1)) +
                                R(v(i,0),v(j,1))*R(v(i,1),v(j,0)));
    }

    // transform stress tensor S to a vector in Voigt notation
    inline void voigtStress(gsVector<T> & Svec, const gsMatrix<T> & S)
    {
        index_t dimTensor = dim*(dim+1)/2;
        Svec.resize(dimTensor);
        if (dim == 2)
        {
            Svec(0) = S(0,0);
            Svec(1) = S(1,1);
            Svec(2) = S(0,1);
        }
        if (dim == 3)
        {
            Svec(0) = S(0,0);
            Svec(1) = S(1,1);
            Svec(2) = S(2,2);
            Svec(3) = S(0,1);
            Svec(4) = S(1,2);
            Svec(5) = S(0,2);
        }
    }

    // auxiliary matrix B (see Bernal, Calo, Collier, et. at., "PetIGA", ICCS 2013, p. 1610)
    inline void setB( gsMatrix<T> & B, const gsMatrix<T> & F, const gsVector<T> & bGrad)
    {
        index_t dimTensor = dim*(dim+1)/2;
        B.resize(dimTensor,dim);

        for (short_t j = 0; j < dim; ++j)
        {
            for (short_t i = 0; i < dim; ++i)
                B(i,j) = F(j,i) * bGrad(i);

            if (dim == 2)
                B(2,j) = F(j,0) * bGrad(1) + F(j,1) * bGrad(0);

            if (dim == 3)
                for (short_t i = 0; i < dim; ++i)
                {
                    short_t k = (i+1)%dim;
                    B(i+dim,j) = F(j,i) * bGrad(k) + F(j,k) * bGrad(i);
                }
        }
    }

protected:
    // general problem info
    short_t dim;
    const gsPoissonPde<T> * pde_ptr;
    // geometry mapping
    gsMapData<T> md;

    // Lame coefficients, density and force scaling factor
    T lambda, mu, forceScaling;

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
