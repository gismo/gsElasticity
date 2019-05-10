/** @file gsVisitorBaseElasticity.h

    @brief Abstract parent class for different elasticity visitor classes.

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s):
        A.Shamanskiy (2016 - ...., TU Kaiserslautern)
*/

#pragma once

#include <gsAssembler/gsQuadrature.h>

namespace gismo
{

template <class T>
class gsVisitorBaseElasticity
{
public:

    gsVisitorBaseElasticity(const gsPde<T> & pde_, bool assembleMatrix_)
        : pde_ptr(static_cast<const gsPoissonPde<T>*>(&pde_)),
          assembleMatrix(assembleMatrix_) { }

    void initialize(const gsBasisRefs<T> & basisRefs,
                    const index_t patchIndex,
                    const gsOptionList & options,
                    gsQuadRule<T> & rule)
    {
        // parametric dimension of the first displacement component
        dim = basisRefs.front().dim();
        // a quadrature rule is defined by the basis for the first displacement component.
        // the same rule is used for every unknown.
        // for the displacement formulation this is not a probmlem (all displacement components have the same basis).
        // for a mixed formulation with pressure, the displacement basis is of higher polynomial order or finer than the pressure basis,
        // so the same quadrature rule is at least as accurate for pressure as well
        rule = gsQuadrature::get(basisRefs.front(), options);
    }

    virtual inline void evaluate(const gsBasisRefs<T> & basisRefs,
                                 const gsGeometry<T> & geo,
                                 const gsMatrix<T> & quNodes) = 0;

    virtual inline void assemble(gsDomainIterator<T> & element,
                                 const gsVector<T> & quWeights) = 0;


    inline void localToGlobal(const int patchIndex,
                                      const std::vector<gsMatrix<T> > & eliminatedDofs,
                                      gsSparseSystem<T> & system)
    {
        // dim displacement components plus any other variables
        short_t numUnknows = dim + localIndices.size() - 1;
        std::vector< gsMatrix<unsigned> > globalIndices(numUnknows);
        gsVector<size_t> blockNumbers(numUnknows);
        // computes global indices for displacement components
        for (short_t d = 0; d < dim; ++d)
        {
           system.mapColIndices(localIndices[0], patchIndex, globalIndices[d], d);
           blockNumbers.at(d) = d;
        }
        // computes global indices for other variables
        for (short d = dim; d < numUnknows; ++d)
        {
            system.mapColIndices(localIndices[d-dim+1], patchIndex, globalIndices[d], d);
            blockNumbers.at(d) = d;
        }
        // push to global system
        system.pushToRhs(localRhs,globalIndices,blockNumbers);
        if (assembleMatrix)
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
        // componentwise definition
        for (index_t i = 0; i < dimTensor; ++i)
            for (index_t j = 0; j < dimTensor; ++j)
                C(i,j) = a*R(v(i,0),v(i,1))*R(v(j,0),v(j,1)) +
                         b*(R(v(i,0),v(j,0))*R(v(i,1),v(j,1)) +
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

    // auxiliary matrix B such that E:S = B*Svec in the weak form
    // (see Bernal, Calo, Collier, et. at., "PetIGA", ICCS 2013, p. 1610)
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
    const bool assembleMatrix; // true: assemble matrix and rhs; false: assemble only rhs
    // geometry mapping
    gsMapData<T> md;

    // local components of the global linear system
    gsMatrix<T> localMat;
    gsMatrix<T> localRhs;

    // local indices (at the current patch) of basis functions active at the current element
    // first matrix is always for the first displacement component (which is the same for all components)
    // any subsequent is for any additional variable (for example, pressure in mixed formulation)
    std::vector<gsMatrix<unsigned> > localIndices;
    // number of displacement basis functions active at the current element
    index_t N_D;
    // values and derivatives of displacement basis functions at quadrature points at the current element
    // values are stored as a N_D x numQuadPoints matrix; not sure about derivatives, must be smth like N_D*dim x numQuadPoints
    std::vector<gsMatrix<T> > basisValuesDisp;
    // RHS values at quadrature points at the current element; stored as a dim x numQuadPoints matrix
    gsMatrix<T> forceValues;
};

} // namespace gismo
