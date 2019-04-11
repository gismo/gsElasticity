/** @file gsVisitorElThermo.h

    @brief Visitor class for volumetric integration of the thermal stress.

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): A. Shamanskiy
*/

#pragma once

#include <gsAssembler/gsQuadrature.h>
#include <gsCore/gsFuncData.h>
#include <gsPde/gsPoissonPde.h>

namespace gismo
{

template <class T>
class gsVisitorElThermo
{
public:

    gsVisitorElThermo(const gsPde<T> & pde_,
                      const gsFunctionSet<T> & temperatureField_) :
        temperatureField(temperatureField_)
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
        patch = patchIndex;
        dim = basis.dim();

        paramTemp = options.getSwitch("ParamTemp");
        thermalExpCoef = options.getReal("ThExpCoef");

        T E = options.getReal("YoungsModulus");
        T pr = options.getReal("PoissonsRatio");
        lambda = E * pr / ( ( 1. + pr ) * ( 1. - 2. * pr ) );
        mu     = E / ( 2. * ( 1. + pr ) );
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
        basis.evalAllDers_into(quNodes, 1, basisValues);
        // Compute image of the quadrature points plus gradient, jacobian and other necessary data
        geo.computeMap(md);
        // Compute temperature gradients
        if (paramTemp) // evaluate gradients in the parametric domain
            temperatureField.piece(patch).deriv_into(quNodes,tempGrads);
        else           // evaluate gradients in the physical domain
            temperatureField.piece(patch).deriv_into(md.values[0],tempGrads);
        // Initialize local matrix/rhs
        localRhs.setZero(dim*numActiveFunctions, 1);
    }

    inline void assemble(gsDomainIterator<T> & element,
                         const gsVector<T> & quWeights)
    {
        // values of basis functions, 1 x numActiveFunctions
        gsMatrix<T> & basisVals = basisValues[0];
        index_t N = numActiveFunctions;

        // Loop over the quadrature nodes
        for (index_t q = 0; q < quWeights.rows(); ++q)
        {
            // Multiply quadrature weight by the geometry measure
            const T weight = thermalExpCoef*(2*mu+dim*lambda)*quWeights[q]*md.measure(q);

            if (paramTemp) // transform temperature gradients to the physical domain
            {
                // temperature gradient at one point in the physical domain, dim x 1
                gsMatrix<T> physGrad;
                transformGradients(md,q,tempGrads,physGrad);
                for (index_t d = 0; d < dim; ++d)
                    localRhs.middleRows(d*N,N).noalias() -= weight * physGrad(d,0) * basisVals.col(q);
            }
            else  // use temperature gradients as they are
                for (index_t d = 0; d < dim; ++d)
                    localRhs.middleRows(d*N,N).noalias() -= weight * tempGrads(d,q) * basisVals.col(q);
        }
    }

    inline void localToGlobal(const int patchIndex,
                              const std::vector<gsMatrix<T> > & eliminatedDofs,
                              gsSparseSystem<T> & system)
    {
        std::vector< gsMatrix<unsigned> > globalIndices(dim,localIndices);
        gsVector<size_t> blockNumbers(dim);
        for (index_t d = 0; d < dim; ++d)
        {
           system.mapColIndices(localIndices, patchIndex, globalIndices[d], d);
           blockNumbers.at(d) = d;
        }
        system.pushToRhs(localRhs,globalIndices,blockNumbers);
    }

protected:
    // general problem info
    index_t dim;
    index_t patch;
    const gsPoissonPde<T> * pde_ptr;
    // geometry mapping
    gsMapData<T> md;

    // Lame coefficients
    T lambda, mu;

    // Temperature info
    const gsFunctionSet<T> & temperatureField;
    bool paramTemp;
    T thermalExpCoef;

    // local components of the global linear system
    gsMatrix<T> localRhs;

    // Stored evaluations
    index_t numActiveFunctions;
    gsMatrix<unsigned> localIndices;
    std::vector<gsMatrix<T> > basisValues;
    gsMatrix<T> tempGrads;

}; //class definition ends

} // namespace ends
