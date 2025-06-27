/** @file gsVisitorThermoBoundary.h

    @brief Visitor class for surface integration of the thermal stress.

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s):
        D. Fusseder  (2012 - 2015, TU Kaiserslautern),
        A.Shamanskiy (2016 - ...., TU Kaiserslautern)
*/

#pragma once

#include <gsAssembler/gsQuadrature.h>
#include <gsCore/gsFuncData.h>

namespace gismo
{

template <class T>
class gsVisitorThermoBoundary
{
public:
    gsVisitorThermoBoundary(boxSide s,
                            const gsFunctionSet<T> & temperatureField_)
        : side(s),
          temperatureField(temperatureField_) {}

    void initialize(const gsBasisRefs<T>& basisRefs,
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
        patch = patchIndex;
        paramTemp = options.getSwitch("ParamTemp");
        initTemp = options.getReal("InitTemp");
        thermalExpCoef = options.getReal("ThExpCoef");
        T E = options.getReal("YoungsModulus");
        T pr = options.getReal("PoissonsRatio");
        lambda = E * pr / ( ( 1. + pr ) * ( 1. - 2. * pr ) );
        mu     = E / ( 2. * ( 1. + pr ) );
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
        // NEED_VALUE to get points in the physical domain for evaluation of the temperature gradients
        // NEED_MEASURE to get the Jacobian determinant values for integration
        if (paramTemp)
            md.flags = NEED_MEASURE;
        else
            md.flags = NEED_VALUE | NEED_MEASURE;
        // Compute image of the quadrature points plus gradient, jacobian and other necessary data
        geo.computeMap(md);
        // Evaluate temperature
        temperatureField.piece(patch).eval_into(quNodes,tempValues);
        // find local indices of the displacement basis functions active on the element
        basisRefs.front().active_into(quNodes.col(0),localIndicesDisp);
        N_D = localIndicesDisp.rows();
        // Evaluate displacement basis functions on the element
        basisRefs.front().eval_into(quNodes,basisValuesDisp);
    }

    inline void assemble(gsDomainIteratorWrapper<T> & element,
                         const gsVector<T> & quWeights)
    {
        GISMO_UNUSED(element);
        // Initialize local matrix/rhs
        localRhs.setZero(dim*N_D,1);
        // loop over the quadrature nodes
        for (index_t q = 0; q < quWeights.rows(); ++q)
        {
            // Compute the outer normal vector on the side
            // normal length equals to the local area measure
            outerNormal(md, q, side, unormal);
            // Collect the factors here: quadrature weight, geometry measure and time factor
            const T weight = thermalExpCoef*(2*mu+dim*lambda)*quWeights[q] * (tempValues.at(q) - initTemp);

            for (index_t d = 0; d < dim; ++d)
                localRhs.middleRows(d*N_D,N_D).noalias() += weight * unormal(d,0) * basisValuesDisp.col(q) ;
        }
    }

    inline void localToGlobal(const int patchIndex,
                              const std::vector<gsMatrix<T> > & eliminatedDofs,
                              gsSparseSystem<T> & system)
    {
        GISMO_UNUSED(eliminatedDofs);
        // computes global indices for displacement components
        for (short_t d = 0; d < dim; ++d)
        {
            system.mapColIndices(localIndicesDisp, patchIndex, globalIndices[d], d);
            blockNumbers.at(d) = d;
        }
        // push to global system
        system.pushToRhs(localRhs,globalIndices,blockNumbers);
    }

protected:
    // problem info
    short_t dim;
    // Lame and thermal expansion coefficients, initial temperature
    T lambda, mu, thermalExpCoef, initTemp;
    // geometry mapping
    gsMapData<T> md;
    // local components of the global linear system
    gsMatrix<T> localRhs;
    // local indices (at the current patch) of the displacement basis functions active at the current element
    gsMatrix<index_t> localIndicesDisp;
    // number of displacement basis functions active at the current element
    index_t N_D;
    // values of displacement basis functions at quadrature points at the current element stored as a N_D x numQuadPoints matrix;
    gsMatrix<T> basisValuesDisp;
    // Temperature info
    index_t patch;
    boxSide side;
    const gsFunctionSet<T> & temperatureField;
    // true if temperature field is defined in the parametric domain; false if in the physical
    bool paramTemp;
    // temperature gradient evaluated at the quadrature points or at their images in the physical domain;
    // stored as a dim x numQuadPoints matrix
    gsMatrix<T> tempGrads;
    // temperature evaluated at the quadrature points or at their images in the physical domain;
    // stored as a 1 x numQuadPoints matrix
    gsMatrix<T> tempValues;

    // all temporary matrices defined here for efficiency
    gsVector<T> unormal;
    // containers for global indices
    std::vector< gsMatrix<index_t> > globalIndices;
    gsVector<index_t> blockNumbers;

}; //class definition ends

} // namespace ends

