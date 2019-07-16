/** @file gsVisitorElThermo.h

    @brief Visitor class for volumetric integration of the thermal stress.

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s):
        D. Fusseder  (2012 - 2015, TU Kaiserslautern),
        A.Shamanskiy (2016 - ...., TU Kaiserslautern)
*/

#pragma once

#include <gsElasticity/gsVisitorBaseElasticity.h>

namespace gismo
{

template <class T>
class gsVisitorElThermo : public gsVisitorBaseElasticity<T>
{
public:
    typedef gsVisitorBaseElasticity<T> Base;

    gsVisitorElThermo(const gsPde<T> & pde_,
                      const gsFunctionSet<T> & temperatureField_)
        : Base(pde_, false), // false = no matrix assembly
          temperatureField(temperatureField_) {}

    void initialize(const gsBasisRefs<T> & basisRefs,
                    const index_t patchIndex,
                    const gsOptionList & options,
                    gsQuadRule<T> & rule)
    {
        Base::initialize(basisRefs,patchIndex,options,rule);
        // storing necessary info
        patch = patchIndex;
        paramTemp = options.getSwitch("ParamTemp");
        thermalExpCoef = options.getReal("ThExpCoef");
        T E = options.getReal("YoungsModulus");
        T pr = options.getReal("PoissonsRatio");
        lambda = E * pr / ( ( 1. + pr ) * ( 1. - 2. * pr ) );
        mu     = E / ( 2. * ( 1. + pr ) );
    }

    inline void evaluate(const gsBasisRefs<T> & basisRefs,
                         const gsGeometry<T> & geo,
                         const gsMatrix<T> & quNodes)
    {        
        // store quadrature points of the element for geometry evaluation
        md.points = quNodes;
        // NEED_VALUE to get points in the physical domain for evaluation of the temperature gradients
        // NEED_MEASURE to get the Jacobian determinant values for integration
        // NEED_GRAD_TRANSFORM to get the Jacobian matrix to transform temperature gradient from the parametric to physical domain
        if (paramTemp)
            md.flags = NEED_MEASURE | NEED_GRAD_TRANSFORM;
        else
            md.flags = NEED_VALUE | NEED_MEASURE;
        // Compute image of the quadrature points plus gradient, jacobian and other necessary data
        geo.computeMap(md);
        // Compute temperature gradients
        if (paramTemp) // evaluate gradients in the parametric domain
            temperatureField.piece(patch).deriv_into(quNodes,tempGrads);
        else           // evaluate gradients in the physical domain
            temperatureField.piece(patch).deriv_into(md.values[0],tempGrads);
        // find local indices of the displacement basis functions active on the element
        localIndices.resize(1);
        basisRefs.front().active_into(quNodes.col(0),localIndices[0]);
        N_D = localIndices[0].rows();
        // Evaluate displacement basis functions on the element
        basisRefs.front().evalAllDers_into(quNodes,0,basisValuesDisp);
    }

    inline void assemble(gsDomainIterator<T> & element,
                         const gsVector<T> & quWeights)
    {
        // Initialize local matrix/rhs
        localRhs.setZero(dim*N_D, 1);
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
                    localRhs.middleRows(d*N_D,N_D).noalias() -= weight * physGrad(d,0) * basisValuesDisp[0].col(q);
            }
            else  // use temperature gradients as they are
                for (index_t d = 0; d < dim; ++d)
                    localRhs.middleRows(d*N_D,N_D).noalias() -= weight * tempGrads(d,q) * basisValuesDisp[0].col(q);
        }
    }

protected:
    //------ inherited ------//
    using Base::dim;
    using Base::md;
    using Base::localRhs;
    using Base::localIndices;
    using Base::N_D;
    using Base::basisValuesDisp;

    //------ class specific ----//
    // Lame coefficients
    T lambda, mu;
    // Temperature info
    index_t patch;
    const gsFunctionSet<T> & temperatureField;
    // true if temperature field is defined in the parametric domain; false if in the physical
    bool paramTemp;
    T thermalExpCoef;
    // temperature gradient evaluated at the quadrature points or at their images in the physical domain;
    // stored as a dim x numQuadPoints matrix
    gsMatrix<T> tempGrads;

}; //class definition ends

} // namespace ends
