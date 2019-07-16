/** @file gsVisitorElThermoBoundary.h

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

#include <gsElasticity/gsVisitorBaseElasticity.h>

#include <gsCore/gsFuncData.h>

namespace gismo
{

template <class T>
class gsVisitorElThermoBoundary : public gsVisitorBaseElasticity<T>
{
public:
    typedef gsVisitorBaseElasticity<T> Base;

    gsVisitorElThermoBoundary(const gsPde<T> & pde_, boxSide s,
                              const gsFunctionSet<T> & temperatureField_)
        : Base (pde_,false), // false = no matrix assembly
          side(s),
          temperatureField(temperatureField_) {}

    void initialize(const gsBasisRefs<T>& basisRefs,
                    const index_t patchIndex,
                    const gsOptionList & options,
                    gsQuadRule<T> & rule)
    {
        Base::initialize(basisRefs,patchIndex,options,rule);
        // storing necessary info
        patch = patchIndex;
        paramTemp = options.getSwitch("ParamTemp");
        initTemp = options.getReal("InitTemp");
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
        if (paramTemp)
            md.flags = NEED_MEASURE;
        else
            md.flags = NEED_VALUE | NEED_MEASURE;
        // Compute image of the quadrature points plus gradient, jacobian and other necessary data
        geo.computeMap(md);
        // Evaluate temperature
        if (paramTemp) // evaluate temperature in the parametric domain
            temperatureField.piece(patch).eval_into(quNodes,tempValues);
        else           // evaluate temperature in the physical domain
            temperatureField.eval_into(md.values[0],tempValues);
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
        localRhs.setZero(dim*N_D,1);
        // loop over the quadrature nodes
        for (index_t q = 0; q < quWeights.rows(); ++q)
        {
            // Compute the outer normal vector on the side
            // normal length equals to the local area measure
            gsVector<T> unormal;
            outerNormal(md, q, side, unormal);
            // Collect the factors here: quadrature weight, geometry measure and time factor
            const T weight = thermalExpCoef*(2*mu+dim*lambda)*quWeights[q] * (tempValues.at(q) - initTemp);

            for (index_t d = 0; d < dim; ++d)
                localRhs.middleRows(d*N_D,N_D).noalias() += weight * unormal(d,0) * basisValuesDisp[0].col(q) ;
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
    boxSide side;
    const gsFunctionSet<T> & temperatureField;
    // true if temperature field is defined in the parametric domain; false if in the physical
    bool paramTemp;
    T thermalExpCoef;
    T initTemp; // fixed initial temperature
    // temperature evaluated at the quadrature points or at their images in the physical domain;
    // stored as a 1 x numQuadPoints matrix
    gsMatrix<T> tempValues;

}; //class definition ends

} // namespace ends

