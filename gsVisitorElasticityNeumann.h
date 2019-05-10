/** @file gsVisitorElasticityNeumann.h

    @brief Visitor class for the surface load integration.

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
#include <gsCore/gsFuncData.h>

namespace gismo
{

template <class T>
class gsVisitorElasticityNeumann : public gsVisitorBaseElasticity<T>
{
public:
    typedef gsVisitorBaseElasticity<T> Base;

    gsVisitorElasticityNeumann(const gsPde<T> & pde_,
                               const boundary_condition<T> & s)
        : Base(pde_,false), // false = no matrix assembly
          neumannFunction_ptr(s.function().get()),
          patchSide(s.side()) {}

    void initialize(const gsBasisRefs<T> & basisRefs,
                    const index_t patchIndex,
                    const gsOptionList & options,
                    gsQuadRule<T> & rule)
    {
        Base::initialize(basisRefs,patchIndex,options,rule);
        // storing necessary info
        forceScaling = options.getReal("ForceScaling");
    }

    inline void evaluate(const gsBasisRefs<T> & basisRefs,
                         const gsGeometry<T> & geo,
                         const gsMatrix<T> & quNodes)
    {
        // store quadrature points of the element for geometry evaluation
        md.points = quNodes;
        // NEED_VALUE to get points in the physical domain for evaluation of the load
        // NEED_MEASURE to get the Jacobian determinant values for integration
        md.flags = NEED_VALUE | NEED_MEASURE;
        // Compute image of the quadrature points plus gradient, jacobian and other necessary data
        geo.computeMap(md);
        // Evaluate the Neumann functon on the images of the quadrature points
        neumannFunction_ptr->eval_into(md.values[0], neumannValues);
        // find local indices of the displacement basis functions active on the element
        localIndices.resize(1);
        basisRefs.front().active_into(quNodes.col(0),localIndices[0]);
        N_D = localIndices[0].rows();
        // Evaluate basis functions on element
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
            outerNormal(md, q, patchSide, unormal);
            // Collect the factors here: quadrature weight and geometry measure
            const T weight = quWeights[q] * unormal.norm() * forceScaling;

            for (short_t d = 0; d < dim; ++d)
                localRhs.middleRows(d*N_D,N_D).noalias() += weight * neumannValues(d,q) * basisValuesDisp[0].col(q);
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
    T forceScaling;

    // Neumann BC info
    const gsFunction<T> * neumannFunction_ptr;
    boxSide patchSide;
    gsMatrix<T> neumannValues;
};

} // namespace gismo
