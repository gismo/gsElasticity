/** @file gsVisitorElMass.h

    @brief Visitor class for the mass matrix assembly for elasticity problems.

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

namespace gismo
{

template <class T>
class gsVisitorElMass : public gsVisitorBaseElasticity<T>
{
public:
    typedef gsVisitorBaseElasticity<T> Base;

    gsVisitorElMass(const gsPde<T> & pde_)
        : Base(pde_,true,false) {}

    void initialize(const gsBasisRefs<T> & basisRefs,
                    const index_t patchIndex,
                    const gsOptionList & options,
                    gsQuadRule<T> & rule)
    {
        Base::initialize(basisRefs,patchIndex,options,rule);
        density = options.getReal("Density");
    }

    inline void evaluate(const gsBasisRefs<T> & basisRefs,
                         const gsGeometry<T> & geo,
                         const gsMatrix<T> & quNodes)
    {
        // store quadrature points of the element for geometry evaluation
        md.points = quNodes;
        // NEED_MEASURE to get the Jacobian determinant values for integration
        md.flags = NEED_MEASURE;
        // Compute the geometry mapping at the quadrature points
        geo.computeMap(md);
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
        // initialize local matrix and rhs
        localMat.setZero(dim*N_D,dim*N_D);
        for (index_t q = 0; q < quWeights.rows(); ++q)
        {
            // Multiply quadrature weight by the geometry measure
            const T weight = density * quWeights[q] * md.measure(q);
            for (index_t i = 0; i < N_D; ++i)
                for (index_t j = 0; j < N_D; ++j)
                    for (short_t d = 0; d < dim; ++d)
                        localMat(d*N_D+i,d*N_D+j) += weight*basisValuesDisp[0](i,q)*basisValuesDisp[0](j,q);
        }
    }

protected:
    //------ inherited ------//
    using Base::dim;
    using Base::pde_ptr;
    using Base::md;
    using Base::localMat;
    using Base::localIndices;
    using Base::N_D;
    using Base::basisValuesDisp;

    //------ class specific ----//
    //density
    T density;
};

} // namespace gismo

