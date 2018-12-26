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

#include <gsAssembler/gsQuadrature.h>
#include <gsCore/gsFuncData.h>

namespace gismo
{

template <class T>
class gsVisitorElasticityNeumann
{
public:

    gsVisitorElasticityNeumann(const gsPde<T> & , const boundary_condition<T> & s)
        : neumannFunction_ptr( s.function().get() ),
          patchSide( s.side() )
    {}

    void initialize(const gsBasis<T>   & basis,
                    const index_t ,
                    const gsOptionList & options,
                    gsQuadRule<T>      & rule)
    {
        rule = gsQuadrature::get(basis, options, patchSide.direction());
        md.flags = NEED_VALUE|NEED_MEASURE;

        dim = basis.dim();
        timeFactor = options.getReal("TimeFactor");
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
        basis.eval_into(quNodes, basisValues);
        // Compute image of the quadrature points plus gradient, jacobian and other necessary data
        geo.computeMap(md);
        // Evaluate the Neumann functon on the images of the quadrature points
        neumannFunction_ptr->eval_into(md.values[0], neumannValues);
        // Initialize local matrix/rhs
        localRhs.setZero(dim*numActiveFunctions,1);
    }

    inline void assemble(gsDomainIterator<T>    & element,
                         gsVector<T> const      & quWeights)
    {

        for (index_t q = 0; q < quWeights.rows(); ++q) // loop over quadrature nodes
        {
			// Compute the outer normal vector on the side
            // normal length equals to the local area measure
            gsVector<T> unormal;
            outerNormal(md, q, patchSide, unormal);

            // Coolect of the factors here: quadrature weight, geometry measure and time factor
            const T weight = quWeights[q] * unormal.norm() * timeFactor;

            for (index_t d = 0; d < dim; ++d)
                localRhs.middleRows(d*numActiveFunctions,numActiveFunctions).noalias() +=
                    weight * neumannValues(d,q) * basisValues.col(q) ;
        }
    }

    inline void localToGlobal(const int patchIndex,
                              const std::vector<gsMatrix<T> > & ,
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
    const gsFunction<T> * neumannFunction_ptr;
    boxSide patchSide;
    // geometry mapping
    gsMapData<T> md;

    // Time factor
    T timeFactor;

    // local components of the global linear system
    gsMatrix<T> localMat;
    gsMatrix<T> localRhs;

    // Stored evaluations
    index_t numActiveFunctions;
    gsMatrix<unsigned> localIndices;
    gsMatrix<T> basisValues;
    gsMatrix<T> neumannValues;
};

} // namespace gismo
