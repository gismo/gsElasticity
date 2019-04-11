/** @file gsVisitorElThermoBoundary.h

    @brief Visitor class for surface integration of the thermal stress.

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
class gsVisitorElThermoBoundary
{
public:

    gsVisitorElThermoBoundary(const gsPde<T> & pde_,
                              boxSide s,
                              const gsFunctionSet<T> & temperatureField_) :
        side(s),
        temperatureField(temperatureField_)
    {
        pde_ptr = static_cast<const gsPoissonPde<T>*>(&pde_);
    }

    void initialize(const gsBasis<T>& basis,
                    const index_t patchIndex,
                    const gsOptionList & options,
                    gsQuadRule<T> & rule)
    {
        rule = gsQuadrature::get(basis, options, side.direction());
        md.flags = NEED_VALUE | NEED_MEASURE;
        patch = patchIndex;
        dim = basis.dim();

        paramTemp = options.getSwitch("ParamTemp");
        initTemp = options.getReal("InitTemp");
        thermalExpCoef = options.getReal("ThExpCoef");

        T E = options.getReal("YoungsModulus");
        T pr = options.getReal("PoissonsRatio");
        lambda = E * pr / ( ( 1. + pr ) * ( 1. - 2. * pr ) );
        mu     = E / ( 2. * ( 1. + pr ) );
    }

    inline void evaluate(const gsBasis<T> & basis, // to do: more unknowns
                         const gsGeometry<T> & geo,
                         const gsMatrix<T> & quNodes)
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
        // Evaluate temperature
        if (paramTemp) // evaluate temperature in the parametric domain
            temperatureField.piece(patch).eval_into(quNodes,tempValues);
        else           // evaluate temperature in the physical domain
            temperatureField.eval_into(md.values[0],tempValues);
        // Initialize local matrix/rhs
        localRhs.setZero(dim*numActiveFunctions,1);
    }

    inline void assemble(gsDomainIterator<T> & element,
                         const gsVector<T> & quWeights)
    {
        index_t N = numActiveFunctions;

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
                localRhs.middleRows(d*N,N).noalias() += weight * unormal(d,0) * basisValues.col(q) ;
        }
    }

    inline void localToGlobal(const int patchIndex,
                              const std::vector<gsMatrix<T> > & eliminatedDofs,
                              gsSparseSystem<T>     & system)
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
    boxSide side;
    const gsPoissonPde<T> * pde_ptr;
    // geometry mapping
    gsMapData<T> md;

    // Lame coefficients
    T lambda, mu;

    // temperature info
    const gsFunctionSet<T> & temperatureField;
    bool paramTemp;
    T initTemp, thermalExpCoef;

    // local components of the global linear system
    gsMatrix<T> localRhs;

    // Stored evaluations
    index_t numActiveFunctions;
    gsMatrix<unsigned> localIndices;
    gsMatrix<T> basisValues;
    gsMatrix<T> tempValues;

}; //class definition ends

} // namespace ends

