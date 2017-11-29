/** @file gsVisitorElPoisson.h

   @brief Allows to introduce a scaling constant in front of the stiffness matrix.

   This file is part of the G+Smo library.

   This Source Code Form is subject to the terms of the Mozilla Public
   License, v. 2.0. If a copy of the MPL was not distributed with this
   file, You can obtain one at http://mozilla.org/MPL/2.0/.

   Author(s): A. Shamanskiy
*/

#pragma once

#include <gsAssembler/gsVisitorPoisson.h>
#include <gsThermoElasticity/gsElPoissonPde.h>

namespace gismo
{

template <class T>
class gsVisitorElPoisson : public gsVisitorPoisson<T>
{
public:

    gsVisitorElPoisson(const gsPde<T> & pde)
        : gsVisitorPoisson<T>(pde)
    {
        pde_ptr = static_cast<const gsElPoissonPde<T>*>(&pde);
    }

    inline void assemble(gsDomainIterator<T>    & element,
                         gsGeometryEvaluator<T> & geoEval,
                         gsVector<T> const      & quWeights)
    {
        gsMatrix<T> & bVals  = basisData[0];
        gsMatrix<T> & bGrads = basisData[1];

        for (index_t k = 0; k < quWeights.rows(); ++k) // loop over quadrature nodes
        {
            // Multiply weight by the geometry measure
            const T weight = quWeights[k] * geoEval.measure(k);

            // Compute physical gradients at k as a Dim x NumActive matrix
            geoEval.transformGradients(k, bGrads, physGrad);

            localRhs.noalias() += weight * ( bVals.col(k) * rhsVals.col(k).transpose() ) ;
            localMat.noalias() += pde_ptr->k() * weight * (physGrad.transpose() * physGrad);
        }
    }

protected:

    const gsElPoissonPde<T> * pde_ptr;
    using gsVisitorPoisson<T>::basisData;
    using gsVisitorPoisson<T>::physGrad;
    using gsVisitorPoisson<T>::localRhs;
    using gsVisitorPoisson<T>::localMat;
    using gsVisitorPoisson<T>::rhsVals;

};

} // namespace ends;
