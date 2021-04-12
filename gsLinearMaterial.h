/** @file gsLinearMaterial.h

    @brief

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s):
    O. Weeger    (2012 - 2015, TU Kaiserslautern),
    A.Shamanskiy (2016 - 2020, TU Kaiserslautern),
    H.M.Verhelst (2019 - ...., TU Delft)
*/

#pragma once

#include <gsElasticity/gsMaterialBase.h>
#include <gsElasticity/gsVisitorElUtils.h>
#include <gsUtils/gsThreaded.h>

namespace gismo
{

template <class T>
class gsLinearMaterial : public gsMaterialBase<T>
{

public:
    using Base = gsMaterialBase<T>;

    gsLinearMaterial(   const T E,
                        const T nu)
    : m_homogeneous(true)
    {
        m_Emodulus.push_back(E);
        m_PoissonRatio.push_back(nu);
    }

    gsLinearMaterial(std::vector<T> E,
                     std::vector<T> nu)
    : m_Emodulus(give(E)), m_PoissonRatio(give(nu)), m_homogeneous(false)
    {
    }

    mutable util::gsThreaded<gsMatrix<T> > C, Ctemp;

    /// \a u points, \a k patch index
    void eval_matrix_into(const gsMatrix<T>& u, gsMatrix<T>& result, index_t k = 0) const
    {
        if (m_homogeneous) k = 0;
        GISMO_ASSERT(m_Emodulus.size() > static_cast<size_t>(k), "Invalid patch index");

        const index_t sz = (u.rows()+1)*u.rows()/2;
        result.resize(sz*sz,u.cols());

        T lambda = m_Emodulus.at(k) * m_PoissonRatio.at(k) / ( ( 1. + m_PoissonRatio.at(k) ) * ( 1. - 2. * m_PoissonRatio.at(k) ) );
        T mu     = m_Emodulus.at(k) / ( 2. * ( 1. + m_PoissonRatio.at(k) ) );
        // linear elasticity tensor
        gsMatrix<T> I = gsMatrix<T>::Identity(u.rows(),u.rows());
        matrixTraceTensor<T>(C,I,I);
        C.mine() *= lambda;
        symmetricIdentityTensor<T>(Ctemp,I);
        C.mine() += mu*Ctemp.mine();

        for (index_t j=0; j!=u.cols(); j++)
            result.reshapeCol(j,sz,sz) = C.mine();
    }

        /// \a u points, \a k patch index
        void eval_vector_into(const gsMatrix<T>& u, gsMatrix<T>& result, index_t k = 0) const
    {
        GISMO_NO_IMPLEMENTATION;
    }

protected:
    bool m_homogeneous;
    std::vector<T> m_Emodulus, m_PoissonRatio;

};

}
