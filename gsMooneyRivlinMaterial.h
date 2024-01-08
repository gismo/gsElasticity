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
class gsMooneyRivlinMaterial : public gsMaterialBase<T>
{

public:
    using Base = gsMaterialBase<T>;

    gsMooneyRivlinMaterial( const T c1,
                            const T c2,
                            const gsFunctionSet<T> * patches,
                            const gsFunctionSet<T> * deformed)
    :
    Base(patches,deformed),
    m_homogeneous(true)
    {
        m_Emodulus.push_back(c1);
        m_PoissonRatio.push_back(c2);
    }

    gsMooneyRivlinMaterial( std::vector<T> c1,
                            std::vector<T> c2,
                            const gsFunctionSet<T> * patches,
                            const gsFunctionSet<T> * deformed)
    :
    m_Emodulus(give(c1)),
    m_PoissonRatio(give(c2)),
    Base(patches,deformed),
    m_homogeneous(false)
    {
    }

    // mutable util::gsThreaded<gsMatrix<T> > C, Ctemp;

    /// \a u points, \a k patch index
    void eval_matrix_into(const gsMatrix<T>& u, gsMatrix<T>& result, index_t k = 0) const
    {
        // if (m_homogeneous) k = 0;
        // GISMO_ASSERT(m_Emodulus.size() > static_cast<size_t>(k), "Invalid patch index");

        // gsMatrix<T> Fres, Eres;
        // Base::allStrains_into(u,Fres,Eres,k);

        // T lambda = m_Emodulus.at(k) * m_PoissonRatio.at(k) / ( ( 1. + m_PoissonRatio.at(k) ) * ( 1. - 2. * m_PoissonRatio.at(k) ) );
        // T mu     = m_Emodulus.at(k) / ( 2. * ( 1. + m_PoissonRatio.at(k) ) );

        // const index_t sz = (u.rows()+1)*u.rows()/2;
        // result.resize(sz,u.cols());
        // gsMatrix<T> RCG, RCGinv;
        // T J;
        // for (index_t j=0; j!=u.cols(); j++)
        // {
        //     gsAsMatrix<T, Dynamic, Dynamic> F = Fres.reshapeCol(j,u.cols(),u.cols());
        //     gsAsMatrix<T, Dynamic, Dynamic> E = Eres.reshapeCol(j,u.cols(),u.cols());
        //     J = F.determinant();
        //     RCG = F.transpose() * F;
        //     RCGinv = RCG.cramerInverse();

        //     C.mine() *= lambda*J*J;
        //     symmetricIdentityTensor<T>(Ctemp.mine(),RCGinv);
        //     C.mine() += (mu-lambda*(J*J-1)/2)*Ctemp.mine();
        //     result.reshapeCol(j,sz,sz) = C.mine();
        // }
    }

    /// \a u points, \a k patch index
    void eval_vector_into(const gsMatrix<T>& u, gsMatrix<T>& result, index_t k = 0) const
    {
        // if (m_homogeneous) k = 0;
        // GISMO_ASSERT(m_Emodulus.size() > static_cast<size_t>(k), "Invalid patch index");

        // gsMatrix<T> Fres, Eres;
        // Base::allStrains_into(u,Fres,Eres,k);

        // T lambda = m_Emodulus.at(k) * m_PoissonRatio.at(k) / ( ( 1. + m_PoissonRatio.at(k) ) * ( 1. - 2. * m_PoissonRatio.at(k) ) );
        // T mu     = m_Emodulus.at(k) / ( 2. * ( 1. + m_PoissonRatio.at(k) ) );
        // gsMatrix<T> I = gsMatrix<T>::Identity(u.rows(),u.rows());

        // const index_t sz = (u.rows()+1)*u.rows()/2;
        // result.resize(sz,u.cols());
        // gsMatrix<T> RCG, RCGinv;
        // gsMatrix<T> S;
        // T J;
        // for (index_t i=0; i!=u.cols(); i++)
        // {
        //     gsAsMatrix<T, Dynamic, Dynamic> F = Fres.reshapeCol(i,u.rows(),u.rows());
        //     gsAsMatrix<T, Dynamic, Dynamic> E = Eres.reshapeCol(i,u.rows(),u.rows());
        //     J = F.determinant();
        //     RCG = F.transpose() * F;
        //     RCGinv = RCG.cramerInverse();

        //     S = (m_c1.at(k))*RCGinv + (m_c1.at(k)+3*m_c1.at(k))*I;
                // result.reshapeCol(i,sz,sz) = S;
        // }
    }

protected:
    bool m_homogeneous;
    std::vector<T> m_c1, m_c2;

};

}
