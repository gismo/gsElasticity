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

#define WITHEIGEN
#include <muesli/muesli.h>
#include <muesli/tensor.h>
#include <muesli/finitestrain.h>
#undef matrix

namespace gismo
{

template <class T>
class gsMuesliMaterial : public gsMaterialBase<T>
{

public:
    using Base = gsMaterialBase<T>;

    gsMuesliMaterial(   const muesli::finiteStrainMaterial * material,
                        const gsFunctionSet<T> * patches,
                        const gsFunctionSet<T> * defpatches)
    :
    gsMuesliMaterial(material,patches)
    {
        this->updateCurrentState(defpatches);
        this->commitCurrentState();
    }

    gsMuesliMaterial(   const muesli::finiteStrainMaterial * material,
                        const gsFunctionSet<T> * patches)
    :
    Base(patches),
    m_material(material),
    m_convergedTime(0),
    m_currentTime(0)
    {
    }

    // mutable util::gsThreaded<gsMatrix<T> > C, Ctemp;

    /// \a u points, \a k patch index
    void eval_matrix_into(const gsMatrix<T>& u, gsMatrix<T>& result, index_t k = 0) const
    {
        muesli::finiteStrainMP* p = m_material->createMaterialPoint();
        gsMatrix<T> Fresn, Fn, Fresc, Fc;
        gsVector<T> Stmp;
        Base::gradient_into(m_converged,u,Fresn);
        Base::gradient_into(m_deformed,u,Fresc);
        itensor4 tangent;

        const index_t sz = (u.rows()+1)*u.rows()/2;
        GISMO_ASSERT(sz==6,"Now works for 6 only (due to conversion of Ctemp to C).");
        result.resize(sz*sz,u.cols());
        for (index_t k=0; k!= u.cols(); k++)
        {
            // ẗodo: hardcoded sizes
            // Set the converged (previous) deformation gradient
            Fn = Fresn.reshapeCol(k,3,3);
            p->updateCurrentState(m_convergedTime,Fn);
            p->commitCurrentState();

            // Set the current deformation gradient
            Fc = Fresc.reshapeCol(k,3,3);
            p->updateCurrentState(m_currentTime,Fc);

            // compute tangent
            p->convectedTangent(tangent);
            T Ctemp[6][6];
            muesli::tensorToMatrix(tangent,Ctemp);
            gsAsMatrix<T, Dynamic, Dynamic> C = result.reshapeCol(k,sz,sz);
            for (index_t i=0; i!=sz; i++)
                for (index_t j=0; j!=sz; j++)
                    C(i,j) = Ctemp[i][j];
            // result.col(k) = Stmp;
        }

        // if (m_homogeneous) k = 0;
        // GISMO_ASSERT(m_Emodulus.size() > static_cast<size_t>(k), "Invalid patch index");

        // gsMatrix<T> Fresn, Eres;
        // Base::allStrains_into(u,Fresn,Eres,k);

        // T lambda = m_Emodulus.at(k) * m_PoissonRatio.at(k) / ( ( 1. + m_PoissonRatio.at(k) ) * ( 1. - 2. * m_PoissonRatio.at(k) ) );
        // T mu     = m_Emodulus.at(k) / ( 2. * ( 1. + m_PoissonRatio.at(k) ) );

        // const index_t sz = (u.rows()+1)*u.rows()/2;
        // result.resize(sz,u.cols());
        // gsMatrix<T> RCG, RCGinv;
        // T J;
        // for (index_t j=0; j!=u.cols(); j++)
        // {
        //     gsAsMatrix<T, Dynamic, Dynamic> F = Fresn.reshapeCol(j,u.rows(),u.rows());
        //     gsAsMatrix<T, Dynamic, Dynamic> E = Eres.reshapeCol(j,u.rows(),u.rows());
        //     J = F.determinant();
        //     RCG = F.transpose() * F;
        //     RCGinv = RCG.cramerInverse();

        //     C.mine() *= lambda*J*J;
        //     symmetricIdentityTensor<T>(Ctemp.mine(),RCGinv);
        //     C.mine() += (mu-lambda*(J*J-1)/2)*Ctemp.mine();
        //     result.reshapeCol(j,sz,sz) = C.mine();
        // }
        // 
        delete p;
    }

    /// \a u points, \a k patch index
    void eval_vector_into(const gsMatrix<T>& u, gsMatrix<T>& result, index_t k = 0) const
    {
        muesli::finiteStrainMP* p = m_material->createMaterialPoint();
        gsMatrix<T> Fresn, Fn, Fresc, Fc;
        Base::gradient_into(m_converged,u,Fresn);
        Base::gradient_into(m_deformed,u,Fresc);
        istensor stress;

        const index_t sz = u.rows();
        result.resize(sz*sz,u.cols());
        for (index_t k=0; k!= u.cols(); k++)
        {
            // ẗodo: hardcoded sizes
            // Set the converged (previous) deformation gradient
            Fn = Fresn.reshapeCol(k,3,3);
            p->updateCurrentState(m_convergedTime,Fn);
            p->commitCurrentState();

            // Set the current deformation gradient
            Fc = Fresc.reshapeCol(k,3,3);
            p->updateCurrentState(m_currentTime,Fc);

            // compute stress
            p->secondPiolaKirchhoffStress(stress);
            result.reshapeCol(k,sz,sz) = stress;
        }

        // if (m_homogeneous) k = 0;
        // GISMO_ASSERT(m_Emodulus.size() > static_cast<size_t>(k), "Invalid patch index");

        // gsMatrix<T> Fresn, Eres;
        // Base::allStrains_into(u,Fresn,Eres,k);

        // T lambda = m_Emodulus.at(k) * m_PoissonRatio.at(k) / ( ( 1. + m_PoissonRatio.at(k) ) * ( 1. - 2. * m_PoissonRatio.at(k) ) );
        // T mu     = m_Emodulus.at(k) / ( 2. * ( 1. + m_PoissonRatio.at(k) ) );
        // gsMatrix<T> I = gsMatrix<T>::Identity(u.rows(),u.rows());

        // const index_t sz = (u.rows()+1)*u.rows()/2;
        // result.resize(sz,u.cols());
        // gsMatrix<T> RCG, RCGinv;
        // gsVector<T> Stmp;
        // gsMatrix<T> S;
        // T J;
        // for (index_t i=0; i!=u.cols(); i++)
        // {
        //     gsAsMatrix<T, Dynamic, Dynamic> F = Fresn.reshapeCol(i,u.rows(),u.rows());
        //     gsAsMatrix<T, Dynamic, Dynamic> E = Eres.reshapeCol(i,u.rows(),u.rows());
        //     J = F.determinant();
        //     RCG = F.transpose() * F;
        //     RCGinv = RCG.cramerInverse();

        //     S = (m_c1.at(k))*RCGinv + (m_c1.at(k)+3*m_c1.at(k))*I;
        //     voigtStress<T>(Stmp,S);
        //     result.col(i) = Stmp;
        // }
    }

public:

    void updateCurrentState(const gsFunctionSet<T> * deformed, T time = 0.0)
    {
        m_deformed = deformed;
        m_currentTime = time;
    }

    void commitCurrentState()
    {
        m_converged = m_deformed;
        m_convergedTime = m_currentTime;
    }

// protected:
//     void

protected:
    const muesli::finiteStrainMaterial * m_material;
    // m_deformed is the current deformation
    using Base::m_deformed;
    const gsFunctionSet<T> * m_converged;

    // gsFunctionSet<T> m_current;


    T m_currentTime, m_convergedTime;

};

}
