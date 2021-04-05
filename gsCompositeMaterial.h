/** @file gsCompositeMaterial.h

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
class gsCompositeMaterial : public gsMaterialBase<T>
{

public:
    using Base = gsMaterialBase<T>;

    gsCompositeMaterial(    const gsMatrix<T> G,
                            const T alpha,
                            const T beta,
                            const T gamma)
    : m_homogeneous(true)
    {
        gsMatrix<T> Phi;
        _makePhi(alpha,beta,gamma,Phi);
        m_G.push_back( Phi.transpose()*G*Phi );
    }

    gsCompositeMaterial(std::vector<gsMatrix<T>> G,
                        std::vector<T> alpha,
                        std::vector<T> beta,
                        std::vector<T> gamma)
    : m_homogeneous(false)
    {
        gsMatrix<T> Phi;
        for (index_t k=0; k!=G.size(); k++)
        {
            _makePhi(alpha.at(k),beta.at(k),gamma.at(k),Phi);
            m_G.push_back( Phi.transpose()*G.at(k)*Phi );
        }

    }

    gsCompositeMaterial(    const T Exx,
                            const T Eyy,
                            const T Ezz,
                            const T Gxy,
                            const T Gyz,
                            const T Gxz,
                            const T nuxy,
                            const T nuyz,
                            const T nuxz,
                            const T alpha = 0,
                            const T beta = 0,
                            const T gamma = 0)
    : m_homogeneous(true)
    {
        gsMatrix<T> G,Phi;
        _makeG(Exx,Eyy,Ezz,Gxy,Gyz,Gxz,nuxy,nuyz,nuxz,G);
        _makePhi(alpha,beta,gamma,Phi);
        m_G.push_back( Phi.transpose()*G*Phi );
    }

    gsCompositeMaterial(    std::vector<T> Exx,
                            std::vector<T> Eyy,
                            std::vector<T> Ezz,
                            std::vector<T> Gxy,
                            std::vector<T> Gyz,
                            std::vector<T> Gxz,
                            std::vector<T> nuxy,
                            std::vector<T> nuyz,
                            std::vector<T> nuxz,
                            std::vector<T> alpha,
                            std::vector<T> beta,
                            std::vector<T> gamma)
    : m_homogeneous(false)
    {
        gsMatrix<T> G,Phi;
        for (index_t k=0; k!=Exx.size(); k++)
        {
            _makeG(Exx.at(k),Eyy.at(k),Ezz.at(k),Gxy.at(k),Gyz.at(k),Gxz.at(k),nuxy.at(k),nuyz.at(k),nuxz.at(k),G);
            _makePhi(alpha.at(k),beta.at(k),gamma.at(k),Phi);
            m_G.push_back( Phi.transpose()*G*Phi );
        }
    }

    mutable util::gsThreaded<gsMatrix<T> > C, Ctemp;

    /// \a u points, \a k patch index
    void eval_matrix_into(const gsMatrix<T>& u, gsMatrix<T>& result, index_t k = 0) const
    {
        if (m_homogeneous) k = 0;
        GISMO_ASSERT(m_G.size()<static_cast<size_t>(k), "Invalid patch index");

        const index_t sz = (u.rows()+1)*u.rows()/2;
        result.resize(sz*sz,u.cols());

        // composite linear elasticity tensor
        for (index_t j=0; j!=u.cols(); j++)
            result.reshapeCol(j,sz,sz) = m_G.at(k);
    }

        /// \a u points, \a k patch index
        void eval_vector_into(const gsMatrix<T>& u, gsMatrix<T>& result, index_t k = 0) const
    {
        GISMO_NO_IMPLEMENTATION;
    }

protected:
    void _makeG(const T Exx, const T Eyy, const T Ezz,
                const T Gxy, const T Gxz, const T Gyz,
                const T nuxy,const T nuxz,const T nuyz,
                gsMatrix<T> & G)
    {
        // Note: The stress and strain tensors are ordered as:
        // 2D: [Sxx,Syy,Sxy]
        // 3D: [Sxx,Syy,Szz,Sxy,Syz,Sxz]
        /// Via https://www.efunda.com/formulae/solid_mechanics/mat_mechanics/hooke_orthotropic.cfm

        G.resize(6,6);
        G.setZero();
        T D = (1-nuxy*nuxy-nuyz*nuyz-nuxz*nuxz-2*nuxy*nuyz*nuxz) / (Exx*Eyy*Ezz);
        G(0,0) = (1-nuyz*nuyz) / (Eyy*Ezz); // Gzz
        G(1,1) = (1-nuxz*nuxz) / (Exx*Ezz); // Gyy
        G(2,2) = (1-nuxy*nuxy) / (Exx*Eyy); // Gzz

        // Gij = nuij + nuik*nukj / Ejj*Ekk
        G(0,1) = G(1,0) = (nuxy+nuxz*nuyz) / (Eyy*Ezz); // Gxy
        G(0,2) = G(2,0) = (nuxz+nuxy*nuyz) / (Eyy*Ezz); // Gxz
        G(1,2) = G(2,1) = (nuyz+nuxy*nuxz) / (Exx*Ezz); // Gyz

        G *= 1.0/D;

        G(3,3) = 2.*Gxy;
        G(4,4) = 2.*Gyz;
        G(5,5) = 2.*Gxz;
    }


    /**
     * @brief      Makes the rotation matrix (Tiat-Bryan angle)
     *
     * @param[in]  alpha  The rotation around the x axis
     * @param[in]  beta   The rotation around the y axis
     * @param[in]  gamma  The rotation around the z axis
     * @param[in]  Phi    The rotation matrix
     *
     */
    void _makePhi(  const T alpha,
                    const T beta,
                    const T gamma,
                    gsMatrix<T> & Phi)
    {
        gsMatrix<T,3,3> Ax,Ay,Az,A;
        Ax <<                   1,                  0,                  0,
                                0,   math::cos(alpha),  -math::sin(alpha),
                                0,   math::sin(alpha),   math::cos(alpha);

        Ay <<     math::cos(beta),                  0,    math::sin(beta),
                                0,                  1,                  0,
                 -math::sin(beta),                  0,    math::cos(beta);

        Az <<    math::cos(gamma),  -math::sin(gamma),                  0,
                 math::sin(gamma),   math::cos(gamma),                  0,
                                0,                  0,                  1;

        // Tiat-Bryan angle (X1,Y2,Z3) --> https://en.wikipedia.org/wiki/Euler_angles
        Ax *= Ay;
        Ax *= Az;
        A.swap(Ax);

        /*
              a11^2     ,a12^2      ,a13^2      ,2*a12*a11          ,2*a12*a13          ,2*a11*a13
              a21^2     ,a22^2      ,a23^2      ,2*a21*a22          ,2*a22*a23          ,2*a21*a23
              a31^2     ,a32^2      ,a33^2      ,2*a31*a32          ,2*a32*a33          ,2*a31*a33
              a11*a21   ,a12*a22    ,a13*a23    ,a11*a22+a12*a21    ,a12*a23+a13*a22    ,a11*a23+a13*a21
              a21*a31   ,a22*a32    ,a23*a33    ,a21*a32+a22*a31    ,a22*a33+a23*a32    ,a21*a33+a23*a31
              a11*a31   ,a12*a32    ,a13*a33    ,a11*a32+a31*a12    ,a12*a33+a13*a32    ,a11*a33+a31*a13
         */

        Phi.resize(6,6);
        Phi.setZero();

        Phi(0,0) = A(0,0)*A(0,0);
        Phi(0,1) = A(0,1)*A(0,1);
        Phi(0,2) = A(0,2)*A(0,2);
        Phi(0,3) = 2*A(0,0)*A(0,1);
        Phi(0,4) = 2*A(0,1)*A(0,2);
        Phi(0,5) = 2*A(0,0)*A(0,2);

        Phi(1,0) = A(1,0)*A(1,0);
        Phi(1,1) = A(1,1)*A(1,1);
        Phi(1,2) = A(1,2)*A(1,2);
        Phi(1,3) = 2*A(1,0)*A(1,1);
        Phi(1,4) = 2*A(1,1)*A(1,2);
        Phi(1,5) = 2*A(1,0)*A(1,2);

        Phi(2,0) = A(2,0)*A(2,0);
        Phi(2,1) = A(2,1)*A(2,1);
        Phi(2,2) = A(2,2)*A(2,2);
        Phi(2,3) = 2*A(2,0)*A(2,1);
        Phi(2,4) = 2*A(2,1)*A(2,2);
        Phi(2,5) = 2*A(2,0)*A(2,2);

        Phi(3,0) = A(0,0)*A(1,0);
        Phi(3,1) = A(0,1)*A(1,1);
        Phi(3,2) = A(0,2)*A(1,2);// !
        Phi(3,3) = A(0,0)*A(1,1)+A(0,1)*A(1,0);
        Phi(3,4) = A(0,1)*A(1,2)+A(0,2)*A(1,1);
        Phi(3,5) = A(0,0)*A(1,2)+A(0,2)*A(1,0);

        Phi(4,0) = A(1,0)*A(2,0);
        Phi(4,1) = A(1,1)*A(2,1);
        Phi(4,2) = A(1,2)*A(2,2);
        Phi(4,3) = A(1,0)*A(2,1)+A(1,1)*A(2,0);//////////////////
        Phi(4,4) = A(1,1)*A(2,2)+A(1,2)*A(2,1);
        Phi(4,5) = A(1,0)*A(2,2)+A(1,2)*A(2,0);

        Phi(5,0) = A(0,0)*A(2,0);
        Phi(5,1) = A(0,1)*A(2,1);
        Phi(5,2) = A(0,2)*A(2,2);// !
        Phi(5,3) = A(0,0)*A(2,1)+A(2,0)*A(0,1);
        Phi(5,4) = A(0,1)*A(2,2)+A(0,2)*A(2,1);
        Phi(5,5) = A(0,0)*A(2,2)+A(2,0)*A(0,2);
    }

protected:
    bool m_homogeneous;
    std::vector<gsMatrix<T>> m_G;

};

}
