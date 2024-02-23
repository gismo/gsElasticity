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

/*
    To do:
        - G will be a gsFunctionSet or gsFunction!
        - alpha, beta, gamma --> in one gsFunctionSet m_angles

        -   if m_angles.rows()==1 --> 2D elasticity
            if m_angles.rows()==3 --> 3D elasticity
        -   m_G.rows()==4 OR m_G.rows()==9
        - _makePhi will have input gsMatrix<T>. If 3 rows, then phi3D and if 1 row, then phi2D

*/

    gsCompositeMaterial(    const gsFunctionSet<T> & G,
                            const gsFunctionSet<T> & alpha)
    :   m_G(&G),
        m_alpha(&alpha)
    {

    }

    mutable util::gsThreaded<gsMatrix<T> > C, Ctemp;

    /// \a u points, \a k patch index
    void eval_matrix_into(const gsMatrix<T>& u, gsMatrix<T>& result, index_t k = 0) const
    {
        const index_t sz = (u.rows()+1)*u.rows()/2;
        result.resize(sz*sz,u.cols());


        gsMatrix<T> angles, Phi, G;
        angles = m_alpha->piece(k).eval(u);
        G = m_G->piece(k).eval(u);
        // composite linear elasticity tensor
        for (index_t j=0; j!=u.cols(); j++)
        {
           _makePhi(angles.col(k),Phi);
           result.reshapeCol(j,sz,sz) = Phi.transpose() * G.reshapeCol(j,sz,sz) * Phi;
        }
    }

        /// \a u points, \a k patch index
        void eval_vector_into(const gsMatrix<T>& u, gsMatrix<T>& result, index_t k = 0) const
    {
        GISMO_NO_IMPLEMENTATION;
    }

protected:
    /**
     * @brief      Makes the rotation matrix (Tiat-Bryan angle)
     *
     * @param[in]  alpha  The rotation around the x axis
     * @param[in]  beta   The rotation around the y axis
     * @param[in]  gamma  The rotation around the z axis
     * @param[in]  Phi    The rotation matrix
     *
     */
    void _makePhi(  const gsMatrix<T> angles,
                    gsMatrix<T> & Phi) const
    {
        if (angles.rows()==1) // 2d elasticity
        {
            Phi.resize(3,3);
            Phi.setZero();

            const T & alpha = angles(0,0);

            // Make transformation matrix
            Phi(0,0) = Phi(1,1) = math::pow(math::cos(alpha),2);
            Phi(0,1) = Phi(1,0) = math::pow(math::sin(alpha),2);
            Phi(2,0) = Phi(0,2) = Phi(2,1) = Phi(1,2) = math::sin(alpha) * math::cos(alpha);
            Phi(2,0) *= -2.0;
            Phi(2,1) *= 2.0;
            Phi(1,2) *= -1.0;
            Phi(2,2) = math::pow(math::cos(alpha),2) - math::pow(math::sin(alpha),2);
        }
        else if (angles.rows()==3)
        {
            Phi.resize(6,6);
            Phi.setZero();

            const T & alpha = angles(0,0);
            const T & beta  = angles(1,0);
            const T & gamma = angles(2,0);

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
                Phi:
                  a11^2     ,a12^2      ,a13^2      ,2*a12*a11          ,2*a12*a13          ,2*a11*a13
                  a21^2     ,a22^2      ,a23^2      ,2*a21*a22          ,2*a22*a23          ,2*a21*a23
                  a31^2     ,a32^2      ,a33^2      ,2*a31*a32          ,2*a32*a33          ,2*a31*a33
                  a11*a21   ,a12*a22    ,a13*a23    ,a11*a22+a12*a21    ,a12*a23+a13*a22    ,a11*a23+a13*a21
                  a21*a31   ,a22*a32    ,a23*a33    ,a21*a32+a22*a31    ,a22*a33+a23*a32    ,a21*a33+a23*a31
                  a11*a31   ,a12*a32    ,a13*a33    ,a11*a32+a31*a12    ,a12*a33+a13*a32    ,a11*a33+a31*a13
             */

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
            Phi(3,2) = A(0,2)*A(1,2);
            Phi(3,3) = A(0,0)*A(1,1)+A(0,1)*A(1,0);
            Phi(3,4) = A(0,1)*A(1,2)+A(0,2)*A(1,1);
            Phi(3,5) = A(0,0)*A(1,2)+A(0,2)*A(1,0);

            Phi(4,0) = A(1,0)*A(2,0);
            Phi(4,1) = A(1,1)*A(2,1);
            Phi(4,2) = A(1,2)*A(2,2);
            Phi(4,3) = A(1,0)*A(2,1)+A(1,1)*A(2,0);
            Phi(4,4) = A(1,1)*A(2,2)+A(1,2)*A(2,1);
            Phi(4,5) = A(1,0)*A(2,2)+A(1,2)*A(2,0);

            Phi(5,0) = A(0,0)*A(2,0);
            Phi(5,1) = A(0,1)*A(2,1);
            Phi(5,2) = A(0,2)*A(2,2);
            Phi(5,3) = A(0,0)*A(2,1)+A(2,0)*A(0,1);
            Phi(5,4) = A(0,1)*A(2,2)+A(0,2)*A(2,1);
            Phi(5,5) = A(0,0)*A(2,2)+A(2,0)*A(0,2);
        }
        else
            GISMO_ERROR("angles have wrong dimension ("<<angles.rows()<<")");
    }

protected:
    const gsFunctionSet<T> * m_G;
    const gsFunctionSet<T> * m_alpha;

};

}
