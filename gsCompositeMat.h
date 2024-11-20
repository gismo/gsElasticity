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

#include <gsElasticity/gsMaterialBaseDim.h>
#include <gsElasticity/gsVisitorElUtils.h>
#include <gsUtils/gsThreaded.h>

namespace gismo
{

template <short_t d, class T>
class gsCompositeMat : public gsMaterialBaseDim<d,T>
{

public:
    using Base = gsMaterialBaseDim<d,T>;

    gsCompositeMat( const T Exx,
                    const T Eyy,
                    const T Ezz,
                    const T Gxy,
                    const T Gyz,
                    const T Gxz,
                    const T nuxy,
                    const T nuyz,
                    const T nuxz,
                    const gsFunctionSet<T> & alpha,
                    const gsFunctionSet<T> & patches)
    :   
    m_alpha(&alpha),
    Base(&patches)
    {
        // gsMatrix<T> m_C(6,6);

        m_C.resize(6,6);
        std::cout << "m_C rows: " << m_C.rows() << ", cols: " << m_C.cols() << std::endl;

        m_C.setZero();
        T D = (1-nuxy*nuxy-nuyz*nuyz-nuxz*nuxz-2*nuxy*nuyz*nuxz) / (Exx*Eyy*Ezz); // determinant of the compliance matrix (S)
        
        m_C(0,0) = (1-nuyz*nuyz) / (Eyy*Ezz); // Gxx
        m_C(1,1) = (1-nuxz*nuxz) / (Exx*Ezz); // Gyy
        m_C(2,2) = (1-nuxy*nuxy) / (Exx*Eyy); // Gzz

        // Gij = nuij + nuik*nukj / Ejj*Ekk
        m_C(0,1) = m_C(1,0) = (nuxy+nuxz*nuyz) / (Eyy*Ezz); // Gxy
        m_C(0,2) = m_C(2,0) = (nuxz+nuxy*nuyz) / (Eyy*Ezz); // Gxz
        m_C(1,2) = m_C(2,1) = (nuyz+nuxy*nuxz) / (Exx*Ezz); // Gyz

        m_C *= 1.0/D;

        m_C(3,3) = Gxy; // Factor 2? ??????? CHECK!
        m_C(4,4) = Gyz; // Factor 2?
        m_C(5,5) = Gxz; // Factor 2?

        // Set parameter becuase they could vary (heterogeneous materials? short fiber?)
        // Not sure about this ?
        // this -> setParameter(0,m_C(0,0)); //????? comprobar!
        // this -> setParameter(1,m_C(0,1));
        // this -> setParameter(2,m_C(0,2));
        // this -> setParameter(3,m_C(1,1));
        // this -> setParameter(4,m_C(1,2));
        // this -> setParameter(5,m_C(2,2));
        // this -> setParameter(6,m_C(3,3));
        // this -> setParameter(7,m_C(4,4));
        // this -> setParameter(8,m_C(5,5));

        gsDebugVar(m_C);
    }

    gsCompositeMat( const T Exx,
                    const T Eyy,
                    const T Gxy,
                    const T nuxy,
                    const gsFunctionSet<T> & alpha,
                    const gsFunctionSet<T> & patches)
    :   
    m_alpha(&alpha),
    Base(&patches,nullptr,nullptr)
    // Base(&patches,nullptr,nullptr), in this case it will check the dimension of mp_def
    {
        // gsMatrix<T> m_C(3,3);
        m_C.resize(3,3);
        m_C.setZero();

        m_C(0,0) = Exx*nuxy / ((1+nuxy)*(1-2*nuxy)) + Exx / ((1+nuxy)); // Gxx
        m_C(1,1) = Eyy*nuxy / ((1+nuxy)*(1-2*nuxy)) + Eyy / ((1+nuxy)); // Gxx
        m_C(0,1) = m_C(1,0) = (nuxy*Exx) / ((1+nuxy)*(1-2*nuxy)); // Gxy

        m_C(2,2) = Gxy;

        // this -> setParameter(0,m_C(0,0));
        // this -> setParameter(1,m_C(0,1));
        // this -> setParameter(2,m_C(1,0));
        // this -> setParameter(3,m_C(1,1));
        // this -> setParameter(4,m_C(2,2));
    }

    void eval_matrix_into(const index_t patch, const gsMatrix<T> & u, gsMatrix<T> & result) const
    {
        gsDebugVar(patch);
        gsDebugVar(u);
        // Compute the parameter data (writes into m_data.mine().m_parMat)
        this->_computeParameterData(patch,u);
 
        // Voigt-size of the tensor
        const index_t sz = (d+1)*d/2;

        // Resize the result
        result.resize(sz*sz,u.cols());

        gsMatrix<T> angles, Phi;
        angles = m_alpha->piece(patch).eval(u);
        // composite linear elasticity tensor

        gsInfo<<"Matrix unrotated\n";
        gsDebugVar(m_C);
        gsInfo<<"================\n";

        for (index_t j=0; j!=u.cols(); j++)
        {
           _makePhi(angles.col(patch),Phi);
        //    result.reshapeCol(j,sz,sz) = Phi.transpose() * G.reshapeCol(j,sz,sz) * Phi;
           result.reshapeCol(j,sz,sz) = Phi.transpose() * m_C * Phi; // m_C is sz*sz?
           gsInfo<<"Matrix per point\n";
           gsDebugVar(Phi.transpose() * m_C * Phi);
           gsInfo<<"================\n";
        }
        //gsDebugVar(result);
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
    void _makePhi(  const gsMatrix<T> angles, gsMatrix<T> & Phi) const
    {
        if (angles.rows()==1) // 2D elasticity
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
    using Base::m_data;
    const gsFunctionSet<T> * m_alpha;
    gsMatrix<T> m_C; // make stiffness matrix a member of the class-- constant over the whole domain


};

}
