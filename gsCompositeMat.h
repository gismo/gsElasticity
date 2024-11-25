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
        m_C.resize(6,6);
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

        //gsDebugVar(m_C);

        // Compliance matrix
        m_S.resize(6,6);
        m_S.setZero();
        
        // Diagonal components
        m_S(0,0) = 1/Exx; 
        m_S(1,1) = 1/Eyy; 
        m_S(2,2) = 1/Ezz; 
        m_S(3,3) = 1.0/Gxy; // Factor 2? ??????? CHECK!
        m_S(4,4) = 1.0/Gyz; // Factor 2?
        m_S(5,5) = 1.0/Gxz; // Factor 2?

        // Off-diagonal components
        m_S(0,1) = m_S(1,0) = -nuxy/Exx; // Gxy
        m_S(0,2) = m_S(2,0) = -nuxz/Exx; // Gxz
        m_S(1,2) = m_S(2,1) = -nuyz/Eyy;
        gsDebugVar(m_S);
    }

    gsCompositeMat( const T Exx,
                    const T Eyy,
                    const T Gxy,
                    const T nuxy,
                    const gsFunctionSet<T> & alpha,
                    const gsFunctionSet<T> & patches)
    :   
    m_alpha(&alpha),
    Base(&patches)
    // Base(&patches,nullptr,nullptr), in this case it will check the dimension of mp_def
    {
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

    // void computeHomogenizedStiffness(const index_t patch,  const gsMatrix<T> & u)
    // {
    //     // input a gsMultiPatch??????????
    //     // input numpatches?
        
    //     // Compute the parameter data (writes into m_data.mine().m_parMat)
    //     this->_computeParameterData(patch,u);

    //     // Voigt-size of the tensor
    //     const index_t sz = (d+1)*d/2;

    //     // Resize the matrix
    //     m_C_homogenized.resize(sz*sz,u.cols());
    //     m_C_homogenized.setZero();

    //     // const index_t k = numPatches; //?????
    //     const index_t k = 3; //?????
    //     //const real_t t_k = 1/k;
    //     const real_t t_k = 1./k;
    //     // const real_t Delta, term_44, term_55, Delta_k, sum_33, sum_44, sum_55;
    //     real_t Delta_k, sum_33 = 0.0, sum_44 = 0.0, sum_55 = 0.0;
    //     real_t term_44 = 0.0, term_55 = 0.0;

    //     gsMatrix<T> material_mat_1, material_mat_k; //???? I would need to do a resize?
        
    //     gsCompositeMat::eval_matrix_into(0,u,material_mat_1); // ????????????

    //     for (index_t j=0; j!=u.cols(); j++)
    //     {   
    //         // I would need to change the indices for the matrix!

    //         for (index i=0; i<k; i++) // loop over all layers ????
    //         {
    //             gsCompositeMat::eval_matrix_into(i,u,material_mat_k); // get material matrix of the container

    //             // First term of the sums
    //             m_C_homogenized(0,0) += t_k*material_mat_k(0,0);
    //             m_C_homogenized(0,1) += t_k*material_mat_k(0,1);
    //             m_C_homogenized(0,2) += t_k*material_mat_k(0,2);
    //             m_C_homogenized(1,1) += t_k*material_mat_k(1,1);
    //             m_C_homogenized(1,2) += t_k*material_mat_k(1,2);

    //             Delta_k = material_mat_k(3,3) * material_mat_k(4,4);
                
    //             // To calculate the other terms in the diagonal
    //             sum_33 += t_k / material_mat_k(2,2);
    //             sum_44 += t_k * material_mat_k(3,3) / Delta_k;
    //             sum_55 += t_k * material_mat_k(4,4) / Delta_k;

    //             // To calculate C_66_bar
    //             m_C_homogenized(5,5) += t_k*material_mat_k(5,5);

    //             // To calculate Delta
    //             term_44  += t_k * material_mat_k(3,3)/Delta_k;
    //             term_55  += t_k * material_mat_k(4,4)/Delta_k;
    //         }

    //         const real_t Delta = term_44*term_55;
    //         m_C_homogenized(2,2) = 1./sum_33;
    //         m_C_homogenized(3,3) = sum_44/Delta;
    //         m_C_homogenized(4,4) = sum_55/Delta;
            
    //         for (index q=2; q<k; q++) // loop over all layers
    //         {
    //             gsCompositeMat::eval_matrix_into(q,u,material_mat_k); // ????
    //             // Second term of the sums
    //             m_C_homogenized(0,0) += (material_mat_k(0,2)-m_C_homogenized(0,2))*t_k*(material_mat_1(0,2)-material_mat_k(0,2))/material_mat_k(2,2);
    //             m_C_homogenized(0,1) += (material_mat_k(0,2)-m_C_homogenized(0,2))*t_k*(material_mat_1(1,2)-material_mat_k(1,2))/material_mat_k(2,2);
    //             m_C_homogenized(0,2) += (material_mat_k(2,2)-m_C_homogenized(2,2))*t_k*(material_mat_1(0,2)-material_mat_k(0,2))/material_mat_k(2,2);
    //             m_C_homogenized(1,1) += (material_mat_k(1,2)-m_C_homogenized(1,2))*t_k*(material_mat_1(1,2)-material_mat_k(1,2))/material_mat_k(2,2);
    //             m_C_homogenized(1,2) += (material_mat_k(2,2)-m_C_homogenized(2,2))*t_k*(material_mat_1(1,2)-material_mat_k(1,2))/material_mat_k(2,2);
    //         } 
    //     }
    // }

    // void eval_stress_into(const index_t patch, const gsMatrix<T> & u, gsMatrix<T> & result) const
    // {
    //     // Compute the parameter data (writes into m_data.mine().m_parMat)
    //     this->_computeParameterData(patch,u);
    //     // Compute the strain
    //     gsMatrix<T> Eres,  stress_tensor, stress_voigt;

    //     Base::eval_strain_into(patch,u,Eres); // Eres in tensor notation

    //      // Resize the result
    //     result.resize(d*d,u.cols());
    //     stress_tensor.resize(d,d); 

    //     // Voigt-size of the tensor
    //     const index_t sz = (d+1)*d/2;

    //     gsMatrix<T> angles, Phi;
    //     angles = m_alpha->piece(patch).eval(u);
    //     // composite linear elasticity tensor

    //     for (index_t j=0; j!=u.cols(); j++)
    //     {
    //         gsMatrix<T> E_full = Eres.reshapeCol(j,d,d);

    //         gsVector<T> E_voigt;
    //         voigtStrain(E_voigt,E_full); //E_vect in voigt

    //         _makePhi(angles.col(patch),Phi);

    //         stress_voigt = Phi.transpose() * m_C * Phi * E_voigt; // this gives the stress in Voigt notation
    //         tensorStress(d,stress_voigt,stress_tensor);  // to get stress in tensor notation
    //         result.reshapeCol(j,d,d) = stress_tensor; // this gives the stress in tensor notation
    //     }
    // }

    void eval_matrix_into(const index_t patch, const gsMatrix<T> & u, gsMatrix<T> & result) const
    {
        // Compute the parameter data (writes into m_data.mine().m_parMat)
        this->_computeParameterData(patch,u);

        // Voigt-size of the tensor
        const index_t sz = (d+1)*d/2;

        // Resize the result
        result.resize(sz*sz,u.cols());

        gsMatrix<T> angles, Phi;
        angles = m_alpha->piece(patch).eval(u);
        // composite linear elasticity tensor

        //gsDebugVar(u);

        for (index_t j=0; j!=u.cols(); j++)
        {
           _makePhi(angles.col(patch),Phi);
           result.reshapeCol(j,sz,sz) = Phi.transpose() * m_C * Phi; // m_C is sz*sz?
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
    gsMatrix<T> m_S; // make stiffness matrix a member of the class-- constant over the whole domain
    gsMatrix<T> m_C_homogenized; // make stiffness matrix a member of the class-- constant over the whole domain


};

}
