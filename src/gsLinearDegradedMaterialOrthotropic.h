/** @file gsLinearDegradedMaterialOrthotropic.h

    @brief Linear *orthotropic* elastic material with a scalar phase-field
           degradation g(d) = (1-d)^2 applied to the (full) elastic energy.

    This is the anisotropic counterpart of gsLinearDegradedMaterial. Following
    "Calibration and validation of a phase-field model of brittle fracture",
    the anisotropy enters ONLY through the elastic strain-energy density (the
    orthotropic Hooke matrix). The fracture toughness / phase-field part is
    handled elsewhere and remains fully isotropic.

    No tension/compression (vol/dev) split is used here: the degradation
    function multiplies the complete orthotropic stress / stiffness, i.e.
        S   = g(d) * C_ortho : E
        C   = g(d) * C_ortho
        Psi = 0.5 * E : C_ortho : E      (undegraded crack-driving energy)

    Parameter layout (set via gsMaterialBase::setParameter):
        3D: 0:E1 1:E2 2:E3 3:G12 4:G13 5:G23 6:nu12 7:nu13 8:nu23 9:damage
        2D (plane stress): uses 0:E1 1:E2 3:G12 6:nu12 and 9:damage

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
#include <gsElasticity/gsMaterialUtils.h>
#include <gsCore/gsConstantFunction.h>

namespace gismo
{

/**
 * @brief Orthotropic linear elastic material degraded by a phase field.
 * @ingroup Elasticity
 * @tparam T Real type
 */
template <class T>
class gsLinearDegradedMaterialOrthotropic : public gsMaterialBase<T>
{

public:
    using Base = gsMaterialBase<T>;

    GISMO_CLONE_FUNCTION(gsLinearDegradedMaterialOrthotropic);

    // Parameter indices
    enum
    {
        iE1 = 0, iE2 = 1, iE3 = 2,
        iG12 = 3, iG13 = 4, iG23 = 5,
        inu12 = 6, inu13 = 7, inu23 = 8,
        iDamage = 9
    };

    /**
     * @brief Constructor with constant orthotropic parameters (no density).
     * @param E1,E2,E3    Young's moduli in the material directions
     * @param G12,G13,G23 Shear moduli
     * @param nu12,nu13,nu23 Poisson's ratios (convention nu_ij/E_i = nu_ji/E_j)
     * @param damage      Phase-field damage field
     * @param d           Domain dimension (2 or 3)
     */
    gsLinearDegradedMaterialOrthotropic(const T E1,  const T E2,  const T E3,
                                        const T G12, const T G13, const T G23,
                                        const T nu12,const T nu13,const T nu23,
                                        const gsFunctionSet<T> & damage,
                                        short_t d)
    :
    Base()
    {
        this->setParameter(iE1 ,gsConstantFunction<T>(E1 ,d));
        this->setParameter(iE2 ,gsConstantFunction<T>(E2 ,d));
        this->setParameter(iE3 ,gsConstantFunction<T>(E3 ,d));
        this->setParameter(iG12,gsConstantFunction<T>(G12,d));
        this->setParameter(iG13,gsConstantFunction<T>(G13,d));
        this->setParameter(iG23,gsConstantFunction<T>(G23,d));
        this->setParameter(inu12,gsConstantFunction<T>(nu12,d));
        this->setParameter(inu13,gsConstantFunction<T>(nu13,d));
        this->setParameter(inu23,gsConstantFunction<T>(nu23,d));
        this->setParameter(iDamage,damage,true);
    }

    /**
     * @brief Constructor with constant orthotropic parameters and density.
     */
    gsLinearDegradedMaterialOrthotropic(const T E1,  const T E2,  const T E3,
                                        const T G12, const T G13, const T G23,
                                        const T nu12,const T nu13,const T nu23,
                                        const T rho,
                                        const gsFunctionSet<T> & damage,
                                        short_t d)
    :
    gsLinearDegradedMaterialOrthotropic(E1,E2,E3,G12,G13,G23,nu12,nu13,nu23,damage,d)
    {
        this->setDensity(gsConstantFunction<T>(rho,d));
    }

    /// See \ref gsMaterialBase.h for more details
    void compute_stress_into(const gsMaterialData<T> & data, gsMatrix<T> & Sresult) const override
    {
        this->eval_stress_into(data, Sresult);
    }

    static void eval_stress_into(const gsMaterialData<T> & data, gsMatrix<T> & Sresult)
    {
        const short_t dim = data.dim;
        const index_t N = data.size;

        Sresult.resize(dim*dim,N);

        gsMatrix<T> E_vec;
        calculate_voigt_strain(data.strain,dim,E_vec); // engineering strain in Voigt

        for (index_t i=0; i!=N; i++)
        {
            const T damage = data.parameters[iDamage](0,i);
            const T omega  = math::pow((1. - damage),2);

            gsMatrix<T> C = _stiffnessVoigt(dim,data,i);

            gsMatrix<T> Svec = omega * (C * E_vec.col(i));
            gsMatrix<T> Smat;
            tensorStress(dim,Svec,Smat);
            Sresult.reshapeCol(i,dim,dim) = Smat;
        }
    }

    /// See \ref gsMaterialBase.h for more details
    void compute_matrix_into(const gsMaterialData<T> & data, gsMatrix<T> & Cresult) const override
    {
        this->eval_matrix_into(data, Cresult);
    }

    static void eval_matrix_into(const gsMaterialData<T> & data, gsMatrix<T> & Cresult)
    {
        const short_t dim = data.dim;
        const index_t N = data.size;
        const index_t sz = (dim+1)*dim/2;

        Cresult.resize(sz*sz,N);

        for (index_t i=0; i!=N; i++)
        {
            const T damage = data.parameters[iDamage](0,i);
            const T omega  = math::pow((1. - damage),2);

            gsMatrix<T> C = _stiffnessVoigt(dim,data,i);
            Cresult.reshapeCol(i,sz,sz) = omega * C;
        }
    }

    /// See \ref gsMaterialBase.h for more details
    void compute_energy_into(const gsMaterialData<T> & data, gsMatrix<T> & Presult) const override
    {
        this->eval_energy_into(data, Presult);
    }

    static void eval_energy_into(const gsMaterialData<T> & data, gsMatrix<T> & Presult)
    {
        const short_t dim = data.dim;
        const index_t N = data.size;

        Presult.resize(1,N);

        gsMatrix<T> E_vec;
        calculate_voigt_strain(data.strain,dim,E_vec); // engineering strain in Voigt

        for (index_t i=0; i!=N; i++)
        {
            // Undegraded (full) orthotropic elastic energy: drives the phase field.
            gsMatrix<T> C = _stiffnessVoigt(dim,data,i);
            Presult(0,i) = 0.5 * (E_vec.col(i).transpose() * C * E_vec.col(i)).value();
        }
    }

protected:

    /**
     * @brief Builds the orthotropic Hooke (stiffness) matrix in Voigt notation
     *        at evaluation point i, using the parameters stored in @a data.
     *
     * 3D Voigt order: [xx, yy, zz, xy, yz, xz] with engineering shear strains.
     * 2D Voigt order: [xx, yy, xy] (plane stress).
     */
    static gsMatrix<T> _stiffnessVoigt(const short_t dim,
                                       const gsMaterialData<T> & data,
                                       const index_t i)
    {
        const index_t sz = (dim+1)*dim/2;
        gsMatrix<T> C = gsMatrix<T>::Zero(sz,sz);

        const T E1  = data.parameters[iE1 ](0,i);
        const T E2  = data.parameters[iE2 ](0,i);
        const T nu12= data.parameters[inu12](0,i);
        const T G12 = data.parameters[iG12](0,i);

        if (dim==3)
        {
            const T E3  = data.parameters[iE3 ](0,i);
            const T G13 = data.parameters[iG13](0,i);
            const T G23 = data.parameters[iG23](0,i);
            const T nu13= data.parameters[inu13](0,i);
            const T nu23= data.parameters[inu23](0,i);

            // Compliance matrix S (engineering Voigt), then C = S^{-1}
            gsMatrix<T> S = gsMatrix<T>::Zero(6,6);
            S(0,0) = 1./E1;  S(1,1) = 1./E2;  S(2,2) = 1./E3;
            S(0,1) = S(1,0) = -nu12/E1;
            S(0,2) = S(2,0) = -nu13/E1;
            S(1,2) = S(2,1) = -nu23/E2;
            // Voigt shear indices: 3->xy, 4->yz, 5->xz
            S(3,3) = 1./G12; S(4,4) = 1./G23; S(5,5) = 1./G13;
            C = S.inverse();
        }
        else if (dim==2)
        {
            // Plane stress orthotropic
            const T nu21 = nu12 * E2 / E1;
            const T den  = 1. - nu12*nu21;
            C(0,0) = E1/den;        C(1,1) = E2/den;
            C(0,1) = C(1,0) = nu12*E2/den;
            C(2,2) = G12;
        }
        else
            GISMO_ERROR("gsLinearDegradedMaterialOrthotropic: dimension must be 2 or 3.");

        return C;
    }
};

}
