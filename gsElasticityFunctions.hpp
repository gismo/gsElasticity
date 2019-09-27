/** @file gsElasticityFunctions.hpp

    @brief Provides useful classes derived from gsFunction which can be used
    for visualization or coupling.

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s):
        A.Shamanskiy (2016 - ...., TU Kaiserslautern)
*/

#pragma once

#include <gsElasticity/gsElasticityFunctions.h>

namespace gismo
{
template <class T>
void gsCauchyStressFunction<T>::eval_into(const gsMatrix<T> & u, gsMatrix<T> & result) const
{
    result.setZero(targetDim(),u.cols());

    // Each column of dispGrad contains components of the displacement Jacabian.
    // The order is u_xi u_eta v_xi v_eta (or similar in 3D)
    gsMatrix<T> dispGrad;
    m_displacement.patch(m_patch).deriv_into(u,dispGrad);

    for (index_t p = 0; p < u.cols(); ++p)
    {
        if (domainDim() == 2)
        {
            T s11 = m_lambda * (dispGrad(0,p) + dispGrad(3,p)) + 2 * m_mu * dispGrad(0,p);
            T s22 = m_lambda * (dispGrad(0,p) + dispGrad(3,p)) + 2 * m_mu * dispGrad(3,p);
            T s12 = m_mu * (dispGrad(1,p) + dispGrad(2,p));

            if (m_type == stress_type::von_mises)
                result(0,p) = math::sqrt(s11*s11 + s22*s22 - s11*s22 + 3.*s12*s12);
            if (m_type == stress_type::all_2D)
            {
                result(0,p) = s11;
                result(1,p) = s22;
                result(2,p) = s12;
            }
        }

        if (domainDim() == 3)
        {
            T s11 = m_lambda * (dispGrad(0,p) + dispGrad(4,p) + dispGrad(8,p)) + 2 * m_mu * dispGrad(0,p);
            T s22 = m_lambda * (dispGrad(0,p) + dispGrad(4,p) + dispGrad(8,p)) + 2 * m_mu * dispGrad(4,p);
            T s33 = m_lambda * (dispGrad(0,p) + dispGrad(4,p) + dispGrad(8,p)) + 2 * m_mu * dispGrad(8,p);
            T s12 = m_mu * (dispGrad(1,p) + dispGrad(3,p));
            T s13 = m_mu * (dispGrad(2,p) + dispGrad(6,p));
            T s23 = m_mu * (dispGrad(5,p) + dispGrad(7,p));
            if (m_type == stress_type::von_mises)
                result(0,p) = math::sqrt(0.5 * ( (s11-s22)*(s11-s22) + (s11-s33)*(s11-s33) + (s22-s33)*(s22-s33) +
                                                 6 * (s12*s12 + s13*s13 + s23*s23) ) );
            if (m_type == stress_type::normal_3D)
            {
                result(0,p) = s11;
                result(1,p) = s22;
                result(2,p) = s33;
            }
            if (m_type == stress_type::shear_3D)
            {
                result(0,p) = s12;
                result(1,p) = s13;
                result(2,p) = s23;
            }
        }
    }
}

template <class T>
void gsDetFunction<T>::eval_into(const gsMatrix<T> & u, gsMatrix<T> & result) const
{
    result.resize(1,u.cols());

    gsMapData<T> mappingData;
    mappingData.points = u;
    mappingData.flags = NEED_DERIV;
    m_geo.patch(m_patch).computeMap(mappingData);

    for (index_t i = 0; i < u.cols(); ++i)
        result(0,i) = mappingData.jacobian(i).determinant();
}

template <class T>
void gsFsiLoad<T>::eval_into(const gsMatrix<T> & u, gsMatrix<T> & result) const
{
    result.setZero(targetDim(),u.cols());
    // mapping points back to the parameter space via the reference configuration
    gsMatrix<T> paramPoints;
    m_geo.patch(m_patchGeo).invertPoints(u,paramPoints);
    // evaluate reference geometry mapping at the param points
    // NEED_GRAD_TRANSFORM for velocity gradients transformation from parametric to reference domain
    gsMapData<T> mdGeo(NEED_GRAD_TRANSFORM);
    mdGeo.points = paramPoints;
    m_geo.patch(m_patchGeo).computeMap(mdGeo);
    // evaluate velocity at the param points
    // NEED_DERIV for velocity gradients
    gsMapData<T> mdVel(NEED_DERIV);
    mdVel.points = paramPoints;
    m_vel.patch(m_patchVP).computeMap(mdVel);
    // evaluate pressure at the quad points
    gsMatrix<T> pressureValues;
    m_pres.patch(m_patchVP).eval_into(paramPoints,pressureValues);
    // evaluate ALE dispacement at the param points
    // NEED_DERIV for gradients
    gsMapData<T> mdALE(NEED_DERIV);
    mdALE.points = paramPoints;
    m_ale.patch(m_patchGeo).computeMap(mdALE);

    for (index_t p = 0; p < paramPoints.cols(); ++p)
    {
        // transform velocity gradients from parametric to reference
        gsMatrix<T> physGradVel = mdVel.jacobian(p)*(mdGeo.jacobian(p).cramerInverse());
        // ALE jacobian (identity + physical displacement gradient)
        gsMatrix<T> physJacALE = gsMatrix<T>::Identity(targetDim(),targetDim()) +
                mdALE.jacobian(p)*(mdGeo.jacobian(p).cramerInverse());
        // inverse ALE jacobian
        gsMatrix<T> invJacALE = physJacALE.cramerInverse();
        // ALE stress tensor
        gsMatrix<T> sigma = pressureValues.at(p)*gsMatrix<T>::Identity(targetDim(),targetDim())
                            - m_density*m_viscosity*(physGradVel*invJacALE +
                                           invJacALE.transpose()*physGradVel.transpose());
        // stress tensor pull back
        gsMatrix<T> sigmaALE = physJacALE.determinant()*sigma*(invJacALE.transpose());

        // normal length is the local measure
        gsVector<T> normal;
        outerNormal(mdGeo,p,m_sideGeo,normal);
        result.col(p) = sigmaALE * normal / normal.norm();
    }
}

} // namespace gismo ends
