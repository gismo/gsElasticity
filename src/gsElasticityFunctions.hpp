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

#include <gsElasticity/src/gsElasticityFunctions.h>
#include <gsCore/gsFuncData.h>
#include <gsAssembler/gsAssembler.h>

namespace gismo
{

template <class T>
void gsCauchyStressFunction<T>::linearElastic(const gsMatrix<T> & u, gsMatrix<T> & result) const
{
    result.setZero(targetDim(),outputCols(u.cols()));
    // evaluating the fields
    gsMapData<T> mdGeo(NEED_GRAD_TRANSFORM);
    mdGeo.points = u;
    m_geometry->patch(m_patch).computeMap(mdGeo);
    gsMapData<T> mdDisp(NEED_DERIV);
    mdDisp.points = u;
    m_displacement->patch(m_patch).computeMap(mdDisp);
    // define temporary matrices here for efficieny
    gsMatrix<T> I = gsMatrix<T>::Identity(m_dim,m_dim);
    gsMatrix<T> sigma,eps,dispGrad;
    // material parameters
    T YM = m_options.getReal("YoungsModulus");
    T PR = m_options.getReal("PoissonsRatio");
    T lambda = YM * PR / ( ( 1. + PR ) * ( 1. - 2. * PR ) );
    T mu     = YM / ( 2. * ( 1. + PR ) );

    for (index_t q = 0; q < u.cols(); ++q)
    {
        // linear strain tensor eps = (gradU+gradU^T)/2
        if (mdGeo.jacobian(q).determinant() <= 0)
            gsInfo << "Invalid domain parametrization: J = " << mdGeo.jacobian(q).determinant() <<
                      " at point (" << u.col(q).transpose() << ") of patch " << m_patch << std::endl;
        if (abs(mdGeo.jacobian(q).determinant()) > 1e-20)
           dispGrad = mdDisp.jacobian(q)*(mdGeo.jacobian(q).cramerInverse());
        else
            dispGrad = gsMatrix<T>::Zero(m_dim,m_dim);
        eps = (dispGrad + dispGrad.transpose())/2;
        // Cauchy stress tensor
        sigma = lambda*eps.trace()*I + 2*mu*eps;
        saveStress(sigma,result,q);
    }
}

template <class T>
void gsCauchyStressFunction<T>::nonLinearElastic(const gsMatrix<T> & u, gsMatrix<T> & result) const
{
    result.setZero(targetDim(),outputCols(u.cols()));
    // evaluating the fields
    gsMapData<T> mdGeo(NEED_GRAD_TRANSFORM);
    mdGeo.points = u;
    m_geometry->patch(m_patch).computeMap(mdGeo);
    gsMapData<T> mdDisp(NEED_DERIV);
    mdDisp.points = u;
    m_displacement->patch(m_patch).computeMap(mdDisp);
    // define temporary matrices here for efficieny
    gsMatrix<T> I = gsMatrix<T>::Identity(m_dim,m_dim);
    gsMatrix<T> S,sigma,F,C,E;
    // material parameters
    T YM = m_options.getReal("YoungsModulus");
    T PR = m_options.getReal("PoissonsRatio");
    T lambda = YM * PR / ( ( 1. + PR ) * ( 1. - 2. * PR ) );
    T mu     = YM / ( 2. * ( 1. + PR ) );

    for (index_t q = 0; q < u.cols(); ++q)
    {
        // deformation gradient F = I + gradU*gradGeo^-1
        if (mdGeo.jacobian(q).determinant() <= 0)
            gsInfo << "Invalid domain parametrization: J = " << mdGeo.jacobian(q).determinant() <<
                      " at point (" << u.col(q).transpose() << ") of patch " << m_patch << std::endl;
        if (abs(mdGeo.jacobian(q).determinant()) > 1e-20)
            F = I + mdDisp.jacobian(q)*(mdGeo.jacobian(q).cramerInverse());
        else
            F = I;
        T J = F.determinant();
        if (J <= 0)
            gsInfo << "Invalid displacement field: J = " << J <<
                      " at point (" << u.col(q).transpose() << ") of patch " << m_patch << std::endl;
        // Second Piola-Kirchhoff stress tensor
        if (material_law::law(m_options.getInt("MaterialLaw")) == material_law::saint_venant_kirchhoff)
        {
            // Green-Lagrange strain tensor E = 0.5(F^T*F-I)
            E = (F.transpose() * F - I)/2;
            S = lambda*E.trace()*I + 2*mu*E;
        }
        if (material_law::law(m_options.getInt("MaterialLaw")) == material_law::neo_hooke_ln)
        {
            // Right Cauchy Green strain, C = F'*F
            C = F.transpose() * F;
            S = (lambda*log(J)-mu)*(C.cramerInverse()) + mu*I;
        }
        if (material_law::law(m_options.getInt("MaterialLaw")) == material_law::neo_hooke_quad)
        {
            // Right Cauchy Green strain, C = F'*F
            C = F.transpose() * F;
            S = (lambda*(J*J-1)/2-mu)*(C.cramerInverse()) + mu*I;
        }
        // transformation to Cauchy stress
        sigma = F*S*F.transpose()/J;
        saveStress(sigma,result,q);
    }
}

template <class T>
void gsCauchyStressFunction<T>::mixedLinearElastic(const gsMatrix<T> & u, gsMatrix<T> & result) const
{
    result.setZero(targetDim(),outputCols(u.cols()));
    // evaluating the fields
    gsMapData<T> mdGeo(NEED_GRAD_TRANSFORM);
    mdGeo.points = u;
    m_geometry->patch(m_patch).computeMap(mdGeo);
    gsMapData<T> mdDisp(NEED_DERIV);
    mdDisp.points = u;
    m_displacement->patch(m_patch).computeMap(mdDisp);
    gsMatrix<T> presVals;
    m_pressure->patch(m_patch).eval_into(u,presVals);
    // define temporary matrices here for efficieny
    gsMatrix<T> I = gsMatrix<T>::Identity(m_dim,m_dim);
    gsMatrix<T> sigma,eps,dispGrad;
    // material parameters
    T YM = m_options.getReal("YoungsModulus");
    T PR = m_options.getReal("PoissonsRatio");
    T mu     = YM / ( 2. * ( 1. + PR ) );

    for (index_t q = 0; q < u.cols(); ++q)
    {
        // linear strain tensor eps = (gradU+gradU^T)/2
        if (mdGeo.jacobian(q).determinant() <= 0)
            gsInfo << "Invalid domain parametrization: J = " << mdGeo.jacobian(q).determinant() <<
                      " at point (" << u.col(q).transpose() << ") of patch " << m_patch << std::endl;
        if (abs(mdGeo.jacobian(q).determinant()) > 1e-20)
           dispGrad = mdDisp.jacobian(q)*(mdGeo.jacobian(q).cramerInverse());
        else
            dispGrad = gsMatrix<T>::Zero(m_dim,m_dim);
        eps = (dispGrad + dispGrad.transpose())/2;
        // Cauchy stress tensor
        sigma = presVals.at(q)*I + 2*mu*eps;
        saveStress(sigma,result,q);
    }
}

template <class T>
void gsCauchyStressFunction<T>::mixedNonLinearElastic(const gsMatrix<T> & u, gsMatrix<T> & result) const
{
    result.setZero(targetDim(),outputCols(u.cols()));
    // evaluating the fields
    gsMapData<T> mdGeo(NEED_GRAD_TRANSFORM);
    mdGeo.points = u;
    m_geometry->patch(m_patch).computeMap(mdGeo);
    gsMapData<T> mdDisp(NEED_DERIV);
    mdDisp.points = u;
    m_displacement->patch(m_patch).computeMap(mdDisp);
    gsMatrix<T> presVals;
    m_pressure->patch(m_patch).eval_into(u,presVals);
    // define temporary matrices here for efficieny
    gsMatrix<T> I = gsMatrix<T>::Identity(m_dim,m_dim);
    gsMatrix<T> S,sigma,F,C,E;
    // material parameters
    T YM = m_options.getReal("YoungsModulus");
    T PR = m_options.getReal("PoissonsRatio");
    T mu     = YM / ( 2. * ( 1. + PR ) );

    for (index_t q = 0; q < u.cols(); ++q)
    {
        if (mdGeo.jacobian(q).determinant() <= 0)
            gsInfo << "Invalid domain parametrization: J = " << mdGeo.jacobian(q).determinant() <<
                      " at point (" << u.col(q).transpose() << ") of patch " << m_patch << std::endl;
        // deformation gradient F = I + gradU*gradGeo^-1
        if (abs(mdGeo.jacobian(q).determinant()) > 1e-20)
            F = I + mdDisp.jacobian(q)*(mdGeo.jacobian(q).cramerInverse());
        else
            F = I;
        T J = F.determinant();
        if (J <= 0)
            gsInfo << "Invalid displacement field: J = " << J <<
                      " at point (" << u.col(q).transpose() << ") of patch " << m_patch << std::endl;
        // Second Piola-Kirchhoff stress tensor
        C = F.transpose() * F;
        S = (presVals.at(q)-mu)*(C.cramerInverse()) + mu*I;
        // transformation to Cauchy stress
        sigma = F*S*F.transpose()/J;
        saveStress(sigma,result,q);
    }
}

template <class T>
void gsCauchyStressFunction<T>::saveStress(const gsMatrix<T> & S, gsMatrix<T> & result, index_t q) const
{
    switch (m_type)
    {
    case stress_components::von_mises :
    {
        if (m_geometry->parDim() == 2)
            result(0,q) = sqrt( S(0,0)*S(0,0) + S(1,1)*S(1,1) - S(0,0)*S(1,1) + 3*S(0,1)*S(0,1) );
        if (m_geometry->parDim() == 3)
            result(0,q) = sqrt(0.5*( pow(S(0,0)-S(1,1),2) + pow(S(0,0)-S(2,2),2) + pow(S(1,1)-S(2,2),2) +
                                     6 * (pow(S(0,1),2) + pow(S(0,2),2) + pow(S(1,2),2) ) ) );
        return;
    }
    case stress_components::all_2D_vector : result(0,q) = S(0,0); result(1,q) = S(1,1); result(2,q) = S(0,1); return;
    case stress_components::all_2D_matrix : result.middleCols(q*2,2) = S; return;
    case stress_components::normal_3D_vector : result(0,q) = S(0,0); result(1,q) = S(1,1); result(2,q) = S(2,2); return;
    case stress_components::shear_3D_vector : result(0,q) = S(0,1); result(1,q) = S(0,2); result(2,q) = S(1,2); return;
    case stress_components::all_3D_matrix : result.middleCols(q*3,3) = S; return;
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

    gsMatrix<T> I  = gsMatrix<T>::Identity(targetDim(),targetDim());
    for (index_t p = 0; p < paramPoints.cols(); ++p)
    {
        // transform velocity gradients from parametric to reference
        gsMatrix<T> physGradVel = mdVel.jacobian(p)*(mdGeo.jacobian(p).cramerInverse());
        // ALE jacobian (identity + physical displacement gradient)
        gsMatrix<T> physJacALE = I + mdALE.jacobian(p)*(mdGeo.jacobian(p).cramerInverse());
        // inverse ALE jacobian
        gsMatrix<T> invJacALE = physJacALE.cramerInverse();
        // ALE stress tensor
        gsMatrix<T> sigma = pressureValues.at(p)*I
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
