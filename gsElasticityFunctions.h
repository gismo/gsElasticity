/** @file gsElasticityFunctions.h

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

#include <gsCore/gsFunction.h>
#include <gsElasticity/gsBaseUtils.h>
#include <gsCore/gsBoundary.h>

namespace gismo
{


/** @brief Compute Cauchy stresses for a previously computed/defined displacement field.
 *         Can be pushed into gsPiecewiseFunction to construct gsField for visualization in Paraview.
*/
template <class T>
class gsCauchyStressFunction : public gsFunction<T>
{
public:

    gsCauchyStressFunction(const gsMultiPatch<T> & geometry,
                           const gsMultiPatch<T> & displacement,
                           index_t patch, stress_type::type type,
                           T lambda, T mu)
        : m_geometry(geometry),
          m_displacement(displacement),
          m_patch(patch),
          m_type(type),
          m_lambda(lambda),
          m_mu(mu)
    {}

    virtual short_t domainDim() const
    {
        return m_displacement.dim();
    }

    virtual short_t targetDim() const
    {
        switch(m_type)
        {
            case stress_type::von_mises:
                return 1;
            case stress_type::all_2D:
                return 3;
            case stress_type::normal_3D:
                return 3;
            case stress_type::shear_3D:
                return 3;
            default:
                gsWarn << "Invalid stress type\n.";
                return 1;
        };
    }

    /** @brief Each column of the input matrix (u) corresponds to one evaluation point.
     *         Each column of the output matrix (result) corresponds to a set of stress components:
     *         all_2D:    s_11 s_22 s_12
     *         normal_3D: s_11 s_22 s_33
     *         shear_3D:  s_12 s_13 s_23
     *         or is one number in case if stress_type::von_mises is chosen.
     */
    virtual void eval_into(const gsMatrix<T> & u, gsMatrix<T> & result) const;

protected:

    const gsMultiPatch<T>& m_geometry;
    const gsMultiPatch<T>& m_displacement;
    index_t m_patch;
    stress_type::type m_type;
    T m_lambda;
    T m_mu;

}; // class definition ends

/** @brief Compute jacobian determinant of the geometry mapping.
 *         Can be pushed into gsPiecewiseFunction to construct gsField for visualization in Paraview.
*/
template <class T>
class gsDetFunction : public gsFunction<T>
{
public:

    gsDetFunction(const gsMultiPatch<T> & geo, index_t patch)
        : m_geo(geo),
          m_patch(patch)
    {}

    virtual short_t domainDim() const { return m_geo.domainDim(); }

    virtual short_t targetDim() const { return 1; }

    /** @brief Each column of the input matrix (u) corresponds to one evaluation point.
     *         Each column of the output matrix is the jacobian determinant of the mapping at this point.
     */
    virtual void eval_into(const gsMatrix<T> & u, gsMatrix<T> & result) const;

protected:

    gsMultiPatch<T> const & m_geo;
    index_t m_patch;

}; // class definition ends

template <class T>
class gsStiffFunction : public gsFunction<T>
{
public:

    gsStiffFunction(const gsMultiPatch<T> & geo,const gsMultiPatch<T> & disp, index_t patch)
        : m_geo(geo),
          m_disp(disp),
          m_patch(patch)
    {}

    virtual short_t domainDim() const { return m_geo.domainDim(); }

    virtual short_t targetDim() const { return 1; }

    /** @brief Each column of the input matrix (u) corresponds to one evaluation point.
     *         Each column of the output matrix is the jacobian determinant of the mapping at this point.
     */
    virtual void eval_into(const gsMatrix<T> & u, gsMatrix<T> & result) const;

protected:

    gsMultiPatch<T> const & m_geo;
    gsMultiPatch<T> const & m_disp;
    index_t m_patch;

}; // class definition ends

/** @brief Loading function to transfer fluid action to the solid.
 * Used in Fluid-Structure Interaction simulation.
 * Different parametrizations can be used for the geometry+ALE and velocity+pressure
*/
template <class T>
class gsFsiLoad : public gsFunction<T>
{
public:

    gsFsiLoad(const gsMultiPatch<T> & geoRef, const gsMultiPatch<T> & ALEdisplacement,
              index_t patchGeo, boxSide sideGeo,
              const gsMultiPatch<T> & velocity, const gsMultiPatch<T> & pressure,
              index_t patchVelPres, T viscosity, T density)
        : m_geo(geoRef),
          m_ale(ALEdisplacement),
          m_patchGeo(patchGeo),
          m_sideGeo(sideGeo),
          m_vel(velocity),
          m_pres(pressure),
          m_patchVP(patchVelPres),
          m_viscosity(viscosity),
          m_density(density)
    {}

    virtual short_t domainDim() const { return m_geo.domainDim(); }

    virtual short_t targetDim() const { return m_geo.domainDim(); }

    /** @brief Each column of the input matrix (u) corresponds to one evaluation point.
     *         Each column of the output matrix is the jacobian determinant of the mapping at this point.
     */
    virtual void eval_into(const gsMatrix<T> & u, gsMatrix<T> & result) const;

protected:

    gsMultiPatch<T> const & m_geo;
    gsMultiPatch<T> const & m_ale;
    index_t m_patchGeo;
    boxSide m_sideGeo;
    gsMultiPatch<T> const & m_vel;
    gsMultiPatch<T> const & m_pres;
    index_t m_patchVP;
    T m_viscosity;
    T m_density;

}; // class definition ends

} // namespace ends


#ifndef GISMO_BUILD_LIB
#include GISMO_HPP_HEADER(gsElasticityFunctions.hpp)
#endif
