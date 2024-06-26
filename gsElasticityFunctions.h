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

#include <gsCore/gsMultiPatch.h>
#include <gsElasticity/gsBaseUtils.h>
#include <gsIO/gsOptionList.h>

namespace gismo
{

/** @brief Compute Cauchy stresses for a previously computed/defined displacement field.
 *         Can be pushed into gsPiecewiseFunction to construct gsField for visualization in Paraview.
*/
template <class T>
class gsCauchyStressFunction : public gsFunction<T>
{
public:

    gsCauchyStressFunction(index_t patch, stress_components::components comp,
                           const gsOptionList & options,
                           const gsMultiPatch<T> * geometry,
                           const gsMultiPatch<T> * displacement,
                           const gsMultiPatch<T> * pressure = nullptr,
                           const gsMultiPatch<T> * velocity = nullptr)
        : m_geometry(geometry),
          m_displacement(displacement),
          m_pressure(pressure),
          m_velocity(velocity),
          m_patch(patch),
          m_dim(m_geometry->patch(m_patch).parDim()),
          m_options(options),
          m_type(comp)
    {}



    virtual short_t domainDim() const
    {
        return m_geometry->patch(m_patch).parDim();
    }

    virtual short_t targetDim() const
    {
        switch(m_type)
        {
        case stress_components::von_mises: return 1;
        case stress_components::all_2D_vector: return 3;
        case stress_components::all_2D_matrix: return 2;
        case stress_components::normal_3D_vector: return 3;
        case stress_components::shear_3D_vector: return 3;
        case stress_components::all_3D_matrix: return 3;
        default: return 0;
        };
    }

    /** @brief Each column of the input matrix (u) corresponds to one evaluation point.
     *         Columns of the output matrix (result) correspond to a set of stress components for vector cases,
     *         or consist of one number in case of stress_type::von_mises,
     *         or form square matrices concatinated in the col-direction.
     */
    virtual void eval_into(const gsMatrix<T> & u, gsMatrix<T> & result) const
    {
        switch (material_law::law(m_options.getInt("MaterialLaw")))
        {
        case material_law::hooke : linearElastic(u,result); return;
        case material_law::saint_venant_kirchhoff : nonLinearElastic(u,result); return;
        case material_law::neo_hooke_ln : nonLinearElastic(u,result); return;
        case material_law::neo_hooke_quad : nonLinearElastic(u,result); return;
        case material_law::mixed_hooke : mixedLinearElastic(u,result); return;
        case material_law::mixed_neo_hooke_ln : mixedNonLinearElastic(u,result); return;
        case material_law::muscle : mixedNonLinearElastic(u,result); return;
        default: return;
        }
    }

protected:

    /// size of the output matrix according to the m_type
    index_t outputCols(index_t inputCols) const
    {
        switch (m_type)
        {
        case stress_components::von_mises: return inputCols;
        case stress_components::all_2D_vector: return inputCols;
        case stress_components::all_2D_matrix: return 2*inputCols;
        case stress_components::normal_3D_vector: return inputCols;
        case stress_components::shear_3D_vector: return inputCols;
        case stress_components::all_3D_matrix: return 3*inputCols;
        default: return 0;
        }
    }
    /// save components of the stress tensor to the output matrix according to the m_type
    void saveStress(const gsMatrix<T> & S, gsMatrix<T> & result, index_t q) const;

    /// computation routines for different material laws
    void linearElastic(const gsMatrix<T> & u, gsMatrix<T> & result) const;
    void nonLinearElastic(const gsMatrix<T> & u, gsMatrix<T> & result) const;
    void mixedLinearElastic(const gsMatrix<T> & u, gsMatrix<T> & result) const;
    void mixedNonLinearElastic(const gsMatrix<T> & u, gsMatrix<T> & result) const;

protected:
    const gsMultiPatch<T> * m_geometry;
    const gsMultiPatch<T> * m_displacement;
    const gsMultiPatch<T> * m_pressure;
    const gsMultiPatch<T> * m_velocity;
    index_t m_patch;
    short_t m_dim;
    const gsOptionList & m_options;
    stress_components::components m_type;

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

    virtual short_t targetDim() const { return m_geo.targetDim(); }

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

/** @brief Loading function to transfer fluid action to the solid boundary.
 * Used in Fluid-Structure Interaction simulation.
 * Different parametrizations can be used for the geometry+ALE and velocity+pressure
*/
template <class T>
class gsFsiBoundaryLoad : public gsFunction<T>
{
public:

    /** @brief  
     * 
     * param[in] geoBoundary    boundary of the geometry matching the ALEboundary
     * param[in] ALEboundary    boundary of the ALE mesh matching the geoBoundary
     * param[in] velocity       patch of the velocity field adjacent to the interface
     * param[in] pressure       patch of the pressure field adjacent to the interface
     * param[in] viscosity      viscosity of the fluid
     * param[in] density        density of the fluid
     * 
    */
    gsFsiBoundaryLoad(const gsMultiPatch<T> & geo, patchSide geoSide,
                      const gsMultiPatch<T> & ALE, patchSide ALESide,
                      const gsMultiPatch<T> & velocity, patchSide velSide,
                      const gsMultiPatch<T> & pressure, patchSide presSide,
                      T viscosity, 
                      T density,
                      bool parametric = false)
        : m_geo(geo),
          m_geoSide(geoSide),
          m_ale(ALE),
          m_aleSide(ALESide),
          m_vel(velocity),
          m_velSide(velSide),
          m_pres(pressure),
          m_presSide(presSide),
          m_viscosity(viscosity),
          m_density(density),
          m_parametric(parametric)
    {}

    virtual short_t domainDim() const { return m_geo.patch(m_geoSide.patch).boundary(m_geoSide.side())->domainDim(); }

    virtual short_t targetDim() const { return m_geo.targetDim(); }

    /** @brief Each column of the input matrix (u) corresponds to one evaluation point.
     *         Each column of the output matrix is the jacobian determinant of the mapping at this point.
     */
    virtual void eval_into(const gsMatrix<T> & u, gsMatrix<T> & result) const;

protected:

    gsMultiPatch<T> const & m_geo;
    patchSide m_geoSide;
    gsMultiPatch<T> const & m_ale;
    patchSide m_aleSide;
    gsMultiPatch<T> const & m_vel;
    patchSide m_velSide;
    gsMultiPatch<T> const & m_pres;
    patchSide m_presSide;
    T m_viscosity;
    T m_density;
    bool m_parametric;

}; // class definition ends

} // namespace ends


#ifndef GISMO_BUILD_LIB
#include GISMO_HPP_HEADER(gsElasticityFunctions.hpp)
#endif
