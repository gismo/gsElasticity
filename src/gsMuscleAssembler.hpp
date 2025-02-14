/** @file gsMuscleAssembler.hpp

    @brief Implements gsMuscleAssembler.h

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s):
        A.Shamanskiy (2016 - ...., TU Kaiserslautern)
*/

#pragma once

#include <gsElasticity/gsMuscleAssembler.h>

#include <gsElasticity/gsGeoUtils.h>

// Element visitors
#include <gsElasticity/gsVisitorMuscle.h>
#include <gsElasticity/gsVisitorElasticityNeumann.h>

namespace gismo
{

template<class T>
gsMuscleAssembler<T>::gsMuscleAssembler(gsMultiPatch<T> const & patches,
                                        gsMultiBasis<T> const & basisDisp,
                                        gsMultiBasis<T> const & basisPres,
                                        gsBoundaryConditions<T> const & bconditions,
                                        gsFunction<T> const & body_force,
                                        gsPiecewiseFunction<T> const & muscleTendonDistribution,
                                        const gsVector<T> & fiberDirection)
    : gsElasticityAssembler<T>(patches,basisDisp,basisPres,bconditions,body_force),
      muscleTendon(muscleTendonDistribution),
      fiberDir(fiberDirection)
{

    m_options.addReal("MuscleYoungsModulus","Youngs modulus of the muscle tissue",3.0e5);
    m_options.addReal("TendonYoungsModulus","Youngs modulus of the tendon tissue",3.0e6);
    m_options.addReal("MusclePoissonsRatio","Poisson's ratio of the muscle tissue",0.5);
    m_options.addReal("TendonPoissonsRatio","Poisson's ratio of the tendon tissue",0.5);
    m_options.addReal("MaxMuscleStress","Maximum stress produced at the optimal fiber stretch",3.0e5);
    m_options.addReal("OptFiberStretch","Optimal fiber stretch",1.3);
    m_options.addReal("DeltaW","Shape parameter of the active reponse function",0.3);
    m_options.addReal("PowerNu","Shape parameter of the active reponse function",4.0);
    m_options.addReal("Alpha","Activation parameter of the active muscle response",0.);

    m_options.setInt("MaterialLaw",material_law::muscle);
}

//--------------------- SYSTEM ASSEMBLY ----------------------------------//


template<class T>
bool gsMuscleAssembler<T>::assemble(const gsMatrix<T> & solutionVector,
                                    const std::vector<gsMatrix<T> > & fixedDoFs)
{
    gsMultiPatch<T> displacement,pressure;
    Base::constructSolution(solutionVector,fixedDoFs,displacement,pressure);
    if (m_options.getSwitch("Check"))
        if (checkDisplacement(m_pde_ptr->patches(),displacement) != -1)
            return false;

    m_system.matrix().setZero();
    Base::reserve();
    m_system.rhs().setZero();

    // Compute volumetric integrals and write to the global linear systemz
    gsVisitorMuscle<T> visitor(*m_pde_ptr,muscleTendon,fiberDir,displacement,pressure);
    Base::template push<gsVisitorMuscle<T> >(visitor);
    // Compute surface integrals and write to the global rhs vector
    // change to reuse rhs from linear system
    Base::template push<gsVisitorElasticityNeumann<T> >(m_pde_ptr->bc().neumannSides());

    m_system.matrix().makeCompressed();

    return true;
}

//--------------------- SPECIALS ----------------------------------//

template<class T>
void gsMuscleAssembler<T>::constructCauchyStresses(const gsMultiPatch<T> & displacement,
                                                   const gsMultiPatch<T> & pressure,
                                                   gsPiecewiseFunction<T> & result,
                                                   stress_components::components comp) const
{
    if (comp == stress_components::all_2D_vector || comp == stress_components::all_2D_matrix)
        GISMO_ENSURE(Base::m_dim == 2, "Invalid stress components for a 2D problem");
    if (comp == stress_components::normal_3D_vector || comp == stress_components::shear_3D_vector ||
        comp == stress_components::all_3D_matrix)
        GISMO_ENSURE(Base::m_dim == 3, "Invalid stress type for a 3D problem");
    result.clear();

    for (size_t p = 0; p < m_pde_ptr->domain().nPatches(); ++p )
        result.addPiecePointer(new gsCauchyStressFunction<T>(p,comp,m_options,
                                                             &(m_pde_ptr->domain()), &displacement, &pressure));
}

}// namespace gismo ends
