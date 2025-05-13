/** @file gsVisitorNonLinearSmallStrainElasticityMM.h

    @brief Element visitor for nonlinear elasticity for 2D plain strain and 3D continua.

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s):
        O. Weeger    (2012 - 2015, TU Kaiserslautern),
        A.Shamanskiy (2016 - ...., TU Kaiserslautern)
*/

#pragma once

#include <gsElasticity/gsVisitorElUtils.h>
#include <gsElasticity/gsBasePde.h>
#include <gsElasticity/gsMaterialBase.h>
#include <gsElasticity/gsMaterialContainer.h>

#include <gsAssembler/gsQuadrature.h>
#include <gsCore/gsFuncData.h>

namespace gismo
{

template <class T>
class gsVisitorLinearSmallStrainElasticityMM : public gsVisitorLinearElasticityMM<T>
{
    using Base = gsVisitorLinearElasticityMM<T>;
public:
    gsVisitorLinearSmallStrainElasticityMM( const gsPde<T> & pde_,
                                            const gsMaterialContainer<T> & materials,
                                            const gsMultiPatch<T> & displacement_,
                                                  gsSparseMatrix<T> * elimMatrix = nullptr)
    :
    Base(pde_, materials, elimMatrix),
    displacement(displacement_)
    {
    }

    inline void evaluate(const gsBasisRefs<T> & basisRefs,
                         const gsGeometry<T> & geo,
                         const gsMatrix<T> & quNodes) override
    {
        GISMO_ASSERT((index_t)geo.id()==patch,"Geometry id ("<<geo.id()<<") does not match with the considered patch index ("<<patch<<")");
        // store quadrature points of the element for geometry evaluation
        md.points = quNodes;
        // NEED_VALUE to get points in the physical domain for evaluation of the RHS
        // NEED_MEASURE to get the Jacobian determinant values for integration
        // NEED_GRAD_TRANSFORM to get the Jacobian matrix to transform gradient from the parametric to physical domain
        md.flags = NEED_VALUE | NEED_MEASURE | NEED_GRAD_TRANSFORM;
        // Compute image of the quadrature points plus gradient, jacobian and other necessary data
        geo.computeMap(md);
        // find local indices of the displacement basis functions active on the element
        basisRefs.front().active_into(quNodes.col(0),localIndicesDisp);
        N_D = localIndicesDisp.rows();
        // Evaluate displacement basis functions and their derivatives on the element
        basisRefs.front().evalAllDers_into(quNodes,1,basisValuesDisp);
        // Evaluate right-hand side at the image of the quadrature points
        pde_ptr->rhs()->eval_into(md.values[0],forceValues);
        // store quadrature points of the element for displacement evaluation
        mdDisplacement.points = quNodes;
        // NEED_DERIV to compute deformation gradient
        mdDisplacement.flags = NEED_DERIV;
        // evaluate displacement gradient
        displacement.patch(patch).computeMap(mdDisplacement);

        mdDeformed = mdDisplacement;
        mdDeformed.values[1] += md.values[1];

        gsMaterialData<T> data;
        gsMaterialBase<T> * material = m_materials.piece(geo.id());
        material->precompute(md,mdDeformed,data,true);
        material->eval_matrix_into(data,matValues);
    }

protected:
    // problem info
    // using Base::dim;
    using Base::patch;
    using Base::pde_ptr;
    using Base::md;
    // local indices (at the current patch) of the displacement basis functions active at the current element
    using Base::localIndicesDisp;
    // number of displacement basis functions active at the current element
    using Base::N_D;
    // values and derivatives of displacement basis functions at quadrature points at the current element
    // values are stored as a N_D x numQuadPoints matrix; not sure about derivatives, must be smth like N_D*dim x numQuadPoints
    using Base::basisValuesDisp;
    // RHS values at quadrature points at the current element; stored as a dim x numQuadPoints matrix
    using Base::forceValues;
    /// Material handler
    using Base::m_materials;
    // current displacement field
    const gsMultiPatch<T> & displacement;
    // evaluation data of the current displacement field
    gsMapData<T> mdDisplacement;
    gsMapData<T> mdDeformed;

    using Base::matValues;
};

} // namespace gismo
