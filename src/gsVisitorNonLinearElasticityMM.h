/** @file gsVisitorNonLinearElasticityMM.h

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

// template <class T>
// class gsAsDeformed : public gsFunction<T>
// {
// public:
//     gsAsDeformed(const gsFunction<T> & undeformed_,
//                  const gsFunction<T> & deformed_)
//     :
//     undeformed(undeformed_),
//     deformed(deformed_)
//     {
//     }

//     short_t targetDim() const override { return undeformed.targetDim(); }
//     short_t domainDim() const override { return undeformed.domainDim(); }

//     void eval_into(const gsMatrix<T> & u, gsMatrix<T> & result) const override
//     {
//         undeformed.eval_into(u,result);
//         result += deformed.eval(u);
//     }

//     void deriv_into(const gsMatrix<T> & u, gsMatrix<T> & result) const override
//     {
//         undeformed.deriv_into(u,result);
//         result += deformed.deriv(u);
//     }

//     void deriv2_into(const gsMatrix<T> & u, gsMatrix<T> & result) const override
//     {
//         undeformed.deriv2_into(u,result);
//         result += deformed.deriv2(u);
//     }

//     void evalAllDers_into(const gsMatrix<T> & u, const int n, std::vector<gsMatrix<T> > & result, bool sameElement) const
//     {
//         std::vector<gsMatrix<T>> defResult;
//         undeformed.evalAllDers_into(u,n,result,sameElement);
//         deformed.evalAllDers_into(u,n,defResult,sameElement);
//         for (index_t i = 0; i <= n; ++i)
//             result[i] += defResult[i];
//     }

//     std::ostream &print(std::ostream &os) const
//     {
//         gsInfo<<"As deformed function\n";
//         undeformed.print(os);
//         deformed.print(os);
//         return os;
//     }

//     GISMO_CLONE_FUNCTION(gsAsDeformed)

// protected:
//     const gsFunction<T> & undeformed;
//     const gsFunction<T> & deformed;
// };



template <class T>
class gsVisitorNonLinearElasticityMM
{
public:
    gsVisitorNonLinearElasticityMM( const gsPde<T> & pde_,
                                    const gsMaterialContainer<T> & materials,
                                    const gsMultiPatch<T> & displacement_)
    :
    pde_ptr(static_cast<const gsBasePde<T>*>(&pde_)),
    m_materials(materials),
    displacement(displacement_)
    {
    }

    void initialize(const gsBasisRefs<T> & basisRefs,
                    const index_t patchIndex,
                    const gsOptionList & options,
                    gsQuadRule<T> & rule)
    {
        // parametric dimension of the first displacement component
        dim = basisRefs.front().dim();
        // a quadrature rule is defined by the basis for the first displacement component.
        rule = gsQuadrature::get(basisRefs.front(), options);
        // saving necessary info
        patch = patchIndex;
        forceScaling = options.getReal("ForceScaling");
        localStiffening = options.getReal("LocalStiff");
        I = gsMatrix<T>::Identity(dim,dim);
        // resize containers for global indices
        globalIndices.resize(dim);
        blockNumbers.resize(dim);
    }

    inline void evaluate(const gsBasisRefs<T> & basisRefs,
                         const gsGeometry<T> & geo,
                         const gsMatrix<T> & quNodes)
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

        // gsAsDeformed<T> def(geo,displacement.patch(patch));

        // // Construct a material evaluator for matrix with id geo.id()
        // gsMaterialEvalSingle<T,gsMaterialOutput::S,false> Seval(0,m_materials.piece(geo.id()),geo,def);
        // gsMaterialEvalSingle<T,gsMaterialOutput::C, true> Ceval(0,m_materials.piece(geo.id()),geo,def);

        // // The evaluator is single-piece, hence we use piece(0)
        // Seval.eval_into(quNodes,vecValues);
        // Ceval.eval_into(quNodes,matValues);

        // The lines below are faster than calling gsMaterialEval(Single),
        // Since we re-use the geometric data and compute the parameter data only once
        gsMaterialData<T> data;
        gsMaterialBase<T> * material = m_materials.piece(geo.id());
        material->precompute(md,mdDeformed,data);
        material->eval_deformation_gradient_into(data,defGradValues); // re-use the pre-computed ones
        material->eval_stress_into(data,vecValues);
        material->eval_matrix_into(data,matValues);
    }

    inline void assemble(gsDomainIterator<T> & element,
                         const gsVector<T> & quWeights)
    {
        // initialize local matrix and rhs
        localMat.setZero(dim*N_D,dim*N_D);
        localRhs.setZero(dim*N_D,1);
        // loop over quadrature nodes
        for (index_t q = 0; q < quWeights.rows(); ++q)
        {
            const T weightForce = quWeights[q] * md.measure(q);
            // Compute physical gradients of basis functions at q as a dim x numActiveFunction matrix
            transformGradients(md,q,basisValuesDisp[1],physGrad);
            // deformation gradient
            F = defGradValues.reshapeCol(q,dim,dim);
            const T weightBody = quWeights[q] * pow(md.measure(q),-1.*localStiffening) * md.measure(q);
            // Elasticity tensor
            C = matValues.reshapeCol(q,math::sqrt(matValues.rows()),math::sqrt(matValues.rows()));
            // Second Piola-Kirchhoff stress tensor
            S = vecValues.reshapeCol(q,math::sqrt(vecValues.rows()),math::sqrt(vecValues.rows()));
            for (index_t i = 0; i < N_D; i++)
            {
                // Material tangent K_tg_mat = B_i^T * C * B_j;
                setB<T>(B_i,F,physGrad.col(i));
                materialTangentTemp = B_i.transpose() * C;
                // Geometric tangent K_tg_geo = gradB_i^T * S * gradB_j;
                geometricTangentTemp = S * physGrad.col(i);
                // loop over active basis functions (v_j)
                for (index_t j = 0; j < N_D; j++)
                {
                    setB<T>(B_j,F,physGrad.col(j));

                    materialTangent = materialTangentTemp * B_j;
                    T geometricTangent =  geometricTangentTemp.transpose() * physGrad.col(j);
                    // K_tg = K_tg_mat + I*K_tg_geo;
                    for (short_t d = 0; d < dim; ++d)
                        materialTangent(d,d) += geometricTangent;

                    for (short_t di = 0; di < dim; ++di)
                        for (short_t dj = 0; dj < dim; ++dj)
                            localMat(di*N_D+i, dj*N_D+j) += weightBody * materialTangent(di,dj);
                }
                // Second Piola-Kirchhoff stress tensor as vector
                voigtStress<T>(Svec,S);
                // rhs = -r = force - B*Svec,
                localResidual = B_i.transpose() * Svec;
                for (short_t d = 0; d < dim; d++)
                    localRhs(d*N_D+i) -= weightBody * localResidual(d);
            }
            // contribution of volumetric load function to residual/rhs
            for (short_t d = 0; d < dim; ++d)
                localRhs.middleRows(d*N_D,N_D).noalias() += weightForce * forceScaling * forceValues(d,q) * basisValuesDisp[0].col(q);
        }
    }

    inline void localToGlobal(const int patchIndex,
                              const std::vector<gsMatrix<T> > & eliminatedDofs,
                              gsSparseSystem<T> & system)
    {
        // computes global indices for displacement components
        for (short_t d = 0; d < dim; ++d)
        {
            system.mapColIndices(localIndicesDisp, patchIndex, globalIndices[d], d);
            blockNumbers.at(d) = d;
        }
        // push to global system
        system.pushToRhs(localRhs,globalIndices,blockNumbers);
        system.pushToMatrix(localMat,globalIndices,eliminatedDofs,blockNumbers,blockNumbers);
    }

protected:
    // problem info
    short_t dim;
    index_t patch; // current patch
    const gsBasePde<T> * pde_ptr;
    index_t materialLaw; // (0: St. Venant-Kirchhoff, 1: ln neo-Hooke, 2: quad neo-Hooke)
    // Lame coefficients and force scaling factor
    T lambda, mu, forceScaling;
    // geometry mapping
    gsMapData<T> md;
    // local components of the global linear system
    gsMatrix<T> localMat;
    gsMatrix<T> localRhs;
    // local indices (at the current patch) of the displacement basis functions active at the current element
    gsMatrix<index_t> localIndicesDisp;
    // number of displacement basis functions active at the current element
    index_t N_D;
    // values and derivatives of displacement basis functions at quadrature points at the current element
    // values are stored as a N_D x numQuadPoints matrix; not sure about derivatives, must be smth like N_D*dim x numQuadPoints
    std::vector<gsMatrix<T> > basisValuesDisp;
    // RHS values at quadrature points at the current element; stored as a dim x numQuadPoints matrix
    gsMatrix<T> forceValues;
    /// Material handler
    const gsMaterialContainer<T> & m_materials;
    // current displacement field
    const gsMultiPatch<T> & displacement;
    // evaluation data of the current displacement field
    gsMapData<T> mdDisplacement;
    gsMapData<T> mdDeformed;

    // all temporary matrices defined here for efficiency
    gsMatrix<T> C, Ctemp, physGrad, physDispJac, F, RCG, E, S, RCGinv, B_i, materialTangentTemp, B_j, materialTangent, I;
    gsVector<T> geometricTangentTemp, Svec, localResidual;
    T localStiffening;
    // containers for global indices
    std::vector< gsMatrix<index_t> > globalIndices;
    gsVector<index_t> blockNumbers;

    gsMatrix<T> matValues;
    gsMatrix<T> vecValues;
    gsMatrix<T> defGradValues;

};

} // namespace gismo
