/** @file gsVisitorMixedNonLinearElasticity.h

    @brief Visitor class for volumetric integration of the mixed nonlinear elasticity system.

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s):
        O. Weeger    (2012 - 2015, TU Kaiserslautern),
        A.Shamanskiy (2016 - ...., TU Kaiserslautern)
*/

#pragma once

#include <gsElasticity/gsVisitorBaseElasticity.h>

namespace gismo
{

template <class T>
class gsVisitorMixedNonLinearElasticity : public gsVisitorBaseElasticity<T>
{
public:
    typedef gsVisitorBaseElasticity<T> Base;

    gsVisitorMixedNonLinearElasticity(const gsPde<T> & pde_, const gsMultiPatch<T> & displacement_,
                                      const gsMultiPatch<T> & pressure_, bool assembleMatrix_)
        : Base(pde_,assembleMatrix_),
          displacement(displacement_),
          pressure(pressure_) {}

    void initialize(const gsBasisRefs<T> & basisRefs,
                    const index_t patchIndex,
                    const gsOptionList & options,
                    gsQuadRule<T> & rule)
    {
        Base::initialize(basisRefs,patchIndex,options,rule);

        patch = patchIndex;

        T E = options.getReal("YoungsModulus");
        T pr = options.getReal("PoissonsRatio");
        lambda_inv = ( 1. + pr ) * ( 1. - 2. * pr ) / E / pr ;
        mu     = E / ( 2. * ( 1. + pr ) );
        forceScaling = options.getReal("ForceScaling");
    }

    inline void evaluate(const gsBasisRefs<T> & basisRefs,
                         const gsGeometry<T> & geo,
                         const gsMatrix<T> & quNodes)
    {
        // store quadrature points of the element for geometry evaluation
        md.points = quNodes;
        // NEED_VALUE to get points in the physical domain for evaluation of the RHS
        // NEED_MEASURE to get the Jacobian determinant values for integration
        // NEED_GRAD_TRANSFORM to get the Jacobian matrix to transform gradient from the parametric to physical domain
        md.flags = NEED_VALUE | NEED_MEASURE | NEED_GRAD_TRANSFORM;
        // Compute image of the quadrature points plus gradient, jacobian and other necessary data
        geo.computeMap(md);
        // find local indices of the displacement and pressure basis functions active on the element
        localIndices.resize(2);
        basisRefs.front().active_into(quNodes.col(0),localIndices[0]);
        N_D = localIndices[0].rows();
        basisRefs.back().active_into(quNodes.col(0), localIndices[1]);
        N_P = localIndices[1].rows();
        // Evaluate displacement basis functions and their derivatives on the element
        basisRefs.front().evalAllDers_into(quNodes,1,basisValuesDisp);
        // Evaluate pressure basis functions on the element
        basisRefs.back().eval_into(quNodes,basisValuesPres);
        // Evaluate right-hand side at the image of the quadrature points
        pde_ptr->rhs()->eval_into(md.values[0],forceValues);

        // store quadrature points of the element for displacement evaluation
        mdDisplacement.points = quNodes;
        // NEED_DERIV to compute deformation gradient
        mdDisplacement.flags = NEED_DERIV;
        // evaluate displacement gradient
        displacement.patch(patch).computeMap(mdDisplacement);

        // evaluate pressure; we use eval_into instead of another gsMapData object
        // because it easier for simple value evaluation
        pressure.patch(patch).eval_into(quNodes,pressureValues);
    }

    inline void assemble(gsDomainIterator<T>    & element,
                         gsVector<T> const      & quWeights)
    {
        // Initialize local matrix/rhs
        if (assembleMatrix)                                     // A | B^T
            localMat.setZero(dim*N_D + N_P, dim*N_D + N_P);     // --|--    matrix structure
        localRhs.setZero(dim*N_D + N_P,1);                      // B | C
        // Loop over the quadrature nodes
        for (index_t q = 0; q < quWeights.rows(); ++q)
        {
            // Multiply quadrature weight by the geometry measure
            const T weight = quWeights[q] * md.measure(q);
            // Compute physical gradients of basis functions at q as a dim x numActiveFunction matrix
            gsMatrix<T> physGradDisp;
            transformGradients(md,q,basisValuesDisp[1],physGradDisp);
            // physical jacobian of displacemnt du/dx = du/dxi * dxi/dx
            gsMatrix<T> physDispJac = mdDisplacement.jacobian(q)*(md.jacobian(q).cramerInverse());
            // deformation gradient F = I + du/dx
            gsMatrix<T> F = gsMatrix<T>::Identity(dim,dim) + physDispJac;
            // deformation jacobian J = det(F)
            T J = F.determinant();
            // Right Cauchy Green strain, C = F'*F
            gsMatrix<T> RCG = F.transpose() * F;
            // Green-Lagrange strain, E = 0.5*(C-I), a.k.a. full geometric strain tensor
            gsMatrix<T> E = 0.5 * (RCG - gsMatrix<T>::Identity(dim,dim));
            // logarithmic neo-Hooke
            GISMO_ENSURE(J>0,"Invalid configuration: J < 0");
            gsMatrix<T> RCGinv = RCG.cramerInverse();
            // Second Piola-Kirchhoff stress tensor
            gsMatrix<T> S = (pressureValues.at(q)-mu)*RCGinv + mu*gsMatrix<T>::Identity(dim,dim);
            gsMatrix<T> C; // elasticity tensor
            if (assembleMatrix)
                Base::setC(C,RCGinv,0.,mu-pressureValues.at(q));

            // Loop over displacement basis functions
            for (index_t i = 0; i < N_D; i++)
            {
                gsMatrix<T> B_i;
                Base::setB(B_i,F,physGradDisp.col(i));
                if (assembleMatrix)
                {
                    gsMatrix<T> materialTangentTemp = B_i.transpose() * C;
                    // Geometric tangent K_tg_geo = gradB_i^T * S * gradB_j;
                    gsVector<T> geometricTangentTemp = S * physGradDisp.col(i);
                    // Loop for A-matrix
                    for (index_t j = 0; j < N_D; j++)
                    {
                        gsMatrix<T> B_j;
                        Base::setB(B_j,F,physGradDisp.col(j));
                        gsMatrix<T> materialTangent = materialTangentTemp * B_j;
                        T geometricTangent =  geometricTangentTemp.transpose() * physGradDisp.col(j);
                        // K_tg = K_tg_mat + I*K_tg_geo;
                        for (short_t d = 0; d < dim; ++d)
                            materialTangent(d,d) += geometricTangent;

                        for (short_t di = 0; di < dim; ++di)
                            for (short_t dj = 0; dj < dim; ++dj)
                                localMat(di*N_D+i, dj*N_D+j) += weight * materialTangent(di,dj);
                    }
                    // Loop for B-matrix
                    for (index_t j = 0; j < N_P; j++)
                    {
                        gsMatrix<T> divV = physGradDisp.col(i).transpose()*F.cramerInverse();
                        for (short_t d = 0; d < dim; ++d)
                        {
                            localMat(dim*N_D+j,d*N_D+i) += weight*divV.at(d)*basisValuesPres(j,q);
                            localMat(d*N_D+i,dim*N_D+j) += weight*divV.at(d)*basisValuesPres(j,q);
                        }

                    }
                    // Loop over pressure basis functions for C-matrix
                    if (abs(lambda_inv) > 0)
                        for (index_t i = 0; i < N_P; ++i)
                            for (index_t j = 0; j < N_P; ++j)
                                localMat(dim*N_D+i,dim*N_D+j) -= weight*lambda_inv*basisValuesPres(i,q)*basisValuesPres(j,q);

                }
                // Second Piola-Kirchhoff stress tensor as vector
                gsVector<T> Svec;
                Base::voigtStress(Svec,S);
                // rhs = -r = force - B*Svec,
                gsVector<T> localResidual = B_i.transpose() * Svec;
                for (short_t d = 0; d < dim; d++)
                    localRhs(d*N_D+i) -= weight * localResidual(d);

            }
            // rhs contribution
            for (short_t d = 0; d < dim; ++d)
                localRhs.middleRows(d*N_D,N_D).noalias() += weight * forceScaling * forceValues(d,q) * basisValuesDisp[0].col(q) ;
        }
    }

protected:
    //------ inherited ------//
    using Base::dim;
    using Base::pde_ptr;
    using Base::assembleMatrix;
    using Base::md;
    using Base::localMat;
    using Base::localRhs;
    using Base::localIndices;
    using Base::N_D;
    using Base::basisValuesDisp;
    using Base::forceValues;

    //------ class specific ----//
    index_t patch; // current patch

    // Lame coefficients, density and force scaling factor
    T lambda_inv, mu, forceScaling;

    // current displacement field
    const gsMultiPatch<T> & displacement;
    // evaluation data of the current displacement field
    gsMapData<T> mdDisplacement;

    // current pressure field
    const gsMultiPatch<T> & pressure;
    // evaluation data of the current pressure field stored as a 1 x numQuadPoints matrix
    gsMatrix<T> pressureValues;
    // number of pressure basis functions active at the current element
    index_t N_P;
    // values of pressure basis functions active at the current element;
    // stores as a N_P x numQuadPoints matrix
    gsMatrix<T> basisValuesPres;
};

} // namespace gismo
