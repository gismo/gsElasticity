/** @file gsVisitorNonLinearElasticity.h

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

#include <gsAssembler/gsQuadrature.h>
#include <gsCore/gsFuncData.h>

namespace gismo
{

template <class T>
class gsVisitorNonLinearElasticity
{
public:
    gsVisitorNonLinearElasticity(const gsPde<T> & pde_, const gsMultiPatch<T> & displacement_)
        : pde_ptr(static_cast<const gsPoissonPde<T>*>(&pde_)),
          displacement(displacement_) { }

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
        materialLaw = options.getInt("MaterialLaw");
        T E = options.getReal("YoungsModulus");
        T pr = options.getReal("PoissonsRatio");
        lambda = E * pr / ( ( 1. + pr ) * ( 1. - 2. * pr ) );
        mu     = E / ( 2. * ( 1. + pr ) );
        forceScaling = options.getReal("ForceScaling");
        localStiffening = options.getReal("LocalStiff");
        // elasticity tensor
        I = gsMatrix<T>::Identity(dim,dim);
        if (materialLaw == 0)
        {
            matrixTraceTensor<T>(C,I,I);
            C *= lambda;
            symmetricIdentityTensor<T>(Ctemp,I);
            C += mu*Ctemp;
        }
        // resize containers for global indices
        globalIndices.resize(dim);
        blockNumbers.resize(dim);
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
            // physical jacobian of displacemnt du/dx = du/dxi * dxi/dx
            physDispJac = mdDisplacement.jacobian(q)*(md.jacobian(q).cramerInverse());
            // deformation gradient F = I + du/dx
            F = I + physDispJac;
            // deformation jacobian J = det(F)
            T J = F.determinant();
            // Right Cauchy Green strain, C = F'*F
            RCG = F.transpose() * F;
            // Green-Lagrange strain, E = 0.5*(C-I), a.k.a. full geometric strain tensor
            E = 0.5 * (RCG - I);
            const T weightBody = quWeights[q] * pow(md.measure(q),-1.*localStiffening) * md.measure(q);
            // Second Piola-Kirchhoff stress tensor
            if (materialLaw == 0) // Saint Venant-Kirchhoff
                S = lambda*E.trace()*I + 2*mu*E;
            if (materialLaw == 1) // neo-Hooke ln(J)
            {
                GISMO_ENSURE(J>0,"Invalid configuration: J < 0");
                RCGinv = RCG.cramerInverse();
                S = (lambda*log(J)-mu)*RCGinv + mu*I;
                // elasticity tensor
                matrixTraceTensor<T>(C,RCGinv,RCGinv);
                C *= lambda;
                symmetricIdentityTensor<T>(Ctemp,RCGinv);
                C += (mu-lambda*log(J))*Ctemp;
            }
            if (materialLaw == 2) // quad neo-Hooke
            {
                RCGinv = RCG.cramerInverse();
                S = (lambda*(J*J-1)/2-mu)*RCGinv + mu*I;
                // elasticity tensor
                matrixTraceTensor<T>(C,RCGinv,RCGinv);
                C *= lambda*J*J;
                symmetricIdentityTensor<T>(Ctemp,RCGinv);
                C += (mu-lambda*(J*J-1)/2)*Ctemp;
            }
            // loop over active basis functions (u_i)
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
    const gsPoissonPde<T> * pde_ptr;
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
    // current displacement field
    const gsMultiPatch<T> & displacement;
    // evaluation data of the current displacement field
    gsMapData<T> mdDisplacement;

    // all temporary matrices defined here for efficiency
    gsMatrix<T> C, Ctemp, physGrad, physDispJac, F, RCG, E, S, RCGinv, B_i, materialTangentTemp, B_j, materialTangent, I;
    gsVector<T> geometricTangentTemp, Svec, localResidual;
    T localStiffening;
    // containers for global indices
    std::vector< gsMatrix<unsigned> > globalIndices;
    gsVector<size_t> blockNumbers;
};

} // namespace gismo
