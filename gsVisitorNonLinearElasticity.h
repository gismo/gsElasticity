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

#include <gsElasticity/gsVisitorBaseElasticity.h>

namespace gismo
{

template <class T>
class gsVisitorNonLinearElasticity : public gsVisitorBaseElasticity<T>
{
public:

    typedef gsVisitorBaseElasticity<T> Base;

    gsVisitorNonLinearElasticity(const gsPde<T> & pde_, const gsMultiPatch<T> & displacement_, bool assembleMatrix_)
        : Base(pde_),
          assembleMatrix(assembleMatrix_),
          displacement(displacement_) { }

    void initialize(const gsBasisRefs<T> & basisRefs,
                    const index_t patchIndex,
                    const gsOptionList & options,
                    gsQuadRule<T> & rule)
    {
        Base::initialize(basisRefs,patchIndex,options,rule);

        patch = patchIndex;
        // NEED_DERIV to compute deformation gradient
        mdDisplacement.flags = NEED_DERIV;

        materialLaw = options.getInt("MaterialLaw");
        T E = options.getReal("YoungsModulus");
        T pr = options.getReal("PoissonsRatio");
        lambda = E * pr / ( ( 1. + pr ) * ( 1. - 2. * pr ) );
        mu     = E / ( 2. * ( 1. + pr ) );
        forceScaling = options.getReal("ForceScaling");
    }

    inline void evaluate(const gsBasisRefs<T> & basisRefs,
                         const gsGeometry<T> & geo,
                         const gsMatrix<T> & quNodes)
    {
        Base::evaluate(basisRefs,geo,quNodes);

        mdDisplacement.points = quNodes;
        displacement.patch(patch).computeMap(mdDisplacement);
    }

    inline void assemble(gsDomainIterator<T> & element,
                         const gsVector<T> & quWeights)
    {     
        // initialize local matrix and rhs
        if (assembleMatrix)
            localMat.setZero(dim*N_D,dim*N_D);
        localRhs.setZero(dim*N_D,1);
        // elasticity tensor
        gsMatrix<T> C;
        if (materialLaw == 0 && assembleMatrix)
            Base::setC(C,gsMatrix<T>::Identity(dim,dim),lambda,mu);
        // loop over quadrature nodes
        for (index_t q = 0; q < quWeights.rows(); ++q)
        {
            // Multiply weight by the geometry measure
            const T weight = quWeights[q] * md.measure(q);
            // Compute physical gradients of basis functions at q as a dim x numActiveFunction matrix
            gsMatrix<T> physGrad;
            transformGradients(md,q,basisValuesDisp[1],physGrad);
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
            // Second Piola-Kirchhoff stress tensor
            gsMatrix<T> S;
            if (materialLaw == 0) // Saint Venant-Kirchhoff
                S = lambda*E.trace()*gsMatrix<T>::Identity(dim,dim) + 2*mu*E;
            if (materialLaw == 1) // neo-Hooke ln(J)
            {
                GISMO_ENSURE(J>0,"Invalid configuration: J < 0");
                gsMatrix<T> RCGinv = RCG.cramerInverse();
                S = (lambda*log(J)-mu)*RCGinv + mu*gsMatrix<T>::Identity(dim,dim);
                if (assembleMatrix)
                    Base::setC(C,RCGinv,lambda,mu-lambda*log(J));
            }
            if (materialLaw == 2) // neo-Hooke J^2
            {
                gsMatrix<T> RCGinv = RCG.cramerInverse();
                S = (lambda/2*(J*J-1)-mu)*RCGinv + mu*gsMatrix<T>::Identity(dim,dim);
                if (assembleMatrix)
                    Base::setC(C,RCGinv,lambda*J*J,mu-lambda*(J*J-1));
            }
            // loop over active basis functions (u_i)
            for (index_t i = 0; i < N_D; i++)
            {
                // Material tangent K_tg_mat = B_i^T * C * B_j;
                gsMatrix<T> B_i;
                Base::setB(B_i,F,physGrad.col(i));
                if (assembleMatrix)
                {
                    gsMatrix<T> materialTangentTemp = B_i.transpose() * C;
                    // Geometric tangent K_tg_geo = gradB_i^T * S * gradB_j;
                    gsVector<T> geometricTangentTemp = S * physGrad.col(i);
                    // loop over active basis functions (v_j)
                    for (index_t j = 0; j < N_D; j++)
                    {
                        gsMatrix<T> B_j;
                        Base::setB(B_j,F,physGrad.col(j));

                        gsMatrix<T> materialTangent = materialTangentTemp * B_j;
                        T geometricTangent =  geometricTangentTemp.transpose() * physGrad.col(j);
                        // K_tg = K_tg_mat + I*K_tg_geo;
                        for (short_t d = 0; d < dim; ++d)
                            materialTangent(d,d) += geometricTangent;

                        for (short_t di = 0; di < dim; ++di)
                            for (short_t dj = 0; dj < dim; ++dj)
                                localMat(di*N_D+i, dj*N_D+j) += weight * materialTangent(di,dj);
                    }
                }
                // Second Piola-Kirchhoff stress tensor as vector
                gsVector<T> Svec;
                Base::voigtStress(Svec,S);
                // rhs = -r = force - B*Svec,
                gsVector<T> localResidual = B_i.transpose() * Svec;
                for (short_t d = 0; d < dim; d++)
                    localRhs(d*N_D+i) -= weight * localResidual(d);
            }
            // contribution of volumetric load function to residual/rhs
            for (short_t d = 0; d < dim; ++d)
                localRhs.middleRows(d*N_D,N_D).noalias() += weight * forceScaling * forceValues(d,q) * basisValuesDisp[0].col(q);
        }
    }

    inline void localToGlobal(const int patchIndex,
                              const std::vector<gsMatrix<T> > & eliminatedDofs,
                              gsSparseSystem<T>     & system)
    {
        std::vector< gsMatrix<unsigned> > globalIndices(dim,localIndicesDisp);
        gsVector<size_t> blockNumbers(dim);
        for (short_t d = 0; d < dim; ++d)
        {
           system.mapColIndices(localIndicesDisp, patchIndex, globalIndices[d], d);
           blockNumbers.at(d) = d;
        }
        system.pushToRhs(localRhs,globalIndices,blockNumbers);
        if (assembleMatrix)
            system.pushToMatrix(localMat,globalIndices,eliminatedDofs,blockNumbers,blockNumbers);
    }

protected:
    //------ inherited ------//
    using Base::dim;
    using Base::pde_ptr;
    using Base::md;
    using Base::localMat;
    using Base::localRhs;
    using Base::N_D;
    using Base::localIndicesDisp;
    using Base::basisValuesDisp;
    using Base::forceValues;

    //------ class specific ----//
    index_t patch; // current patch
    const bool assembleMatrix; // true: assemble matrix and rhs; false: assemble only rhs

    // material law, Lame coefficients, density and force scaling factor
    T lambda, mu, forceScaling;
    index_t materialLaw; // (0: St. Venant-Kirchhoff, 1: ln Neo-Hooke, 2: quadratic Neo-Hooke)

    // current displacement field
    const gsMultiPatch<T> & displacement;
    // evaluation data of the current displacement field
    gsMapData<T> mdDisplacement;
};

} // namespace gismo
