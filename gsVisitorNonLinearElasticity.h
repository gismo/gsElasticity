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

#include <gsElasticity/gsVisitorLinearElasticity.h>

namespace gismo
{

template <class T>
class gsVisitorNonLinearElasticity : public gsVisitorLinearElasticity<T>
{
public:

    typedef gsVisitorLinearElasticity<T> Base;

    gsVisitorNonLinearElasticity(const gsPde<T> & pde_, const gsMultiPatch<T> & deformed)
        : Base(pde_),
          deformedConfiguration(deformed) { }

    void initialize(const gsBasis<T> & basis,
                    const index_t patchIndex,
                    const gsOptionList & options,
                    gsQuadRule<T> & rule)
    {
        Base::initialize(basis,patchIndex,options,rule);

        patch = patchIndex;
        md.flags = NEED_VALUE | NEED_MEASURE | NEED_GRAD_TRANSFORM | NEED_DERIV;
        mdDeformed.flags = NEED_DERIV;
        materialLaw = options.getInt("MaterialLaw");
    }

    inline void evaluate(const gsBasis<T> & basis,
                         const gsGeometry<T> & geo,
                         const gsMatrix<T> & quNodes)
    {
        Base::evaluate(basis,geo,quNodes);

        mdDeformed.points = quNodes;
        deformedConfiguration.patch(patch).computeMap(mdDeformed);
    }

    inline void assemble(gsDomainIterator<T> & element,
                         const gsVector<T> & quWeights)
    {     
        // values of basis functions, 1 x numActiveFunctions
        gsMatrix<T> & basisVals  = basisValues[0];
        // gradients of basis functions, dim x numActiveFunctions
        gsMatrix<T> & basisGrads = basisValues[1];

        index_t dimTensor = dim*(dim+1)/2;
        index_t N = numActiveFunctions;

        // elasticity tensor
        gsMatrix<T> C;
        if (materialLaw == 0)
            Base::setC(C,gsMatrix<T>::Identity(dim,dim),lambda,2*mu);

        // loop over quadrature nodes
        for (index_t q = 0; q < quWeights.rows(); ++q)
        {
            // Multiply weight by the geometry measure
            const T weight = quWeights[q] * md.measure(q);

            // Compute physical gradients of basis functions at q as a dim x numActiveFunction matrix
            gsMatrix<T> physGrad;
            transformGradients(md, q, basisGrads, physGrad);
            // physical jacobian of displacemnt du/dx = du/dxi * dxi/dx
            gsMatrix<T> physDispJac = mdDeformed.jacobian(q)*(md.jacobian(q).cramerInverse());
            // deformation gradient F = I + du/dx
            gsMatrix<T> F = gsMatrix<T>::Identity(dim,dim) + physDispJac;
            // deformation jacobian J = det(F)
            T J = F.determinant();
            // Right Cauchy Green strain, C = F'*F
            gsMatrix<T> RCG = F.transpose() * F;
            // Green-Lagrange strain, E = 0.5*(C-I), a.k.a. full geometric strain tensor
            gsMatrix<T> E = 0.5 * (RCG - gsMatrix<T>::Identity(dim,dim));
            // Green-Lagrange strain as vector
            gsVector<T> Evec(dimTensor);
            // Second Piola-Kirchhoff stress tensor
            gsMatrix<T> S;
            if (materialLaw == 0) // Saint Venant-Kirchhoff
                S = lambda*E.trace()*gsMatrix<T>::Identity(dim,dim) + 2*mu*E;
            if (materialLaw == 1) // neo-Hooke ln(J)
            {
                gsMatrix<T> RCGinv = RCG.cramerInverse();
                S = (lambda*log(J)-mu)*RCGinv + mu*gsMatrix<T>::Identity(dim,dim);
                Base::setC(C,RCGinv,lambda,2*(mu-lambda*log(J)));
            }
            if (materialLaw == 2) // neo-Hooke J^2
            {
                gsMatrix<T> RCGinv = RCG.cramerInverse();
                S = (lambda/2*(J*J-1)-mu)*RCGinv + mu*gsMatrix<T>::Identity(dim,dim);
                Base::setC(C,RCGinv,lambda*J*J,2*(mu-lambda*(J*J-1)));
            }
            // loop over active basis functions (u_i)
            for (index_t i = 0; i < N; i++)
            {
                // Material tangent K_tg_mat = B_i^T * C * B_j;
                gsMatrix<T> B_i;
                Base::setB(B_i,F,physGrad.col(i));
                gsMatrix<T> materialTangentTemp = B_i.transpose() * C;
                // Geometric tangent K_tg_geo = gradB_i^T * S * gradB_j;
                gsVector<T> geometricTangentTemp = S * physGrad.col(i);
                // loop over active basis functions (v_j)
                for (index_t j = 0; j < N; j++)
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
                            localMat(di*N+i, dj*N+j) += weight * materialTangent(di,dj);
                }
                // Second Piola-Kirchhoff stress tensor as vector
                gsVector<T> Svec;
                Base::voigtStress(Svec,S);
                // rhs = -r = force - B*Svec,
                gsVector<T> localResidual = B_i.transpose() * Svec;
                for (short_t d = 0; d < dim; d++)
                    localRhs(d*N+i) -= weight * localResidual(d);
            }
            // contribution of volumetric load function to residual/rhs
            for (short_t d = 0; d < dim; ++d)
                localRhs.middleRows(d*N,N).noalias() += weight * rho * forceValues(d,q) * timeFactor * basisVals.col(q) ;
        }
    }

protected:
    // general problem info
    using Base::dim;
    using Base::pde_ptr;
    const gsMultiPatch<T> & deformedConfiguration;
    index_t patch;
    // geometry mappings
    using Base::md;
    gsMapData<T> mdDeformed;

    // Lame coefficients, material law, density and time factor
    using Base::lambda;
    using Base::mu;
    using Base::rho;
    using Base::timeFactor;
    index_t materialLaw; // (0: St. Venant-Kirchhoff, 1: Neo-Hooke)

    // local components of the global linear system
    using Base::localMat;
    using Base::localRhs;

    // Stored evaluations
    using Base::numActiveFunctions;
    using Base::localIndices;
    using Base::basisValues;
    using Base::forceValues;
};

} // namespace gismo
