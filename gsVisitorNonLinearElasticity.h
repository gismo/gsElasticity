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
        if (materialLaw == 0)
            assembleStVK(element, quWeights);
        else
            assembleNeoH(element, quWeights);
    }

    inline void assembleNeoH(gsDomainIterator<T> & element,
                              const gsVector<T> & quWeights)
    {
        /*
        gsMatrix<T> & bVals  = basisData[0];
        gsMatrix<T> & bGrads = basisData[1];
        //const gsMatrix<T> & defGrads = m_deformation->jacobians();

        for (index_t k = 0; k < quWeights.rows(); ++k) // loop over quadrature nodes
        {
            // Multiply weight by the geometry measure
            weight = quWeights[k] * geoEval.measure(k);

            // compute physical gradients at k as a Dim x NumActive matrix
            geoEval.transformGradients(k, bGrads, physGrad);

            //geoEval.transformGradients(k,
            //gsAsConstMatrix<T,3,3> geoJac(geoEval->jacobian(k).data(),3,3 );
            //gsAsConstMatrix<T,3,3> defDer(m_deformation->jacobian(k).data(),3,3 );
            defDer_k = m_deformation->jacobian(k).transpose();
            //defDer_kv = defDer_k.asVector();
            //geoEval.transformGradients(k, defDer_kv, displGrad);

            // Displacement gradient, H = du = (dx/dxi)'^-1 * (du/dxi)'
            displGrad = geoEval.jacobian(k).transpose().inverse() * defDer_k;

            // Deformation gradient, F = I + H
            defGrad = displGrad.transpose();
            for (size_t di = 0; di < m_dim; di++)
                defGrad(di,di) += 1.;

            // Determinant of deformation gradient, J = det(F)
            detF = defGrad.determinant();
                        logdetF = math::log(detF);
            mulamlogJ = m_mu - m_lambda * logdetF;

            // Inverse of Fi = F^-1
            defGrad_inv = defGrad.inverse();

            // Local internal force vector contribution, mu*(F-Fi') + lambda*log(J)*Fi'
            locResMat = (weight * m_mu) * (defGrad - defGrad_inv.transpose()) + (weight * m_lambda * logdetF) * defGrad_inv.transpose();

            // 1st basis function (U/i)
            for (index_t i = 0; i < numActive; i++)
            {
                locResVec = locResMat * physGrad.col(i);

                // Spatial dimensions of 1st basis function
                for (size_t di = 0; di < m_dim; di++)
                {
                    // Write to Rhs
                    localRhs(di*numActive+i) -= locResVec(di);

                    // Write gradient as matrix
                    gradU.setZero();
                    gradU.row(di) = physGrad.col(i);

                    defGrad_inv_gradU = defGrad_inv * gradU;				// Fi*dU
                    defGrad_inv_gradU_trace = defGrad_inv_gradU.trace();	// tr(Fi*dU) = Fi':dU

                    // 2nd basis function (V/j)
                    //for (index_t j = 0; j < numActive; j++)
                    // Exploit symmetry of K
                    for (index_t j = i; j < numActive; j++)
                    {
                        // Spatial dimensions of 2nd basis function
                        for (size_t dj = 0; dj < m_dim; dj++)
                        {
                            // Write gradient as matrix
                            gradV.setZero();
                            gradV.row(dj) = physGrad.col(j);

                            defGrad_inv_gradV = defGrad_inv * gradV;		// Fi*dV

                            // Local tangent stiffnees matrix contribution
                            locKtgVal = m_mu * ( gradU.transpose()*gradV ).trace()
                                      + mulamlogJ * ( defGrad_inv_gradU * defGrad_inv_gradV ).trace()
                                      + m_lambda * defGrad_inv_gradU_trace * defGrad_inv_gradV.trace();

                            // Write to Mat
                            localMat(di*numActive+i, dj*numActive+j) += weight * locKtgVal;
                        }
                    }
                }
            }

            // Local external force vector contribution
            for (size_t j = 0; j < m_dim; ++j)
                localRhs.middleRows(j*numActive,numActive).noalias() +=
                    weight * m_rho * forceVals(j,k) * m_tfac * bVals.col(k) ;
        }
        //gsDebug<< "local Mat: \n"<< localMat << "\n";*/
    }

    inline void assembleStVK(gsDomainIterator<T> & element,
                              const gsVector<T> & quWeights)
    {        
        // values of basis functions, 1 x numActiveFunctions
        gsMatrix<T> & basisVals  = basisValues[0];
        // gradients of basis functions, dim x numActiveFunctions
        gsMatrix<T> & basisGrads = basisValues[1];

        index_t dimTensor = dim*(dim+1)/2;
        index_t N = numActiveFunctions;

        // elasticity tensor
        gsMatrix<T> C(dimTensor,dimTensor);
        C.setZero();
        if (dim == 2)
        {
            C(0,0) = C(1,1) = 2*mu + lambda;
            C(0,1) = C(1,0) = lambda;
            C(2,2) = mu;
        }
        else if (dim == 3)
        {
            C(0,0) = C(1,1) = C(2,2) = 2*mu + lambda;
            C(0,1) = C(1,0) = C(1,2) = C(2,1) = C(2,0) = C(0,2)= lambda;
            C(3,3) = C(4,4) = C(5,5) = mu;
        }

        // loop over quadrature nodes
        for (index_t q = 0; q < quWeights.rows(); ++q)
        {
            // Multiply weight by the geometry measure
            const T weight = quWeights[q] * md.measure(q);

            // Compute physical gradients of basis functions at q as a dim x numActiveFunction matrix
            gsMatrix<T> physGrad;
            transformGradients(md, q, basisGrads, physGrad);
            // parametric jacobian of displacement du/dxi
            gsMatrix<T> parDispJac = mdDeformed.jacobian(q);
            // physical jacobian of displacemnt du/dx = du/dxi * dxi/dx
            gsMatrix<T> physDispJac = parDispJac*(md.jacobian(q).cramerInverse());
            // deformation gradient F = I + du/dx
            gsMatrix<T> F = gsMatrix<T>::Identity(dim,dim) + physDispJac;
            // Right Cauchy Green strain, C = F'*F
            gsMatrix<T> RCG = F.transpose() * F;
            // Green-Lagrange strain, E = 0.5*(C-I), a.k.a. full geometric strain tensor
            gsMatrix<T> E = 0.5 * (RCG - gsMatrix<T>::Identity(dim,dim));
            // Green-Lagrange strain as vector
            gsVector<T> Evec(dimTensor);
            // Second Piola-Kirchhoff stress tensor
            gsMatrix<T> S(dim,dim);
            // Second Piola-Kirchhoff stress tensor as vector
            gsVector<T> Svec;
            if (dim == 2)
            {
                Evec(0) = E(0,0);
                Evec(1) = E(1,1);
                Evec(2) = 2.*E(0,1);

                Svec = C * Evec;

                S(0,0) = Svec(0);
                S(1,1) = Svec(1);
                S(0,1) = S(1,0) = Svec(2);
            }
            else if (dim == 3)
            {
                Evec(0) = E(0,0);
                Evec(1) = E(1,1);
                Evec(2) = E(2,2);
                Evec(3) = 2.*E(0,1);
                Evec(4) = 2.*E(1,2);
                Evec(5) = 2.*E(0,2);

                Svec = C * Evec;

                S(0,0) = Svec(0);
                S(1,1) = Svec(1);
                S(2,2) = Svec(2);
                S(0,1) = S(1,0) = Svec(3);
                S(1,2) = S(2,1) = Svec(4);
                S(0,2) = S(2,0) = Svec(5);
            }

            // 1st basis function (U/i)
            for (index_t i = 0; i < N; i++)
            {
                // auxiliary matrix B, see doctor thesis of O. Weeger, p.37
                gsMatrix<T> B_i;
                setMatB(B_i,F,physGrad.col(i));
                // Geometric tangent K_tg_geo = gradB_i^T * S * gradB_j;
                gsVector<T> geometricTangentTemp = S * physGrad.col(i);
                // Material tangent K_tg_mat = B_i^T * C * B_j;
                gsMatrix<T> materialTangentTemp = B_i.transpose() * C;
                // index_t j = i - symmetry
                for (index_t j = 0; j < N; j++)
                {
                    gsMatrix<T> B_j;
                    setMatB(B_j,F,physGrad.col(j));

                    gsMatrix<T> materialTangent = materialTangentTemp * B_j;
                    T geometricTangent =  geometricTangentTemp.transpose() * physGrad.col(j);
                    // K_tg = K_tg_mat + I*K_tg_geo;
                    for (index_t d = 0; d < dim; ++d)
                        materialTangent(d,d) += geometricTangent;

                    for (index_t di = 0; di < dim; ++di)
                        for (index_t dj = 0; dj < dim; ++dj)
                            localMat(di*N+i, dj*N+j) += weight * materialTangent(di,dj);
                }
                // rhs = -r = force - B*PK2,
                gsVector<T> localResidual = B_i.transpose() * Svec;
                for (index_t d = 0; d < dim; d++)
                    localRhs(d*N+i) -= weight * localResidual(d);
            }

            // contribution of volumetric load function to residual/rhs
            for (index_t d = 0; d < dim; ++d)
                localRhs.middleRows(d*N,N).noalias() += weight * rho * forceValues(d,q) * timeFactor * basisVals.col(q) ;
        }
    }

protected:

    inline void setMatB( gsMatrix<T> &B, const gsMatrix<T> &F, const gsVector<T> &bGrad)
    {
        B.resize(dim*(dim+1)/2,dim);

        for (index_t j = 0; j < dim; ++j)
        {
            for (index_t i = 0; i < dim; ++i)
                B(i,j) = F(j,i) * bGrad(i);

            if (dim == 2)
                B(2,j) = F(j,0) * bGrad(1) + F(j,1) * bGrad(0);

            if (dim == 3)
                for (index_t i = 0; i < dim; ++i)
                {
                    index_t k = (i+1)%dim;
                    B(i+dim,j) = F(j,i) * bGrad(k) + F(j,k) * bGrad(i);
                }
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
