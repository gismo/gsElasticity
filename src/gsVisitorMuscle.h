/** @file gsVisitorMuscle.h

    @brief Visitor class for the nonlinear elasticity solver with a muscle material.
    The material model is based on the paper by M.H.Gfrerer and B.Simeon
    "Fiber-based modeling and simulation of skeletal muscles"

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s):
        A.Shamanskiy (2016 - ...., TU Kaiserslautern)
*/

#pragma once

#include <gsElasticity/src/gsVisitorElUtils.h>
#include <gsElasticity/src/gsBasePde.h>

#include <gsAssembler/gsQuadrature.h>
#include <gsCore/gsFuncData.h>

namespace gismo
{

template <class T>
class gsVisitorMuscle
{
public:
    gsVisitorMuscle(const gsPde<T> & pde_,
                    const gsPiecewiseFunction<T> & muscleTendon_,
                    const gsVector<T> & fiberDir_,
                    const gsMultiPatch<T> & displacement_,
                    const gsMultiPatch<T> & pressure_)
        : pde_ptr(static_cast<const gsBasePde<T>*>(&pde_)),
          muscleTendon(muscleTendon_),
          fiberDir(fiberDir_),
          displacement(displacement_),
          pressure(pressure_){}

    void initialize(const gsBasisRefs<T> & basisRefs,
                    const index_t patchIndex,
                    const gsOptionList & options,
                    gsQuadRule<T> & rule)
    {
        // parametric dimension of the first displacement component
        dim = basisRefs.front().dim();
        // a quadrature rule is defined by the basis for the first velocity component.
        // the same rule is used for the presure
        rule = gsQuadrature::get(basisRefs.front(), options);
        // saving necessary info
        patch = patchIndex;
        T YM_muscle = options.getReal("MuscleYoungsModulus");
        T PR_muscle = options.getReal("MusclePoissonsRatio");
        T YM_tendon = options.getReal("TendonYoungsModulus");
        T PR_tendon = options.getReal("TendonPoissonsRatio");
        lambda_inv_muscle = ( 1. + PR_muscle ) * ( 1. - 2. * PR_muscle ) / YM_muscle / PR_muscle ;
        mu_muscle     = YM_muscle / ( 2. * ( 1. + PR_muscle ) );
        lambda_inv_tendon = ( 1. + PR_tendon) * ( 1. - 2. * PR_tendon) / YM_tendon / PR_tendon;
        mu_tendon     = YM_tendon / ( 2. * ( 1. + PR_tendon ) );
        forceScaling = options.getReal("ForceScaling");
        maxMuscleStress = options.getReal("MaxMuscleStress");
        optFiberStretch = options.getReal("OptFiberStretch");
        deltaW = options.getReal("DeltaW");
        powerNu = options.getReal("PowerNu");
        alpha = options.getReal("Alpha"); // activation parameter
        I = gsMatrix<T>::Identity(dim,dim);
        // resize containers for global indices
        globalIndices.resize(dim+1);
        blockNumbers.resize(dim+1);
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
        basisRefs.front().active_into(quNodes.col(0),localIndicesDisp);
        N_D = localIndicesDisp.rows();
        basisRefs.back().active_into(quNodes.col(0), localIndicesPres);
        N_P = localIndicesPres.rows();
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
        // evaluate muscle-tendon distribution
        muscleTendon.piece(patch).eval_into(quNodes,muscleTendonValues);
    }

    inline void assemble(gsDomainIteratorWrapper<T> & element,
                         const gsVector<T> & quWeights)
    {
        GISMO_UNUSED(element);
        // Initialize local matrix/rhs                      // A | B^T
        localMat.setZero(dim*N_D + N_P, dim*N_D + N_P);     // --|--    matrix structure
        localRhs.setZero(dim*N_D + N_P,1);                  // B | C
        // Loop over the quadrature nodes
        for (index_t q = 0; q < quWeights.rows(); ++q)
        {
            // Compute material parameters
            const T mu = muscleTendonValues.at(q) * mu_muscle + (1-muscleTendonValues.at(q))*mu_tendon;
            const T lambda_inv = muscleTendonValues.at(q) * lambda_inv_muscle + (1-muscleTendonValues.at(q))*lambda_inv_tendon;
            // Multiply quadrature weight by the geometry measure
            const T weight = quWeights[q] * md.measure(q);
            // Compute physical gradients of basis functions at q as a dim x numActiveFunction matrix
            transformGradients(md,q,basisValuesDisp[1],physGradDisp);
            // physical jacobian of displacemnt du/dx = du/dxi * dxi/dx
            physDispJac = mdDisplacement.jacobian(q)*(md.jacobian(q).cramerInverse());
            // deformation gradient F = I + du/dx
            F = I + physDispJac;
            // deformation jacobian J = det(F)
            T J = F.determinant();
            // Right Cauchy Green strain, C = F'*F
            RCG = F.transpose() * F;
            // logarithmic neo-Hooke
            GISMO_ENSURE(J>0,"Invalid configuration: J < 0");
            RCGinv = RCG.cramerInverse();
            // Second Piola-Kirchhoff stress tensor, passive part
            S = (pressureValues.at(q)-mu)*RCGinv + mu*I;
            /// active stress contribution - start
            // fiber direction in the physical domain
            fiberDirPhys = md.jacobian(q)*fiberDir;
            fiberDirPhys /= fiberDirPhys.norm();
            // dyadic product of the fiber direction
            M = fiberDirPhys * fiberDirPhys.transpose();
            // active stress scaled with the time activation parameter
            T fiberStretch = sqrt((M*RCG).trace());
            T ratioInExp = (fiberStretch/optFiberStretch-1)/deltaW;
            T megaExp = exp(-1*pow(abs(ratioInExp),powerNu));
            S += M * maxMuscleStress * alpha * muscleTendonValues.at(q)/ pow(fiberStretch,2) * megaExp;
            /// active stress contribution - end
            // elasticity tensor
            symmetricIdentityTensor<T>(C,RCGinv);
            C *= mu-pressureValues.at(q);
            /// active stress contribution - start

            matrixTraceTensor<T>(Ctemp,M,M);
            C += -1*Ctemp*alpha*maxMuscleStress*megaExp/pow(fiberStretch,3)* muscleTendonValues.at(q)*
                    (2 + powerNu*pow(ratioInExp,powerNu-1)/deltaW/optFiberStretch);
            /// active stress contribution - end
            // Matrix A and reisdual: loop over displacement basis functions
            for (index_t i = 0; i < N_D; i++)
            {
                setB<T>(B_i,F,physGradDisp.col(i));
                materialTangentTemp = B_i.transpose() * C;
                // Geometric tangent K_tg_geo = gradB_i^T * S * gradB_j;
                geometricTangentTemp = S * physGradDisp.col(i);
                // A-matrix
                for (index_t j = 0; j < N_D; j++)
                {
                    setB<T>(B_j,F,physGradDisp.col(j));
                    materialTangent = materialTangentTemp * B_j;
                    T geometricTangent =  geometricTangentTemp.transpose() * physGradDisp.col(j);
                    // K_tg = K_tg_mat + I*K_tg_geo;
                    for (short_t d = 0; d < dim; ++d)
                        materialTangent(d,d) += geometricTangent;

                    for (short_t di = 0; di < dim; ++di)
                        for (short_t dj = 0; dj < dim; ++dj)
                            localMat(di*N_D+i, dj*N_D+j) += weight * materialTangent(di,dj);
                }

                // Second Piola-Kirchhoff stress tensor as vector
                voigtStress<T>(Svec,S);
                // rhs = -r = force - B*Svec,
                localResidual = B_i.transpose() * Svec;
                for (short_t d = 0; d < dim; d++)
                    localRhs(d*N_D+i) -= weight * localResidual(d);
            }
            // B-matrix
            divV = F.cramerInverse().transpose() * physGradDisp;
            for (short_t d = 0; d < dim; ++d)
            {
                block = weight*basisValuesPres.col(q)*divV.row(d);
                localMat.block(dim*N_D,d*N_D,N_P,N_D) += block.block(0,0,N_P,N_D);
                localMat.block(d*N_D,dim*N_D,N_D,N_P) += block.transpose().block(0,0,N_D,N_P);
            }
            // C-matrix
            if (abs(lambda_inv) > 0)
                localMat.block(dim*N_D,dim*N_D,N_P,N_P) -=
                        (weight*lambda_inv*basisValuesPres.col(q)*basisValuesPres.col(q).transpose()).block(0,0,N_P,N_P);
            // rhs: constraint residual
            localRhs.middleRows(dim*N_D,N_P) += weight*basisValuesPres.col(q)*(lambda_inv*pressureValues.at(q)-log(J));
            // rhs: force
            for (short_t d = 0; d < dim; ++d)
                localRhs.middleRows(d*N_D,N_D).noalias() += weight * forceScaling * forceValues(d,q) * basisValuesDisp[0].col(q) ;
        }
    }

    inline void localToGlobal(const int patchIndex,
                              const std::vector<gsMatrix<T> > & eliminatedDofs,
                              gsSparseSystem<T> & system)
    {
        // computes global indices for displacement components
        for (short_t d = 0; d < dim; ++d)
        {
            system.mapColIndices(localIndicesDisp,patchIndex,globalIndices[d],d);
            blockNumbers.at(d) = d;
        }
        // computes global indices for pressure
        system.mapColIndices(localIndicesPres, patchIndex, globalIndices[dim], dim);
        blockNumbers.at(dim) = dim;
        // push to global system
        system.pushToRhs(localRhs,globalIndices,blockNumbers);
        system.pushToMatrix(localMat,globalIndices,eliminatedDofs,blockNumbers,blockNumbers);
    }

protected:
    // problem info
    short_t dim;
    const gsBasePde<T> * pde_ptr;
    // distribution of the tendon and muscle tissue; given in the parametric domain
    const gsPiecewiseFunction<T> & muscleTendon;
    // orientation of muscle fibers in the parametric domain
    const gsVector<T> & fiberDir;
    index_t patch; // current patch
    // Lame coefficients and force scaling factor
    T lambda_inv_muscle, mu_muscle, lambda_inv_tendon, mu_tendon, forceScaling;
    // Active response parameters
    T maxMuscleStress, optFiberStretch, deltaW, powerNu, alpha;
    // geometry mapping
    gsMapData<T> md;
    // local components of the global linear system
    gsMatrix<T> localMat;
    gsMatrix<T> localRhs;
    // local indices (at the current patch) of basis functions active at the current element
    gsMatrix<index_t> localIndicesDisp;
    gsMatrix<index_t> localIndicesPres;
    // number of displacement and pressure basis functions active at the current element
    index_t N_D, N_P;
    // values and derivatives of displacement basis functions at quadrature points at the current element
    // values are stored as a N_D x numQuadPoints matrix; not sure about derivatives, must be smth like N_D*dim x numQuadPoints
    std::vector<gsMatrix<T> > basisValuesDisp;
    // values of pressure basis functions active at the current element;
    // stores as a N_P x numQuadPoints matrix
    gsMatrix<T> basisValuesPres;
    // RHS values at quadrature points at the current element; stored as a dim x numQuadPoints matrix
    gsMatrix<T> forceValues;
    // current displacement field
    const gsMultiPatch<T> & displacement;
    // evaluation data of the current displacement field
    gsMapData<T> mdDisplacement;
    // current pressure field
    const gsMultiPatch<T> & pressure;
    // evaluation data of the current pressure field stored as a 1 x numQuadPoints matrix
    gsMatrix<T> pressureValues;
    // evaluation data of the muscle-tendon distribution stored as a 1 x numQuadPoints matrix
    gsMatrix<T> muscleTendonValues;

    // all temporary matrices defined here for efficiency
    gsMatrix<T> C, Ctemp, physGradDisp, physDispJac, F, RCG, E, S, RCGinv, B_i, materialTangentTemp, B_j, materialTangent, divV, block, I, M;
    gsVector<T> geometricTangentTemp, Svec, localResidual, fiberDirPhys;
    // containers for global indices
    std::vector< gsMatrix<index_t> > globalIndices;
    gsVector<index_t> blockNumbers;
};

} // namespace gismo
