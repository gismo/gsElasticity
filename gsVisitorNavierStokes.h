/** @file gsVisitorNavierStokes.h

    @brief Visitor class for volumetric integration of thetangential Navier-Stokes system.

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s):
        A.Shamanskiy (2016 - ...., TU Kaiserslautern)
*/

#pragma once

#include <gsAssembler/gsQuadrature.h>
#include <gsCore/gsFuncData.h>
#include <algorithm>

namespace gismo
{

template <class T>
class gsVisitorNavierStokes
{
public:

    gsVisitorNavierStokes(const gsPde<T> & pde_, const gsMultiPatch<T> & velocity_,
                          const gsMultiPatch<T> & pressure_, bool assembleMatrix_)
        : pde_ptr(static_cast<const gsPoissonPde<T>*>(&pde_)),
          assembleMatrix(assembleMatrix_),
          velocity(velocity_),
          pressure(pressure_),
          ALE(false) {}

    gsVisitorNavierStokes(const gsPde<T> & pde_, const gsMultiPatch<T> & velocity_,
                          const gsMultiPatch<T> & pressure_,
                          const gsMultiPatch<T> & ALEdisplacement, bool assembleMatrix_)
        : pde_ptr(static_cast<const gsPoissonPde<T>*>(&pde_)),
          assembleMatrix(assembleMatrix_),
          velocity(velocity_),
          pressure(pressure_),
          dispALE(ALEdisplacement),
          ALE(true) {}

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
        iterationType = options.getInt("Iteration");
        supg = options.getSwitch("SUPG");
        viscosity = options.getReal("Viscosity");
        patch = patchIndex;
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
        // NEED_2ND_DER to transform hessians to physical domain
        if (supg)
            md.flags = NEED_VALUE | NEED_MEASURE | NEED_GRAD_TRANSFORM | NEED_2ND_DER;
        else
            md.flags = NEED_VALUE | NEED_MEASURE | NEED_GRAD_TRANSFORM;
        // Compute image of the quadrature points plus gradient, jacobian and other necessary data
        geo.computeMap(md);
        // find local indices of the velocity and pressure basis functions active on the element
        basisRefs.front().active_into(quNodes.col(0),localIndicesVel);
        N_V = localIndicesVel.rows();
        basisRefs.back().active_into(quNodes.col(0), localIndicesPres);
        N_P = localIndicesPres.rows();
        // Evaluate velocity basis functions and their derivatives on the element (and hessians, if SUPG is used)
        basisRefs.front().evalAllDers_into(quNodes,(supg ? 2 : 1),basisValuesVel);
        // Evaluate pressure basis functions on the element
        basisRefs.back().eval_into(quNodes,basisValuesPres);
        // Evaluate gradients of pressure basis functions if SUPG is used
        if (supg)
            basisRefs.back().deriv_into(quNodes,basisGradsPres);
        // Evaluate right-hand side at the image of the quadrature points
        pde_ptr->rhs()->eval_into(md.values[0],forceValues);
        // store quadrature points of the element for velocity evaluation
        mdVelocity.points = quNodes;
        // NEED_VALUE to compute velocity values
        // NEED_DERIV to compute velocity gradient
        // NEED_2ND_DER to compute velocity hessian
        if (supg)
            mdVelocity.flags = NEED_VALUE | NEED_DERIV | NEED_2ND_DER;
        else
            mdVelocity.flags = NEED_VALUE | NEED_DERIV;
        // evaluate velocity
        velocity.patch(patch).computeMap(mdVelocity);
        // evaluate pressure values; for some reason, pressure evaluation via gsMapData doesn't work
        pressure.patch(patch).eval_into(quNodes,pressureValues);
        // evaluate pressure grads
        if (supg)
            pressure.patch(patch).deriv_into(quNodes,pressureGrads);
        if (ALE)
        {
            mdALE.flags = NEED_DERIV;
            mdALE.points = quNodes;
            dispALE.patch(patch).computeMap(mdALE);
        }
    }

    inline void assemble(gsDomainIterator<T> & element,
                         const gsVector<T> & quWeights)
    {
        if (iterationType == 0 && !ALE)
            assembleOseen(element,quWeights);
        else if (iterationType == 1 && !ALE)
            assembleNewton(element,quWeights);
        else if (ALE)
            assembleNewtonALE(element,quWeights);
    }

    inline void localToGlobal(const int patchIndex,
                              const std::vector<gsMatrix<T> > & eliminatedDofs,
                              gsSparseSystem<T> & system)
    {
        // number of unknowns: dim of velocity + 1 for pressure
        std::vector< gsMatrix<unsigned> > globalIndices(dim+1);
        gsVector<size_t> blockNumbers(dim+1);
        // computes global indices for velocity components
        for (short_t d = 0; d < dim; ++d)
        {
            system.mapColIndices(localIndicesVel, patchIndex, globalIndices[d], d);
            blockNumbers.at(d) = d;
        }
        // computes global indices for pressure
        system.mapColIndices(localIndicesPres, patchIndex, globalIndices[dim], dim);
        blockNumbers.at(dim) = dim;
        // push to global system
        system.pushToRhs(localRhs,globalIndices,blockNumbers);
        if (assembleMatrix)
            system.pushToMatrix(localMat,globalIndices,eliminatedDofs,blockNumbers,blockNumbers);
    }

protected:
    // estimate the element size
    T cellSize(gsDomainIterator<T> & element)
    {
        gsMatrix<T> lowLeft = element.lowerCorner();
        gsMatrix<T> upperRight = element.upperCorner();
        gsMatrix<T> lowRight = lowLeft;
        lowRight.at(0) = upperRight.at(0);
        gsMatrix<T> upperLeft = lowLeft;
        upperLeft.at(1) = upperRight.at(1);

        T diag1 = (pde_ptr->domain().patch(patch).eval(lowLeft) -
                   pde_ptr->domain().patch(patch).eval(upperRight)).norm();
        T diag2 = (pde_ptr->domain().patch(patch).eval(upperLeft) -
                   pde_ptr->domain().patch(patch).eval(lowRight)).norm();
        return (diag1+diag2)/2/sqrt(2);
    }

    void assembleNewton(gsDomainIterator<T> & element,
                        const gsVector<T> & quWeights)
    {
        // Initialize local matrix/rhs
        if (assembleMatrix)                                     // A | D
            localMat.setZero(dim*N_V + N_P, dim*N_V + N_P);     // --|--    matrix structure
        localRhs.setZero(dim*N_V + N_P,1);                      // B | 0
        // roughly estimate h - diameter of the element ( for SUPG)
        // T h = cellSize(element);
        // Loop over the quadrature nodes
        for (index_t q = 0; q < quWeights.rows(); ++q)
        {
            // Multiply quadrature weight by the geometry measure
            const T weight = quWeights[q] * md.measure(q);
            // Compute physical gradients of the velocity basis functions at q as a dim x numActiveFunction matrix
            gsMatrix<T> physGradVel;
            transformGradients(md, q, basisValuesVel[1], physGradVel);
            // Compute physical Jacobian of the current velocity field
            gsMatrix<T> physJacCurVel = mdVelocity.jacobian(q)*(md.jacobian(q).cramerInverse());
            // matrix A: diffusion
            gsMatrix<T> block = weight*viscosity * physGradVel.transpose()*physGradVel;
            for (short_t d = 0; d < dim; ++d)
                localMat.block(d*N_V,d*N_V,N_V,N_V) += block.block(0,0,N_V,N_V);
            // matrix A: advection
            block = weight*basisValuesVel[0].col(q) * (mdVelocity.values[0].col(q).transpose()*physGradVel);
            for (short_t d = 0; d < dim; ++d)
                localMat.block(d*N_V,d*N_V,N_V,N_V) += block.block(0,0,N_V,N_V);
            // matrix A: reaction
            block = weight*basisValuesVel[0].col(q) * basisValuesVel[0].col(q).transpose();
            for (short_t di = 0; di < dim; ++di)
                for (short_t dj = 0; dj < dim; ++dj)
                    localMat.block(di*N_V,dj*N_V,N_V,N_V) += physJacCurVel(di,dj)*block.block(0,0,N_V,N_V);
            // matrices B and D
            for (short_t d = 0; d < dim; ++d)
            {
                block = weight*basisValuesPres.col(q)*physGradVel.row(d);
                localMat.block(dim*N_V,d*N_V,N_P,N_V) -= block.block(0,0,N_P,N_V); // B
                localMat.block(d*N_V,dim*N_V,N_V,N_P) -= block.transpose().block(0,0,N_V,N_P); // D
            }
            // rhs: force
            for (short_t d = 0; d < dim; ++d)
                localRhs.middleRows(d*N_V,N_V).noalias() += weight * forceScaling *
                    forceValues(d,q) * basisValuesVel[0].col(q);
            // rhs: residual diffusion
            for (short_t d = 0; d < dim; ++d)
                localRhs.middleRows(d*N_V,N_V).noalias() -= weight * viscosity *
                    (physJacCurVel.row(d)*physGradVel).transpose();
            // rhs: residual advection
            for (short_t d = 0; d < dim; ++d)
                localRhs.middleRows(d*N_V,N_V).noalias() -= weight *
                    (physJacCurVel.row(d) * mdVelocity.values[0].col(q))(0,0) * basisValuesVel[0].col(q);
            // rhs: residual pressure
            for (short_t d = 0; d < dim; ++d)
                localRhs.middleRows(d*N_V,N_V).noalias() += weight *
                    pressureValues.at(q) * physGradVel.row(d).transpose();
            // rhs: constraint residual
            localRhs.middleRows(dim*N_V,N_P).noalias() += weight *
                basisValuesPres.col(q) * physJacCurVel.trace();
            // SUPG stabilization
            // see D.Kuzmin "A guide to numerical methods for transport equations" pp.67-71
            // and J.Volker " FEM for Incompressible Flow Problems" p. 286
            if (supg)
            {
                /*T velNorm = mdVelocity.values[0].col(q).norm();
                T Pe_h = velNorm * h / viscosity; // Peclet number
                T alpha = std::min(1.,Pe_h/6);
                T tau = abs(velNorm) > 1e-14 ? alpha * h / velNorm / 2. : 0.; // stabilization parameter
                // physical gradients of velocity basis functions multiplied by the current velocity; N_V x 1 matrix
                gsMatrix<T> advectedPhysGradVel = (mdVelocity.values[0].col(q).transpose() * physGradVel).transpose();
                // matrix A: diffusion stabilization
                gsMatrix<T> physLaplVel; // laplacians of the basis functions as a 1 x N_V matrix
                transformLaplaceHgrad(md,q,basisValuesVel[1],basisValuesVel[2],physLaplVel);
                block = weight*viscosity*tau*advectedPhysGradVel*physLaplVel;
                for (short_t d = 0; d < dim; ++d)
                    localMat.block(d*N_V,d*N_V,N_V,N_V) -= block.block(0,0,N_V,N_V);
                // matrix A: advection stabilization
                block = weight*tau* advectedPhysGradVel*advectedPhysGradVel.transpose();
                for (short_t d = 0; d < dim; ++d)
                    localMat.block(d*N_V,d*N_V,N_V,N_V) += block.block(0,0,N_V,N_V);
                // matrix A: reaction stabilization
                block = weight*tau*advectedPhysGradVel* basisValuesVel[0].col(q).transpose();
                for (short_t di = 0; di < dim; ++di)
                    for (short_t dj = 0; dj < dim; ++dj)
                        localMat.block(di*N_V,dj*N_V,N_V,N_V) += physJacCurVel(di,dj)*block.block(0,0,N_V,N_V);
                // matrix D: pressure stabilization
                // Compute physical gradients of the pressure basis functions at q as a dim x numActiveFunction matrix
                gsMatrix<T> physGradPres;
                transformGradients(md, q, basisGradsPres, physGradPres);
                for (short_t d = 0; d < dim; ++d)
                    localMat.block(d*N_V,dim*N_V,N_V,N_P) += (weight*tau*advectedPhysGradVel*
                                                              physGradPres.row(d)).block(0,0,N_V,N_P);
                // rhs: force stabilization
                for (short_t d = 0; d < dim; ++d)
                    localRhs.middleRows(d*N_V,N_V).noalias() += weight * forceScaling * tau *
                        forceValues(d,q) * advectedPhysGradVel;
                // rhs: diffusion stabilization
                gsMatrix<T> physLaplCurVel;
                transformLaplaceHgrad(md,q,mdVelocity.values[1],mdVelocity.values[2],physLaplCurVel);
                for (short_t d = 0; d < dim; ++d)
                    localRhs.middleRows(d*N_V,N_V).noalias() += weight * tau * viscosity *
                        physLaplCurVel.at(d) * advectedPhysGradVel;
                // rhs: nonlinear term stabilization
                for (short_t d = 0; d < dim; ++d)
                    localRhs.middleRows(d*N_V,N_V).noalias() -= weight * tau *
                        (physJacCurVel.row(d) * mdVelocity.values[0].col(q))(0,0) * advectedPhysGradVel;
                // rhs: pressure stabilization
                // Compute physical Jacobian of the current pressure field
                gsMatrix<T> physJacCurPres;
                transformGradients(md, q, pressureGrads, physJacCurPres);
                for (short_t d = 0; d < dim; ++d)
                    localRhs.middleRows(d*N_V,N_V).noalias() -= weight * tau *
                        physJacCurPres.at(d) * advectedPhysGradVel;*/
            }
        }
    }

    void assembleOseen(gsDomainIterator<T> & element,
                       const gsVector<T> & quWeights)
    {
        // Initialize local matrix/rhs
        if (assembleMatrix)                                     // A | D
            localMat.setZero(dim*N_V + N_P, dim*N_V + N_P);     // --|--    matrix structure
        localRhs.setZero(dim*N_V + N_P,1);                      // B | 0
        // roughly estimate h - diameter of the element ( for SUPG)
        //T h = cellSize(element);

        // Loop over the quadrature nodes
        for (index_t q = 0; q < quWeights.rows(); ++q)
        {
            // Multiply quadrature weight by the geometry measure
            const T weight = quWeights[q] * md.measure(q);
            // Compute physical gradients of the velocity basis functions at q as a dim x numActiveFunction matrix
            gsMatrix<T> physGradVel;
            transformGradients(md, q, basisValuesVel[1], physGradVel);
            // Compute physical Jacobian of the current velocity field
            gsMatrix<T> physJacCurVel = mdVelocity.jacobian(q)*(md.jacobian(q).cramerInverse());
            // matrix A: diffusion
            gsMatrix<T> block = weight*viscosity * physGradVel.transpose()*physGradVel;
            for (short_t d = 0; d < dim; ++d)
                localMat.block(d*N_V,d*N_V,N_V,N_V) += block.block(0,0,N_V,N_V);
            // matrix A: advection
            block = weight*basisValuesVel[0].col(q) * (mdVelocity.values[0].col(q).transpose()*physGradVel);
            for (short_t d = 0; d < dim; ++d)
                localMat.block(d*N_V,d*N_V,N_V,N_V) += block.block(0,0,N_V,N_V);
            // matrices B and D
            for (short_t d = 0; d < dim; ++d)
            {
                block = weight*basisValuesPres.col(q)*physGradVel.row(d);
                localMat.block(dim*N_V,d*N_V,N_P,N_V) -= block.block(0,0,N_P,N_V); // B
                localMat.block(d*N_V,dim*N_V,N_V,N_P) -= block.transpose().block(0,0,N_V,N_P); // D
            }
            // rhs: force
            for (short_t d = 0; d < dim; ++d)
                localRhs.middleRows(d*N_V,N_V).noalias() += weight * forceScaling *
                    forceValues(d,q) * basisValuesVel[0].col(q);
            // SUPG stabilization
            // see D.Kuzmin "A guide to numerical methods for transport equations" pp.67-71
            // and J.Volker " FEM for Incompressible Flow Problems" p. 286
            if (supg)
            {
                /*T velNorm = mdVelocity.values[0].col(q).norm();
                T Pe_h = velNorm * h / viscosity; // Peclet number
                T alpha = std::min(1.,Pe_h/6);
                T tau = abs(velNorm) > 1e-14 ? alpha * h / velNorm / 2. : 0.; // stabilization parameter
                // physical gradients of velocity basis functions multiplied by the current velocity; N_V x 1 matrix
                gsMatrix<T> advectedPhysGradVel = (mdVelocity.values[0].col(q).transpose() * physGradVel).transpose();
                // matrix A: diffusion stabilization
                gsMatrix<T> physLaplVel; // laplacians of the basis functions as a 1 x N_V matrix
                transformLaplaceHgrad(md,q,basisValuesVel[1],basisValuesVel[2],physLaplVel);
                block = weight*viscosity*tau*advectedPhysGradVel*physLaplVel;
                for (short_t d = 0; d < dim; ++d)
                    localMat.block(d*N_V,d*N_V,N_V,N_V) -= block.block(0,0,N_V,N_V);
                // matrix A: advection stabilization
                block = weight*tau* advectedPhysGradVel*advectedPhysGradVel.transpose();
                for (short_t d = 0; d < dim; ++d)
                    localMat.block(d*N_V,d*N_V,N_V,N_V) += block.block(0,0,N_V,N_V);
                // matrix D: pressure stabilization
                // Compute physical gradients of the pressure basis functions at q as a dim x numActiveFunction matrix
                gsMatrix<T> physGradPres;
                transformGradients(md, q, basisGradsPres, physGradPres);
                for (short_t d = 0; d < dim; ++d)
                    localMat.block(d*N_V,dim*N_V,N_V,N_P) += (weight*tau*advectedPhysGradVel*
                                                              physGradPres.row(d)).block(0,0,N_V,N_P);
                // rhs: force stabilization
                for (short_t d = 0; d < dim; ++d)
                    localRhs.middleRows(d*N_V,N_V).noalias() += weight * forceScaling * tau *
                        forceValues(d,q) * advectedPhysGradVel;   */
            }
        }
    }

    void assembleNewtonALE(gsDomainIterator<T> & element,
                           const gsVector<T> & quWeights)
    {
        // Initialize local matrix/rhs
        if (assembleMatrix)                                     // A | D
            localMat.setZero(dim*N_V + N_P, dim*N_V + N_P);     // --|--    matrix structure
        localRhs.setZero(dim*N_V + N_P,1);                      // B | 0
        // roughly estimate h - diameter of the element ( for SUPG)
        // T h = cellSize(element);
        // Loop over the quadrature nodes
        for (index_t q = 0; q < quWeights.rows(); ++q)
        {
            // Multiply quadrature weight by the geometry measure
            const T weight = quWeights[q] * md.measure(q);
            // Compute physical gradients of the velocity basis functions at q as a dim x numActiveFunction matrix
            gsMatrix<T> physGradVel;
            transformGradients(md, q, basisValuesVel[1], physGradVel);
            // Compute physical Jacobian of the current velocity field
            gsMatrix<T> physJacCurVel = mdVelocity.jacobian(q)*(md.jacobian(q).cramerInverse());
            // Compute physical Jacobian of the ALE mapping
            gsMatrix<T> physJacALE = gsMatrix<T>::Identity(dim,dim) +
                    mdALE.jacobian(q)*(md.jacobian(q).cramerInverse());
            // ALE mapping determinant
            T J = physJacALE.determinant();
            GISMO_ENSURE(J>0,"Invalid configuration: J < 0");
            // inverse ALE jacobian
            gsMatrix<T> invJacALE = physJacALE.cramerInverse();
            // matrix A: diffusion
            gsMatrix<T> B_i, B_j;
            for (short_t di = 0; di < dim; ++di)
                for (index_t i = 0; i < N_V; ++i)
                {
                    setUFFU(B_i,invJacALE,physGradVel.col(i),di);
                    for (short_t dj = 0; dj < dim; ++dj)
                        for (index_t j = 0; j < N_V; ++j)
                        {
                            setUFFU(B_j,invJacALE,physGradVel.col(j),dj);
                            localMat(di*N_V+i,dj*N_V+j) += 0.5*weight*J*viscosity*
                                    (B_i.array()*B_j.array()).sum();
                        }
                }
            /*gsMatrix<T> block = weight*viscosity * physGradVel.transpose()*physGradVel;
            for (short_t d = 0; d < dim; ++d)
                localMat.block(d*N_V,d*N_V,N_V,N_V) += block.block(0,0,N_V,N_V);*/
            // matrix A: advection
            // transformed advecting velocity
            gsMatrix<T> Fu = invJacALE*mdVelocity.values[0].col(q);
            gsMatrix<T> block = weight*J*basisValuesVel[0].col(q) * (Fu.transpose()*physGradVel);
            for (short_t d = 0; d < dim; ++d)
                localMat.block(d*N_V,d*N_V,N_V,N_V) += block.block(0,0,N_V,N_V);
            // matrix A: reaction
            for (short_t dj = 0; dj < dim; ++dj)
            {
                block = weight*J*physJacCurVel*(invJacALE.col(dj)*basisValuesVel[0].col(q).transpose());
                for (short_t di = 0; di < dim; ++di)
                    localMat.block(di*N_V,dj*N_V,N_V,N_V) += (basisValuesVel[0].col(q)*block.row(di)).block(0,0,N_V,N_V);
            }
            // matrices B and D
            for (short_t d = 0; d < dim; ++d)
            {
                block = weight*J*(basisValuesPres.col(q) * invJacALE.transpose().row(d))*physGradVel;
                localMat.block(dim*N_V,d*N_V,N_P,N_V) -= block.block(0,0,N_P,N_V); // D
                localMat.block(d*N_V,dim*N_V,N_V,N_P) -= block.transpose().block(0,0,N_V,N_P); // B
            }
            // rhs: force
            for (short_t d = 0; d < dim; ++d)
                localRhs.middleRows(d*N_V,N_V).noalias() += weight * forceScaling *
                    J*forceValues(d,q) * basisValuesVel[0].col(q);
            // rhs: residual diffusion
            block = (physJacCurVel*invJacALE + invJacALE.transpose()*physJacCurVel.transpose()) *
                    invJacALE.transpose()*weight*J*viscosity;
            for (short_t d = 0; d < dim; ++d)
                localRhs.middleRows(d*N_V,N_V).noalias() -= (block.row(d)*physGradVel).transpose();
            // rhs: residual nonlinear
            for (short_t d = 0; d < dim; ++d)
                localRhs.middleRows(d*N_V,N_V).noalias() -= weight*J*
                        (physJacCurVel.row(d)*invJacALE*mdVelocity.values[0].col(q))(0,0) * basisValuesVel[0].col(q);
            // rhs: residual pressure
            for (short_t d = 0; d < dim; ++d)
                localRhs.middleRows(d*N_V,N_V).noalias() += weight * J * pressureValues.at(q) *
                       physGradVel.transpose() * invJacALE.col(d);
            // rhs: constraint residual
            localRhs.middleRows(dim*N_V,N_P).noalias() += weight * J * basisValuesPres.col(q) *
                    (physJacCurVel*invJacALE).trace();
        }
    }

    // auxiliary matrix nabla_U*F + F^T*nabla_u^T
    void setUFFU(gsMatrix<T> & UFFU, const gsMatrix<T> & F, const gsVector<T> & bGrad, short_t d)
    {
        UFFU.setZero(dim,dim);
        UFFU.col(d) = bGrad;
        UFFU = F.transpose()*UFFU;
        UFFU = UFFU + UFFU.transpose();
    }

protected:
    // problem info
    short_t dim;
    const gsPoissonPde<T> * pde_ptr;
    bool assembleMatrix;
    // flag to apply SUPG stabilization
    bool supg;
    // switch between assembling Newton (true) or Oseen (false) iteration
    index_t iterationType;
    index_t patch; // current patch
    // constants: viscosity and the force scaling factor
    T viscosity, forceScaling;
    // geometry mapping
    gsMapData<T> md;
    // local components of the global linear system
    gsMatrix<T> localMat;
    gsMatrix<T> localRhs;
    // local indices (at the current patch) of basis functions active at the current element
    gsMatrix<unsigned> localIndicesVel;
    gsMatrix<unsigned> localIndicesPres;
    // number of velocity and pressure basis functions active at the current element
    index_t N_V, N_P;
    // values and derivatives of velocity basis functions at quadrature points at the current element
    // values are stored as a N_V x numQuadPoints matrix; not sure about derivatives, must be smth like N_V*dim x numQuadPoints
    // if supg is true, also stores the hessian
    std::vector<gsMatrix<T> > basisValuesVel;
    // values of pressure basis functions active at the current element;
    // stores as a N_P x numQuadPoints matrix
    gsMatrix<T> basisValuesPres;
    // if supg is true, stores the derivatives
    gsMatrix<T> basisGradsPres;
    // RHS values at quadrature points at the current element; stored as a dim x numQuadPoints matrix
    gsMatrix<T> forceValues;
    // current displacement field
    const gsMultiPatch<T> & velocity;
    // evaluation data of the current velocity field
    gsMapData<T> mdVelocity;
    // current pressure field
    const gsMultiPatch<T> & pressure;
    // pressure values at the current element; stored as a 1 x numQuadPoints matrix
    gsMatrix<T> pressureValues;
    // pressure gradients at the current element (only for supg); stored as a dim x numQuadPoints matrix
    gsMatrix<T> pressureGrads;
    // ALE displacement field (only used for FSI)
    gsMultiPatch<T> dispALE;
    // flag to assemble ALE formulation
    bool ALE;
    // evaluation data of the ALE field
    gsMapData<T> mdALE;
};

} // namespace gismo
