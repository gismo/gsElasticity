/** @file gsElPoissonAssembler.hpp

    @brief Provides assembler implementation for the gsElPoissonAssembler.

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

   Author(s): A. Shamanskiy
*/

#include <gsThermoElasticity/gsElPoissonAssembler.h>

#include <gsCore/gsField.h>

#include <gsThermoElasticity/gsVisitorElPoisson.h> // Stiffness volume integrals
#include <gsThermoElasticity/gsElPoissonPde.h>
#include <gsAssembler/gsVisitorNeumann.h> // Neumann boundary integrals
#include <gsAssembler/gsVisitorNitsche.h> // Nitsche boundary integrals
#include <gsAssembler/gsVisitorDg.h>      // DG interface integrals


namespace gismo
{

template <class T>
gsElPoissonAssembler<T>::gsElPoissonAssembler(const gsMultiPatch<T> & patches,
                                              const gsMultiBasis<T> & bases,
                                              const gsBoundaryConditions<T> & bcInfo,
                                              const gsFunction<T> & force,
                                              T conductivity,
                                              dirichlet::strategy dirStrategy,
                                              iFace::strategy intStrategy)
    : gsPoissonAssembler<T>(patches,bases,bcInfo,force,dirStrategy,intStrategy)
{
    m_options.setInt("DirichletStrategy", dirStrategy);
    m_options.setInt("InterfaceStrategy", intStrategy);

    typename gsPde<T>::Ptr pde( new gsElPoissonPde<T>(patches,bcInfo,force,conductivity) );
    gsAssembler<T>::initialize(pde, bases, m_options);

    m_system.reserve(m_bases[0], m_options, this->pde().numRhs());
    gsAssembler<T>::computeDirichletDofs();

}

template<class T>
void gsElPoissonAssembler<T>::assemble()
{
    GISMO_ASSERT(m_system.initialized(),
                 "Sparse system is not initialized, call initialize() or refresh()");

    // Reserve sparse system
    //m_system.reserve(m_bases[0], m_options, this->pde().numRhs());

    // Compute the Dirichlet Degrees of freedom (if needed by m_options)
    //gsAssembler<T>::computeDirichletDofs();

    // Clean the sparse system
   // m_system.setZero(); //<< this call leads to a quite significant performance degrade!

    // Assemble volume integrals
    gsAssembler<T>::template push<gsVisitorElPoisson<T> >();

    // Enforce Neumann boundary conditions
    gsAssembler<T>::template push<gsVisitorNeumann<T> >(m_pde_ptr->bc().neumannSides() );

    const int dirStr = m_options.getInt("DirichletStrategy");

    // If requested, enforce Dirichlet boundary conditions by Nitsche's method
    if ( dirStr == dirichlet::nitsche )
        gsAssembler<T>::template push<gsVisitorNitsche<T> >(m_pde_ptr->bc().dirichletSides());

     // If requested, enforce Dirichlet boundary conditions by diagonal penalization
     else if ( dirStr == dirichlet::penalize )
         penalizeDirichletDofs();

    if ( m_options.getInt("InterfaceStrategy") == iFace::dg )
        gsAssembler<T>::template pushInterface<gsVisitorDg<T> >();

    // Assembly is done, compress the matrix
    gsAssembler<T>::finalize();
}


template <class T>
void gsElPoissonAssembler<T>::penalizeDirichletDofs()
{
    GISMO_ASSERT( m_options.getInt("DirichletStrategy")
                  == dirichlet::penalize, "Incorrect options");

    const gsMultiBasis<T> & mbasis = m_bases[m_system.colBasis(0)];
    const gsDofMapper     & mapper = m_system.colMapper(0);
    // Note: dofs in the system, however dof values need to be computed
    const gsDofMapper & bmap = mbasis.getMapper(dirichlet::elimination,
                       static_cast<iFace::strategy>(m_options.getInt("InterfaceStrategy")),
                                                m_pde_ptr->bc(), 0) ;

    GISMO_ENSURE( m_ddof[0].rows() != 0 &&
                  m_ddof[0].cols() == m_pde_ptr->numRhs(),
                  "The Dirichlet DoFs were not computed.");

    // BCs
    for ( typename gsBoundaryConditions<T>::const_iterator
          it = m_pde_ptr->bc().dirichletBegin();
          it != m_pde_ptr->bc().dirichletEnd(); ++it )
    {
        const gsBasis<T> & basis = mbasis[it->patch()];

        gsMatrix<unsigned> bnd = basis.boundary(it->side() );
        for (index_t k=0; k!= bnd.size(); ++k)
        {
            // free dof position
            const index_t ii = mapper.index ( bnd(k) , it->patch() );
            // boundary dof position
            const index_t bb = bmap  .bindex( bnd(k) , it->patch() );

            m_system.matrix()(ii,ii) = PP;
            m_system.rhs().row(ii)   = PP * m_ddof[0].row(bb);
        }
    }

    // Corner values
    for ( typename gsBoundaryConditions<T>::const_citerator
          it = m_pde_ptr->bc().cornerBegin();
          it != m_pde_ptr->bc().cornerEnd(); ++it )
    {
        const int i  = mbasis[it->patch].functionAtCorner(it->corner);
        const int ii = mapper.bindex( i , it->patch );
        m_system.matrix()(ii,ii)       = PP;
        m_system.rhs().row(ii).setConstant(PP * it->value);
    }
}

template <class T>
void gsElPoissonAssembler<T>::setDirichletDoFs(const gsMatrix<> & ddofs, int targetPatch, const boxSide & targetSide)
{
    const gsDofMapper & mapper = m_system.colMapper(0);
    const gsMatrix<unsigned> & localBoundaryIndices = (m_bases[0])[targetPatch].boundary(targetSide);

    if (m_options.getInt("DirichletStrategy") == dirichlet::penalize)
    {
        for (index_t i = 0; i != localBoundaryIndices.size(); ++i)
        {
            index_t globalIndex = mapper.index(localBoundaryIndices(i),targetPatch);
            m_system.rhs().row(globalIndex) = PP * ddofs.row(i);
        }
    }
    else if (m_options.getInt("DirichletStrategy") == dirichlet::elimination)
    {
        gsMatrix<unsigned> systemIndices;
        mapper.localToGlobal(localBoundaryIndices,targetPatch,systemIndices);
        for (index_t i = 0; i != systemIndices.size(); ++i)
        {
            index_t dirichletIndex = mapper.global_to_bindex(systemIndices(i));
            m_ddof[0](dirichletIndex,0) = ddofs(i,0);
        }
    }
}

template <class T>
void gsElPoissonAssembler<T>::addDirichletData(const gsMultiPatch<> & sourceGeometry,
                                               const gsMultiPatch<> & sourceSolution,
                                               int sourcePatch, const boxSide & sourceSide,
                                               int targetPatch, const boxSide & targetSide)
{
    int matchingMode = checkMatchingBoundaries(*(sourceGeometry.patch(sourcePatch).boundary(sourceSide)),
                                               *(m_pde_ptr->patches().patch(sourcePatch).boundary(sourceSide)));

    if (matchingMode == 1)
    {
        setDirichletDoFs(sourceSolution.piece(sourcePatch).boundary(sourceSide)->coefs(),targetPatch,targetSide);
    }
    else if (matchingMode == -1)
    {
        setDirichletDoFs(sourceSolution.piece(sourcePatch).boundary(sourceSide)->coefs().colwise().reverse(),targetPatch,targetSide);
    }
    else
        gsInfo << "Doesn't look like matching boundaries!\n"
               << "Source: p " << sourcePatch << " s " << sourceSide << std::endl
               << "Target: p " << targetPatch << " s " << targetSide << std::endl;
}

template <class T>
void gsElPoissonAssembler<T>::addDirichletData(const gsField<> & sourceField,
                                               int sourcePatch, const boxSide & sourceSide,
                                               int targetPatch, const boxSide & targetSide)
{
    int matchingMode = checkMatchingBoundaries(*(sourceField.patches().patch(sourcePatch).boundary(sourceSide)),
                                               *(m_pde_ptr->patches().patch(targetPatch).boundary(targetSide)));
    if (matchingMode == 1)
    {
        setDirichletDoFs(sourceField.igaFunction(sourcePatch).boundary(sourceSide)->coefs(),targetPatch,targetSide);
    }
    else if (matchingMode == -1)
    {
        setDirichletDoFs(sourceField.igaFunction(sourcePatch).boundary(sourceSide)->coefs().colwise().reverse(),targetPatch,targetSide);
    }
    else
        gsInfo << "Doesn't look like matching boundaries!\n"
               << "Source: p " << sourcePatch << " s " << sourceSide << std::endl
               << "Target: p " << targetPatch << " s " << targetSide << std::endl;
}

template <class T>
int gsElPoissonAssembler<T>::checkMatchingBoundaries(const gsGeometry<> & sourceBoundary,
                                                     const gsGeometry<> & targetBoundary)
{
    const T absTol = 1e-6; // another magic number

    int sDim = sourceBoundary.dimensions().first;
    int tDim = targetBoundary.dimensions().first;

    if (sDim == 1 && tDim == 1)
    {
        gsMatrix<> startPoint;
        startPoint.resize(1,1);
        startPoint << 0.;
        gsMatrix<> endPoint;
        endPoint.resize(1,1);
        endPoint << 1.;

        gsMatrix<> sourceStart = sourceBoundary.eval(startPoint);
        gsMatrix<> targetStart = targetBoundary.eval(startPoint);

        gsMatrix<> sourceEnd = sourceBoundary.eval(endPoint);
        gsMatrix<> targetEnd = targetBoundary.eval(endPoint);

        if ((sourceStart-targetStart).norm() + (sourceEnd-targetEnd).norm() < absTol)
            return 1;
        else if ((sourceStart-targetEnd).norm() + (sourceEnd-targetStart).norm() < absTol)
            return -1;
        else
            return 0;
    }
    else if (sDim == 2 && tDim == 2)
    {
        gsInfo << "Accurate orientation check for two 3D is not implemented\n";
        return 1;
    }

    return 0;
}

template <class T>
void gsElPoissonAssembler<T>::setUnitingConstraint(const boxSide & side, bool verbosity)
{
    const gsDofMapper & mapper = m_system.colMapper(0);
    const gsMatrix<unsigned> & bdryIndices = (m_bases[0])[0].boundary(side);

    long long int numConstraints = bdryIndices.rows() - 2;
    long long int size = m_system.rhs().rows();
    m_system.rhs().conservativeResize(size+numConstraints,1);
    m_system.rhs().block(size,0,numConstraints,1).setZero();

    //m_system.matrix().reserve((size+numConstraints)*(size+numConstraints));
    m_system.matrix().conservativeResize(size+numConstraints,size+numConstraints);
    //m_system.matrix().resize(size+numConstraints,size+numConstraints);
    for (int i = 0; i < numConstraints; ++i)
    {
        index_t globalIndex1 = mapper.index(bdryIndices(i+1),0);
        index_t globalIndex2 = mapper.index(bdryIndices(i+2),0);

        if (verbosity)
            gsInfo << "I1 " << bdryIndices(i+1) << " G1 " << globalIndex1 << " I2 " << bdryIndices(i+2) << " G2 " << globalIndex2 << std::endl;

        m_system.matrix().coeffRef(size + i,globalIndex1) = 1.;
        m_system.matrix().coeffRef(globalIndex1,size + i) = 1.;

        m_system.matrix().coeffRef(size + i,globalIndex2) = -1.;
        m_system.matrix().coeffRef(globalIndex2,size + i) = -1.;
    }
}

template <class T>
void gsElPoissonAssembler<T>::addNeummannData(const gsMultiPatch<> & sourceGeometry,
                                              const gsMultiPatch<> & sourceSolution,
                                              int sourcePatch, const boxSide & sourceSide,
                                              int targetPatch, const boxSide & targetSide)
{

}

template <class T>
void gsElPoissonAssembler<T>::addNeummannData(const gsField<> & sourceField,
                                              int sourcePatch, const boxSide & sourceSide,
                                              int targetPatch, const boxSide & targetSide)
{

}

template <class T>
gsFluxFunction<T>::gsFluxFunction(int sourcePatch, boundary::side sourceSide,
                                  gsMultiPatch<T> const & sourceGeo,
                                  gsMultiPatch<T> const & sourceSolution,
                                  T alpha)
    : m_patch(sourcePatch),
      m_side(sourceSide),
      m_geo(sourceGeo),
      m_sol(sourceSolution),
      m_alpha(alpha)
{

}

template <class T>
void gsFluxFunction<T>::eval_into(gsMatrix<T> const & u, gsMatrix<T> & result) const
{
    gsMatrix<T> params;
    m_geo.patch(m_patch).invertPoints(u,params);

    typename gsGeometry<T>::Evaluator geoEval(
                m_geo.patch(m_patch).evaluator(NEED_VALUE|NEED_JACOBIAN|NEED_GRAD_TRANSFORM));
    geoEval->evaluateAt(params);
    gsMatrix<T> grads;
    m_sol.patch(m_patch).deriv_into(params,grads);

    gsMatrix<T> physGrad;
    gsVector<T> normal;
    result.resize(1,params.cols());
    for (int i = 0; i < params.cols(); ++i)
    {
        geoEval->outerNormal(i,m_side,normal);
        geoEval->transformGradients(i,grads,physGrad);
        result.at(i) = (-1*m_alpha*physGrad.transpose()*normal/normal.norm())(0,0);
    }
}

template <class T>
gsGradFunction<T>::gsGradFunction(int sourcePatch,gsMultiPatch<T> const & sourceSolution)
    : m_patch(sourcePatch),
      m_sol(sourceSolution)
{

}

template <class T>
void gsGradFunction<T>::eval_into(gsMatrix<T> const & u, gsMatrix<T> & result) const
{
    m_sol.patch(m_patch).deriv_into(u,result);
}

template <class T>
gsPhysGradFunction<T>::gsPhysGradFunction(int sourcePatch,
                                          gsMultiPatch<T> const & sourceGeo,
                                          gsMultiPatch<T> const & sourceSolution)
    : m_patch(sourcePatch),
      m_geo(sourceGeo),
      m_sol(sourceSolution)
{

}

template <class T>
void gsPhysGradFunction<T>::eval_into(gsMatrix<T> const & u, gsMatrix<T> & result) const
{
    typename gsGeometry<T>::Evaluator geoEval(
                m_geo.patch(m_patch).evaluator(NEED_GRAD_TRANSFORM));
    geoEval->evaluateAt(u);

    gsMatrix<T> grads;
    m_sol.patch(m_patch).deriv_into(u,grads);
    result.resize(targetDim(),u.cols());
    gsMatrix<T> physGrad;
    for (int i = 0; i < u.cols(); ++i)
    {
        geoEval->transformGradients(i,grads,physGrad);
        result(0,i) = physGrad.at(0);
        result(1,i) = physGrad.at(1);
    }
}

template <class T>
gsDetFunction<T>::gsDetFunction(int sourcePatch, gsMultiPatch<T> const & sourceGeo)
    : m_patch(sourcePatch),
      m_geo(sourceGeo)
{

}

template <class T>
void gsDetFunction<T>::eval_into(gsMatrix<T> const & u, gsMatrix<T> & result) const
{
    typename gsGeometry<T>::Evaluator geoEval(
                m_geo.patch(m_patch).evaluator(NEED_MEASURE));
    geoEval->evaluateAt(u);
    result.resize(1,u.cols());
    gsInfo << "GradsTR\n";
    gsMatrix<T> grads;
    m_geo.patch(m_patch).deriv_into(u,grads);
    gsInfo << grads << std::endl;
    gsInfo << jac
    for (int i = 0; i < u.cols(); ++i)
    {
        result(0,i) = geoEval->jacDet(i);
    }

}


}// namespace gismo
