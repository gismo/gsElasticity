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

}

template<class T>
void gsElPoissonAssembler<T>::assemble()
{
    GISMO_ASSERT(m_system.initialized(),
                 "Sparse system is not initialized, call initialize() or refresh()");

    // Reserve sparse system
    m_system.reserve(m_bases[0], m_options, this->pde().numRhs());

    // Compute the Dirichlet Degrees of freedom (if needed by m_options)
    gsAssembler<T>::computeDirichletDofs();

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
    for (index_t i = 0; i != localBoundaryIndices.size(); ++i)
    {
        index_t globalIndex = mapper.index(localBoundaryIndices(i),targetPatch);
        m_system.rhs().row(globalIndex) = PP * ddofs.row(i);
    }
}

template <class T>
void gsElPoissonAssembler<T>::setDirichletDoFs(const gsMultiPatch<> & sourceGeometry,
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
    {
        gsInfo << "Doesn't look like matching boundaries!\n";
        gsInfo << "Source " << sourcePatch << " " << sourceSide << std::endl;
        gsInfo << "Target " << targetPatch << " " << targetSide << std::endl;
    }
}

template <class T>
void gsElPoissonAssembler<T>::setDirichletDoFs(const gsField<> & sourceField,
                                               int sourcePatch, const boxSide & sourceSide,
                                               int targetPatch, const boxSide & targetSide)
{
    int matchingMode = checkMatchingBoundaries(*(sourceField.patch(sourcePatch).boundary(sourceSide)),
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
    {
        gsInfo << "Doesn't look like matching boundaries!\n";
        gsInfo << "Source " << sourcePatch << " " << sourceSide << std::endl;
        gsInfo << "Target " << targetPatch << " " << targetSide << std::endl;
    }
}

template <class T>
int gsElPoissonAssembler<T>::checkMatchingBoundaries(const gsGeometry<> & sourceBoundary,
                                                     const gsGeometry<> & targetBoundary)
{
    const T absTol = 1e-6; // another magic number

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
    {
        return 1;
    }
    else if ((sourceStart-targetEnd).norm() + (sourceEnd-targetStart).norm() < absTol)
    {
        return -1;
    }
    else
    {
        gsInfo << "Source start\n" << sourceStart << std::endl << "Source end\n" << sourceEnd << std::endl;
        gsInfo << "Target start\n" << targetStart << std::endl << "Target end\n" << targetEnd << std::endl;
        return 0;
    }
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


}// namespace gismo
