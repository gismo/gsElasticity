/** @file gsElasticityAssembler.hpp

    @brief Provides nonlinear elasticity system matrices for 2D plain strain and 3D continua.

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): O. Weeger, A. Shamanskiy (TU Kaiserslautern)
*/

#pragma once

#include <gsElasticity/gsElasticityAssembler.h>

#include <gsAssembler/gsGaussRule.h>
#include <gsCore/gsMultiBasis.h>
#include <gsCore/gsDomainIterator.h>
#include <gsCore/gsField.h>
#include <gsUtils/gsPointGrid.h>
#include <gsPde/gsPoissonPde.h>

// Element visitors
#include <gsElasticity/gsVisitorLinearElasticity.h>
//#include <gsElasticity/gsVisitorNonLinElasticity.h>
#include <gsElasticity/gsVisitorElasticityNeumann.h>
//#include <gsElasticity/gsVisitorElasticityPressure.h>
//#include <gsElasticity/gsVisitorMassElasticity.h>

// ---
#include <gsCore/gsBoxTopology.h>
#include <gsNurbs/gsTensorBSplineBasis.h>
#include <gsNurbs/gsBSplineBasis.h>
#include <gsNurbs/gsKnotVector.h>


namespace gismo
{

template<class T>
gsElasticityAssembler<T>::gsElasticityAssembler(gsMultiPatch<T> const & patches,
                                                gsMultiBasis<T> const & basis,
												T E_modulus,
												T poissons_ratio,
												T density_rho,
												gsBoundaryConditions<T> const & bconditions,
                                                const gsFunction<T> & body_force,
                                                dirichlet::values computeStrategy,
                                                dirichlet::strategy enforceStrategy,
                                                gsFunction<T> * E,
                                                gsFunction<T> * pr)
:   Base(),
	m_rho(density_rho),
    m_bConditions(bconditions),
    m_bodyForce(&body_force),
    f_E(E),
    f_pr(pr)
{
    m_options.setInt("DirichletValues"  , computeStrategy);
    m_options.setInt("DirichletStrategy", enforceStrategy);

    // Copy the spline basis
    m_dim = body_force.targetDim();
    for (int d = 0; d < m_dim; ++d)
        m_bases.push_back(basis);

    gsPiecewiseFunction<T> pwf;
    pwf.addPiece(body_force);
    typename gsPde<T>::Ptr pde( new gsPoissonPde<T>(patches,bconditions,pwf) );
    Base::initialize(pde, m_bases, m_options);

    // Refresh .....
	GISMO_ASSERT(m_dim == patches.parDim(), "Dimensions not matching!");
	GISMO_ASSERT(m_dim == m_bases.front().dim(), "Dimensions not matching!");
	GISMO_ASSERT(m_dim == 2 || m_dim == 3, "Only for 2D and 3D!");
    std::vector<gsDofMapper> m_dofMappers(m_dim);
    for (index_t i = 0; i < m_dim; i++)
        m_bases.front().getMapper(enforceStrategy,iFace::glue,m_bConditions,m_dofMappers[i],i,true);

    gsVector<unsigned> dims;
    dims.setOnes(m_dim);
    m_system = gsSparseSystem<T>(m_dofMappers, dims);

    // Initialize material properties
    m_lambda = E_modulus * poissons_ratio / ( (1.+poissons_ratio)*(1.-2.*poissons_ratio)) ;
    m_mu     = E_modulus / (2.*(1.+poissons_ratio)) ;

	// Initialize factors for time-dependent external forces
	m_tfac_neumann = 1.0;
	m_tfac_force = 1.0;

	m_MATERIAL_LAW = 1;
}

template<class T>
void gsElasticityAssembler<T>::set_tfac(const T tfac_neumann,
		                                const T tfac_force)
{
	m_tfac_neumann = tfac_neumann;
	m_tfac_force = tfac_force;
}

/// Sets the \em m_MATERIAL_LAW to 0: St. Venant-Kirchhoff, 1: Neo-Hooke
template<class T>
void gsElasticityAssembler<T>::set_MaterialLaw(const int material)
{
    GISMO_ASSERT( (material >= 0) && (material <= 1), "Only 0 (St. Venant-Kirchhoff) or 1 (Neo-Hooke) allowed!");
	m_MATERIAL_LAW = material;
}

template<class T>
void gsElasticityAssembler<T>::assembleNeumann()
{
    for ( typename gsBoundaryConditions<T>::const_iterator it = m_bConditions.neumannBegin();
          it != m_bConditions.neumannEnd();
		  ++it )
    {
        gsVisitorElasticityNeumann<T> neumann(*it->function(), it->side(), m_tfac_neumann);
        Base::template push<gsVisitorElasticityNeumann<T> >(neumann, *it);
    }
}

template<class T>
void gsElasticityAssembler<T>::assemble()
{
    gsWarn << " USING extensions/gsElasticity.\n";
    //computeDirichletDofs();

    if ( this->numDofs() == 0 ) // Are there any interior dofs ?
    {
        gsWarn << " No internal DOFs. Computed Dirichlet boundary only.\n" <<"\n" ;
        return;
    }

    // Pre-allocate non-zero elements for each column of the
    // sparse matrix
    size_t nonZerosPerCol = m_dim;
    for (index_t i = 0; i < m_dim; ++i) // to do: improve
        nonZerosPerCol *= 2 * m_bases.front().maxDegree(i) + 1;
    m_system.reserve(nonZerosPerCol, 1);

    // Compute the Dirichlet Degrees of freedom (if needed by m_options)
    for (int i = 0; i < m_dim; ++i)
        Base::computeDirichletDofs(i);

    // Assemble volume stiffness and load vector integrals
    gsVisitorLinearElasticity<T> visitor(m_lambda, m_mu, m_rho, *m_bodyForce, m_tfac_force,f_E,f_pr);
    Base::template push<gsVisitorLinearElasticity<T> >(visitor);

    // Enforce Neumann boundary conditions
    assembleNeumann();

    // Assembly is done, compress the matrix
    m_system.matrix().makeCompressed();
}

template<class T>
void gsElasticityAssembler<T>::assemble(const gsMultiPatch<T> & deformed)
{
/*
    std::cout << "Nonlinear elasticity: assemble residual and tangential matrix." << std::endl;

    m_matrix.resize(this->numDofs(),this->numDofs());
    m_matrix.setZero();

    m_system.rhs().resize(this->numDofs(),1);
    m_system.rhs().setZero();


    gsVisitorNonLinElasticity<T> visitor(m_lambda, m_mu, m_rho, *m_bodyForce, deformed.patch(0), m_tfac_force);
	//gsVisitorNonLinElasticity<T> visitor(m_lambda, m_mu, m_rho, *m_bodyForce, deformed.patch(0));

	visitor.set_MaterialLaw(m_MATERIAL_LAW);

    // Assemble volume stiffness and load vector integrals
    for (index_t np=0; np < m_pde_ptr->domain().nPatches(); ++np )
    {
        visitor.setDeformed( deformed.patch(np) );

        //Assemble stiffness matrix and rhs for the local patch
        // with index np and add to m_matrix and m_system.rhs()
        this->apply(visitor, np);
    }

    // Enforce Neumann forces
    assembleNeumann();

    // Assembly is done, compress the matrix
    m_matrix.makeCompressed();


}

size_t sumBoundarySize(size_t a, gsDofMapper b) { return a + b.boundarySize(); }

template<class T>
void gsElasticityAssembler<T>::computeDirichletDofs()
{
    size_t numDDofs = std::accumulate(m_dofMappers.begin(), m_dofMappers.end(), (size_t)0, sumBoundarySize);
    m_ddof.setZero(numDDofs,1);

    if (m_options.dirValues == dirichlet::interpolation)
    {
        computeDirichletDofsIntpl();
    }
    else if (m_options.dirValues == dirichlet::l2Projection)
    {
        computeDirichletDofsL2Proj();
    }
    else
    {
        gsInfo << "This Dirichlet value strategy is not implemented\n";
    }
}


// WARNING: WORKS ONLY FOR TENSOR-PRODUCT-BASES!
template<class T>
void gsElasticityAssembler<T>::computeDirichletDofsIntpl()
{
    // Iterate over all patch-sides with Dirichlet-boundary conditions
    for ( typename gsBoundaryConditions<T>::const_iterator it = m_bConditions.dirichletBegin();
          it != m_bConditions.dirichletEnd();
		  ++it )
    {
        if ( it->isHomogeneous() )
			continue;

		const int unk = it->unknown();
        const int k   = it->patch();
        const gsBasis<T> & basis = (m_bases[unk/m_dim])[k];

		const gsDofMapper &mapper = m_dofMappers[unk];

        // Get dofs on this boundary
        gsMatrix<unsigned> boundary = basis.boundary(it->side()) ;

        // Get the side information
        int dir = it->side().direction( );
        index_t param = (it->side().parameter() ? 1 : 0);

        // Compute grid of points on the face ("face anchors")
        std::vector< gsVector<T> > rr;
        rr.reserve( m_dim );

        for ( int i = 0; i < m_dim; ++i)
        {
            if ( i == dir )
            {
                gsVector<T> b(1);
                b[0] = ( basis.component(i).support() ) (0, param);
                rr.push_back(b);
            }
            else
            {
                rr.push_back( basis.component(i).anchors().transpose() );
            }
        }

        // Compute dirichlet values
        gsMatrix<T> fpts =
            it->function()->eval( m_pde_ptr->domain()[it->patch()].eval(  gsPointGrid<T>( rr ) ) );

        // Interpolate dirichlet boundary
        typename gsBasis<T>::uPtr h = basis.boundaryBasis(it->side());
        typename gsGeometry<T>::uPtr geo = h->interpolateAtAnchors(fpts);
        const gsMatrix<T> & dVals =  geo->coefs();

        // Save corresponding boundary dofs
        for (index_t k2=0; k2!= boundary.size(); ++k2)
        {
            const int ii= mapper.bindex( (boundary)(k2) , it->patch() );
            m_ddof.row(ii) = m_tfac_neumann * dVals.row(k2);
        }
    }
*/
}

template<class T>
void  gsElasticityAssembler<T>::reComputeDirichletDofs(gsMultiPatch<T> &deformed)
{
    /*
	computeDirichletDofsIntpl();

	// -> deformed must be set using updated m_ddof! How?
	for (index_t p=0; p < m_pde_ptr->domain().nPatches(); ++p )
    {
        // Update displacement solution coefficients on patch p
        const int sz  = m_bases[0][p].size();

        gsMatrix<T> & coeffs = deformed.patch(p).coefs();

        for (index_t j = 0; j < m_dim; ++j)
        {
            const gsDofMapper & mapper = m_dofMappers[j];
            for (index_t i = 0; i < sz; ++i)
            {
				if ( !mapper.is_free(i, p) ) // eliminated DoF: fill with Dirichlet data
                {
                    coeffs(i,j) = m_ddof(mapper.bindex(i, p), 0);
                }
            }
        }
	}
    */

}


template<class T>
void  gsElasticityAssembler<T>::setSolution(const gsMatrix<T>& solVector,
                                            gsMultiPatch<T>& result) const
{
    GISMO_ASSERT(this->numDofs() == m_system.rhs().rows(), "Something went wrong, assemble() not called?");

    std::vector<gsFunction<T> * > sols ;

    for (index_t p=0; p < m_pde_ptr->domain().nPatches(); ++p )
    {
        // Update solution coefficients on patch p
        const int sz  = m_bases[0][p].size();

        gsMatrix<T> & coeffs = result.patch(p).coefs();
		coeffs.setZero(); //OWEE

        for (index_t j = 0; j < m_dim; ++j)
        {
            const gsDofMapper & mapper = m_system.colMapper(j);
            for (index_t i = 0; i < sz; ++i)
            {
				if ( mapper.is_free(i, p) ) // DoF value is in the solVector ?
                {
                    coeffs(i,j) += solVector( mapper.index(i, p), 0);
                }
                else // eliminated DoF: fill with Dirichlet data
                {
                    coeffs(i,j) = m_ddof[0](mapper.bindex(i, p), 0);
                }
            }
        }
    }
}


template<class T>
void  gsElasticityAssembler<T>::updateSolution(const gsMatrix<T>& solVector,
                                               gsMultiPatch<T>& result) const
{
    GISMO_ASSERT(this->numDofs() == m_system.rhs().rows(), "Something went wrong, assemble() not called?");

    std::vector<gsFunction<T> * > sols ;

    for (index_t p=0; p < m_pde_ptr->domain().nPatches(); ++p )
    {
        // Update solution coefficients on patch p
        const int sz  = m_bases[0][p].size();

        gsMatrix<T> & coeffs = result.patch(p).coefs();

        for (index_t j = 0; j < m_dim; ++j)
        {
            const gsDofMapper & mapper =  m_system.colMapper(j);
            for (index_t i = 0; i < sz; ++i)
            {
                if ( mapper.is_free(i, p) ) // DoF value is in the solVector
                {
                    coeffs(i,j) += solVector( mapper.index(i, p), 0);
                }
            }
        }
    }
}

/*template<class T>
void  gsElasticityAssembler<T>::constructSolution(const gsMatrix<T>& solVector,
                                                  gsMultiPatch<T>& result) const
{
    GISMO_ASSERT(this->numDofs() == m_system.rhs().rows(), "Something went wrong, assemble() not called?");

    // empty the multipatch container
    result.clear();

	gsMatrix<T> coeffs;

    for (index_t p=0; p < m_pde_ptr->domain().nPatches(); ++p )
    {
        // Update solution coefficients on patch p
        const int sz  = m_bases[0][p].size(); // m_pde_ptr->domain()[p].size();

        coeffs.resize(sz, m_dim);

        for (index_t j = 0; j < m_dim; ++j) // For all components x, y, z
        {
            const gsDofMapper & mapper =  m_system.colMapper(j);

            for (index_t i = 0; i < sz; ++i)
            {
                if ( mapper.is_free(i, p) ) // DoF value is in the solVector ?
                {
                    coeffs(i,j) = solVector( mapper.index(i, p), 0);
                }
                else // eliminated DoF: fill with Dirichlet data
                {
                    coeffs(i,j) = m_ddof[0](mapper.bindex(i, p), 0);
                }
            }
        }

		result.addPatch( m_bases[0][p].makeGeometry( give(coeffs) ) );
    }
}*/


template<class T>
void gsElasticityAssembler<T>::assembleMass()
{
    /*
    std::cout << "Linear Elasticity: assemble mass matrix." << std::endl;

	//computeDirichletDofs();
    index_t numDirichlet = 0;
    for (index_t i = 0; i < m_dim; ++i)
        numDirichlet +=  m_system.colMapper(i).boundarySize();
    m_ddof.setZero(numDirichlet, 1);

    if (this->numDofs() == 0 ) // Are there any interior dofs ?
    {
        gsWarn << " No internal DOFs. Computed Dirichlet boundary only.\n" <<"\n" ;
        return;
    }

    // Pre-allocate non-zero elements for each column of the
    // sparse matrix
    size_t nonZerosPerCol = m_dim;
    for (index_t i = 0; i < m_dim; ++i) // to do: improve
        nonZerosPerCol *= 2 * m_bases.front().maxDegree(i) + 1;

    m_matrix = gsSparseMatrix<T>(this->numDofs(), this->numDofs()); // Clean matrices
    m_matrix.reserve( gsVector<int>::Constant(this->numDofs(), nonZerosPerCol) );

    // Resize the load vector
    m_system.rhs().setZero(this->numDofs(), 1 );

    // Assemble volume stiffness and load vector integrals
    gsVisitorMassElasticity<T> visitor(m_rho);
    for (index_t np=0; np < m_pde_ptr->domain().nPatches(); ++np )
    {
        //Assemble stiffness matrix and rhs for the local patch
        // with index np and add to m_matrix and m_system.rhs()
        this->apply(visitor, np);
    }

    // Assembly is done, compress the matrix
    m_matrix.makeCompressed();
    */
}


template<class T>
void gsElasticityAssembler<T>::computeStresses(
        const gsMatrix<T>& solVector,
        const gsMatrix<T>& u,
        int patchIndex,
        gsMatrix<T>& result,
        bool computeVonMises ) const
{
    /*
    //m_dim = m_basis.dim();
    size_t dimStrain = (m_dim*(m_dim+1))/2;

    result.resize( dimStrain + unsigned( computeVonMises), u.cols() );
    result.setZero();

    unsigned evFlags = NEED_VALUE | NEED_JACOBIAN | NEED_MEASURE | NEED_GRAD_TRANSFORM;
    typename gsGeometryEvaluator<T>::uPtr geoEval(getEvaluator(evFlags, m_pde_ptr->domain()[patchIndex]));

    // Evaluate basis functions on element
    gsMatrix<T> basisGrad( m_dim, 1 );
    gsMatrix<unsigned> actives;

    gsMatrix<T> physGrad( m_dim, 1 );
    gsMatrix<T> uPartDers( m_dim, m_dim );

    // Do this pointwise, since we have no assumptions
    // on the evaluation points.
    for( size_t k=0; k < size_t( u.cols() ); k++ )
    {
        // active functions at the evaluation point
        m_bases[0].basis(patchIndex).active_into( u.col(k), actives);
        // for each component, make global DoFs out of these local DoFs
        std::vector< gsMatrix<unsigned> > ci_actives(m_dim,actives);
        for( index_t ci = 0; ci < m_dim; ++ci )
            m_dofMappers[ci].localToGlobal(actives, patchIndex, ci_actives[ci]);

        // compute the gradients of the basis functions at the eval.pt.
        m_bases[0].basis(patchIndex).deriv_into( u.col(k), basisGrad);
        // and transform them to the physical space
        physGrad.resize(m_dim, actives.rows() );
        physGrad.setZero();
        geoEval->evaluateAt( u.col(k) );
        geoEval->transformGradients( index_t(0) , basisGrad, physGrad);

        // compute the partial derivatives of the
        // displacement field
        uPartDers.setZero();
        for( index_t dc=0; dc < m_dim; dc++) // component
            for( index_t dd=0; dd < m_dim; dd++) // derivative
                for( size_t ai=0; ai < size_t( actives.rows() ); ai++ ) // active function
                {
                    unsigned tmpi = ci_actives[ dc ](ai,0);
                    real_t tmpCoef;
                    if( m_dofMappers[dc].is_free_index( tmpi) )
                        tmpCoef = solVector( tmpi, 0);
                    else // dirichlet boundary
                        tmpCoef = m_ddof( m_dofMappers[dc].global_to_bindex(tmpi) );

                    uPartDers( dc, dd ) += tmpCoef * physGrad( dd, ai );
                }

        if( m_dim == 2 )
        {
            // sigma_{11}:
            result(0,k) = m_lambda * ( uPartDers(0,0)+uPartDers(1,1) ) + 2*m_mu*uPartDers(0,0);
            // sigma_{22}:
            result(1,k) = m_lambda * ( uPartDers(0,0)+uPartDers(1,1) ) + 2*m_mu*uPartDers(1,1);
            // sigma_{12}:
            result(2,k) = m_mu * ( uPartDers(0,1) + uPartDers(1,0) );

            if( computeVonMises )
            {
                result(3,k) = math::sqrt(
                            result(0,k)*result(0,k)
                            - result(0,k)*result(1,k)
                            + result(1,k)*result(1,k)
                            + 3.0 * result(2,k)*result(2,k) );
            }

        }
        else if( m_dim == 3 )
        {
            real_t tmp = m_lambda * ( uPartDers(0,0)+uPartDers(1,1)+uPartDers(2,2) );
            result(0,k) = tmp + 2*m_mu * uPartDers(0,0); // sigma_{11}
            result(1,k) = tmp + 2*m_mu * uPartDers(1,1); // sigma_{22}
            result(2,k) = tmp + 2*m_mu * uPartDers(2,2); // sigma_{33}
            result(3,k) = m_mu * ( uPartDers(0,1) + uPartDers(1,0) ); // sigma_{12}
            result(4,k) = m_mu * ( uPartDers(0,2) + uPartDers(2,0) ); // sigma_{13}
            result(5,k) = m_mu * ( uPartDers(2,1) + uPartDers(1,2) ); // sigma_{23}

            if( computeVonMises )
            {
                T tmp2 =   (result(0,k)-result(1,k))*(result(0,k)-result(1,k))
                        + (result(0,k)-result(2,k))*(result(0,k)-result(2,k))
                        + (result(1,k)-result(2,k))*(result(1,k)-result(2,k))
                        + 6.0 * ( result(3,k)*result(3,k) + result(4,k)*result(4,k) + result(5,k)*result(5,k) );

                result(6,k) = math::sqrt( 0.5 * tmp2 );
            }

        }
    }
    */
}

template <class T>
void gsElasticityAssembler<T>::constructStresses(const gsMatrix<T>& solVector,
                                                 gsMultiFunction<T>& result,
                                                 stress_type::type type) const
{
    result.clear();
    for (index_t i = 0; i < m_pde_ptr->domain().nPatches(); ++i )
    {
        result.addFunction(typename gsFunction<T>::uPtr(new gsStressFunction<T>(solVector,*this,i,type)));
        //gsFunction<T> * temp = new gsStressFunction<T>(solVector,*this,i,type);
        //result.addFunction(temp);
    }
}

template <class T>
const gsMatrix<T> & gsElasticityAssembler<T>::rhs()
{
    if (m_rhsExtra.rows() == 0)
    {
        return m_system.rhs();
    }
    else
    {
        return m_rhsExtra;
    }
}

template <class T>
void gsElasticityAssembler<T>::addNeummannData(const gsMultiPatch<> &sourceGeometry,
                                               const gsMultiPatch<> &sourceSolution,
                                               int sourcePatch, const boxSide &sourceSide,
                                               int targetPatch, const boxSide &targetSide)
{
    /*
    int matchingMode = checkMatchingBoundaries(*(sourceGeometry.patch(sourcePatch).boundary(sourceSide)),
                                               *(m_pde_ptr->domain().patch(targetPatch).boundary(targetSide)));
    if (matchingMode != 0)
    {
        if (m_rhsExtra.rows() == 0 && m_system.rhs().rows() != 0)
            m_rhsExtra = m_system.rhs();

        const gsMatrix<unsigned> & localBoundaryIndices = (m_bases[0])[targetPatch].boundary(targetSide);
        const gsMatrix<T> coefs = sourceSolution.piece(sourcePatch).boundary(sourceSide)->coefs();

        int numLayers = localBoundaryIndices.rows()/coefs.rows();

        std::map<unsigned,T> neumannDoFs;
        for (int l = 0; l < numLayers; ++l)
        {
            for (int i = 0; i < coefs.rows(); ++i)
            {
                if (matchingMode == 1)
                    neumannDoFs[localBoundaryIndices(l*coefs.rows()+i)] = coefs(i);
                if (matchingMode == -1)
                    neumannDoFs[localBoundaryIndices(l*coefs.rows()+i)] = coefs(coefs.rows() - i - 1);
            }

        }

        gsVisitorElasticityPressure<T> visitor(neumannDoFs,targetSide,m_rhsExtra);
        this->apply(visitor,targetPatch,targetSide);
    }
    else
        gsInfo << "Doesn't look like matching boundaries!\n"
               << "Source: p " << sourcePatch << " s " << sourceSide << std::endl
               << "Target: p " << targetPatch << " s " << targetSide << std::endl;
    */
}

template <class T>
void gsElasticityAssembler<T>::addNeummannData(const gsField<> &sourceField,
                                               int sourcePatch, const boxSide &sourceSide,
                                               int targetPatch, const boxSide &targetSide)
{
    /*
    int matchingMode = checkMatchingBoundaries(*(sourceField.patch(sourcePatch).boundary(sourceSide)),
                                               *(m_pde_ptr->domain().patch(targetPatch).boundary(targetSide)));
    if (matchingMode != 0)
    {
        if (m_rhsExtra.rows() == 0 && m_system.rhs().rows() != 0)
            m_rhsExtra = m_system.rhs();

        const gsMatrix<unsigned> & localBoundaryIndices = (m_bases[0])[targetPatch].boundary(targetSide);
        const gsMatrix<T> coefs = sourceField.igaFunction(sourcePatch).boundary(sourceSide)->coefs();

        int numLayers = localBoundaryIndices.rows()/coefs.rows();

        std::map<unsigned,T> neumannDoFs;
        for (int l = 0; l < numLayers; ++l)
        {
            for (int i = 0; i < coefs.rows(); ++i)
            {
                if (matchingMode == 1)
                    neumannDoFs[localBoundaryIndices(l*coefs.rows()+i)] = coefs(i);
                else if (matchingMode == -1)
                    neumannDoFs[localBoundaryIndices(l*coefs.rows()+i)] = coefs(coefs.rows() - i - 1);
            }
        }

        gsVisitorElasticityPressure<T> visitor(neumannDoFs,targetSide,m_rhs Extra);
        this->apply(visitor,targetPatch,targetSide);
    }
    else
        gsInfo << "Doesn't look like matching boundaries!\n"
               << "Source: p " << sourcePatch << " s " << sourceSide << std::endl
               << "Target: p " << targetPatch << " s " << targetSide << std::endl;
    */
 }

template <class T>
int gsElasticityAssembler<T>::checkMatchingBoundaries(const gsGeometry<> & sourceBoundary,
                                                      const gsGeometry<> & targetBoundary)
{
    const T absTol = 1e-6; // another magic number

    int sDim = sourceBoundary.dimensions().first;
    int tDim = targetBoundary.dimensions().first;

    if (sDim == 1 && tDim == 1)
    {
        gsMatrix<> startPoint;
        startPoint.resize(1,1);
        startPoint.at(0) = 0.;
        gsMatrix<> endPoint;
        endPoint.resize(1,1);
        endPoint.at(0) = 1.;

        gsMatrix<> sourceStart = sourceBoundary.eval(startPoint);
        gsMatrix<> sourceEnd = sourceBoundary.eval(endPoint);

        gsMatrix<> targetStart = targetBoundary.eval(startPoint);
        gsMatrix<> targetEnd = targetBoundary.eval(endPoint);

        if ((sourceStart-targetStart).norm() + (sourceEnd-targetEnd).norm() < absTol &&
             sourceBoundary.coefsSize() == targetBoundary.coefsSize())
            return 1;
        else if ((sourceStart-targetEnd).norm() + (sourceEnd-targetStart).norm() < absTol &&
                  sourceBoundary.coefsSize() == targetBoundary.coefsSize())
            return -1;
        else
            return 0;
    }
    else if (sDim == 1 && tDim == 2)
    {
        gsMatrix<> startPoint;
        startPoint.resize(1,1);
        startPoint.at(0) = 0.;
        gsMatrix<> endPoint;
        endPoint.resize(1,1);
        endPoint.at(0) = 1.;

        gsMatrix<> sourceStart = sourceBoundary.eval(startPoint);
        gsMatrix<> sourceEnd = sourceBoundary.eval(endPoint);

        if (targetBoundary.dimensions().first == 2)
        {
            startPoint.conservativeResize(2,1);
            startPoint(1) = 0.;
            endPoint.conservativeResize(2,1);
            endPoint(1) = 0.;

        }
        gsMatrix<> targetStart = targetBoundary.eval(startPoint);
        gsMatrix<> targetEnd = targetBoundary.eval(endPoint);

        real_t alignedCheck = (sourceStart.at(0)-targetStart.at(0))*(sourceStart.at(0)-targetStart.at(0)) +
                              (sourceStart.at(1)-targetStart.at(1))*(sourceStart.at(1)-targetStart.at(1)) +
                              (sourceEnd.at(0)-targetEnd.at(0))*(sourceEnd.at(0)-targetEnd.at(0)) +
                              (sourceEnd.at(1)-targetEnd.at(1))*(sourceEnd.at(1)-targetEnd.at(1));

        real_t reverseCheck = (sourceStart.at(0)-targetEnd.at(0))*(sourceStart.at(0)-targetEnd.at(0)) +
                              (sourceStart.at(1)-targetEnd.at(1))*(sourceStart.at(1)-targetEnd.at(1)) +
                              (sourceEnd.at(0)-targetStart.at(0))*(sourceEnd.at(0)-targetStart.at(0)) +
                              (sourceEnd.at(1)-targetStart.at(1))*(sourceEnd.at(1)-targetStart.at(1));

        if (alignedCheck < absTol)
            return 1;
        else if (reverseCheck < absTol)
            return -1;
        else
            return 0;
    }
    else if (sDim == 2 && tDim == 2)
    {
        gsInfo << "Accurate orientation check for two 3D is not implemented\n";
        return 1;
    }

    /*
    const T absTol = 1e-6; // another magic number

    gsMatrix<> startPoint;
    startPoint.resize(1,1);
    startPoint.at(0) = 0.;
    gsMatrix<> endPoint;
    endPoint.resize(1,1);
    endPoint.at(0) = 1.;

    gsMatrix<> sourceStart = sourceBoundary.eval(startPoint);
    gsMatrix<> sourceEnd = sourceBoundary.eval(endPoint);

    if (targetBoundary.dimensions().first == 2)
    {
        startPoint.conservativeResize(2,1);
        startPoint(1) = 0.;
        endPoint.conservativeResize(2,1);
        endPoint(1) = 0.;

    }
    gsMatrix<> targetStart = targetBoundary.eval(startPoint);
    gsMatrix<> targetEnd = targetBoundary.eval(endPoint);

    real_t alignedCheck = (sourceStart.at(0)-targetStart.at(0))*(sourceStart.at(0)-targetStart.at(0)) +
                          (sourceStart.at(1)-targetStart.at(1))*(sourceStart.at(1)-targetStart.at(1)) +
                          (sourceEnd.at(0)-targetEnd.at(0))*(sourceEnd.at(0)-targetEnd.at(0)) +
                          (sourceEnd.at(1)-targetEnd.at(1))*(sourceEnd.at(1)-targetEnd.at(1));

    real_t reverseCheck = (sourceStart.at(0)-targetEnd.at(0))*(sourceStart.at(0)-targetEnd.at(0)) +
                          (sourceStart.at(1)-targetEnd.at(1))*(sourceStart.at(1)-targetEnd.at(1)) +
                          (sourceEnd.at(0)-targetStart.at(0))*(sourceEnd.at(0)-targetStart.at(0)) +
                          (sourceEnd.at(1)-targetStart.at(1))*(sourceEnd.at(1)-targetStart.at(1));

    if (alignedCheck < absTol)
    {
        return 1;
    }
    else if (reverseCheck < absTol)
    {
        return -1;
    }
    else
    {
        gsInfo << "Doesn't look like matching boundaries!\n";
        gsInfo << "Source start\n" << sourceStart << std::endl << "Source end\n" << sourceEnd << std::endl;
        gsInfo << "Target start\n" << targetStart << std::endl << "Target end\n" << targetEnd << std::endl;
        return 0;
    }
    */
    return 0;
}

template <class T>
void gsElasticityAssembler<T>::addDirichletData(const gsMultiPatch<> & sourceGeometry,
                                                const gsMultiPatch<> & sourceSolution,
                                                int sourcePatch, const boxSide & sourceSide,
                                                int targetPatch, const boxSide & targetSide)
{
    int matchingMode = checkMatchingBoundaries(*(sourceGeometry.patch(sourcePatch).boundary(sourceSide)),
                                               *(m_pde_ptr->domain().patch(targetPatch).boundary(targetSide)));
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
void gsElasticityAssembler<T>::addDirichletData(const gsField<> & sourceField,
                                                int sourcePatch, const boxSide & sourceSide,
                                                int targetPatch, const boxSide & targetSide)
{
    int matchingMode = checkMatchingBoundaries(*(sourceField.patches().patch(sourcePatch).boundary(sourceSide)),
                                               *(m_pde_ptr->domain().patch(targetPatch).boundary(targetSide)));
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
void gsElasticityAssembler<T>::setDirichletDoFs(const gsMatrix<> & ddofs,
                                                int targetPatch,
                                                const boxSide & targetSide)
{
    const gsMatrix<unsigned> & localBoundaryIndices = (m_bases[0])[targetPatch].boundary(targetSide);

    for (index_t d = 0; d < m_dim; ++d)
    {
        gsMatrix<unsigned> systemIndices;
        const gsDofMapper & mapper = m_system.colMapper(d);
        mapper.localToGlobal(localBoundaryIndices,targetPatch,systemIndices);

        for (index_t i = 0; i != systemIndices.size(); ++i)
        {
            index_t dirichletIndex = mapper.global_to_bindex(systemIndices(i));
            m_ddof[0](dirichletIndex,0) = ddofs(i,d);
        }
    }
}


template <class T>
void gsElasticityAssembler<T>::clampPatchCorner(int patch,int corner)
{
    unsigned cornerIndex;
    if (corner == 0 || corner == 1)
    {
        const gsMatrix<unsigned> & boundaryIndices = m_bases[0].basis(patch).boundary(boundary::south);
        if (corner == 0)
            cornerIndex = boundaryIndices(0);
        else
            cornerIndex = boundaryIndices(boundaryIndices.rows()-1);
    }
    else if (corner == 2 || corner == 3)
    {
        const gsMatrix<unsigned> & boundaryIndices = m_bases[0].basis(patch).boundary(boundary::north);
        if (corner == 2)
            cornerIndex = boundaryIndices(0);
        else
            cornerIndex = boundaryIndices(boundaryIndices.rows()-1);
    }
    else
    {
        gsInfo << "Bad corner number!\n"
               << "0 - south-west\n"
               << "1 - south-east\n"
               << "2 - north-west\n"
               << "3 - north-east\n";
        return;
    }

    const T PP = 1e9; // magic number

    for (int i = 0; i < m_dim; ++i)
    {
        index_t globalIndex = m_system.colMapper(i).index(cornerIndex,patch);
        m_system.matrix().coeffRef(globalIndex,globalIndex) = PP;
        m_system.rhs().coeffRef(globalIndex) = 0;
    }
}

template <class T>
void gsElasticityAssembler<T>::deformGeometry(const gsMatrix<T> & solVector,
                                              gsMultiPatch<T> &result)
{
    result.clear();
    gsMatrix<T> newCoeffs;

    for (index_t p = 0; p < m_pde_ptr->domain().nPatches(); ++p)
    {
        index_t pNum = m_bases[0][p].size();
        newCoeffs.resize(pNum,m_dim);

        gsMatrix<T> & coeffs = m_pde_ptr->domain().patch(p).coefs();

        for (index_t d = 0; d < m_dim; ++d)
        {
            const gsDofMapper & mapper = m_system.colMapper(d);

            for (index_t i = 0; i < pNum; ++i )
            {
                if (mapper.is_free(i,p))
                {
                    newCoeffs(i,d) = coeffs(i,d) + solVector(mapper.index(i,p),0);
                }
                else
                {
                    newCoeffs(i,d) = coeffs(i,d) + m_ddof[0](mapper.bindex(i,p),0);
                }
            }
        }
        result.addPatch(m_bases[0][p].makeGeometry(give(newCoeffs)));
    }
}


template <class T>
void gsStressFunction<T>::eval_into(const gsMatrix< T > & u,gsMatrix< T > & result ) const
{
    gsMatrix<T> fullStress;
    result.clear();
    //gsInfo << m_displacement<<std::endl;
    m_assembler.computeStresses(m_displacement,u,m_patch,fullStress,true);
    result.resize(targetDim(),u.cols());

    switch(m_type)
    {
    case stress_type::normal:
    {
        result.block(0,0,targetDim(),u.cols()) = fullStress.block(0,0,targetDim(),u.cols());
        break;
    }
    case stress_type::shear:
    {
        result.block(0,0,targetDim(),u.cols()) = fullStress.block(domainDim(),0,targetDim(),u.cols());
        break;
    }
    case stress_type::von_mises:
    {
        result.block(0,0,1,u.cols()) = fullStress.block(fullStress.rows()-1,0,1,u.cols());
        break;
    }
    default:
        gsInfo << "Bad stress type\n";
        result.setZero();
        break;
    };
}

}// namespace gismo
