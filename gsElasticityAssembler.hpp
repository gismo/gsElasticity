/** @file gsElasticityAssembler.hpp

    @brief Provides linear and nonlinear elasticity systems for 2D plain strain and 3D continua.

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s):
        O. Weeger    (2012 - 2015, TU Kaiserslautern),
        A.Shamanskiy (2016 - ...., TU Kaiserslautern)
*/

#pragma once

#include <gsElasticity/gsElasticityAssembler.h>

#include <gsCore/gsMultiBasis.h>
#include <gsCore/gsField.h>
#include <gsPde/gsPoissonPde.h>

// Element visitors
#include <gsElasticity/gsVisitorLinearElasticity.h>
//#include <gsElasticity/gsVisitorNonLinElasticity.h>
#include <gsElasticity/gsVisitorElasticityNeumann.h>
//#include <gsElasticity/gsVisitorElasticityPressure.h>
//#include <gsElasticity/gsVisitorMassElasticity.h>


namespace gismo
{

template<class T>
gsElasticityAssembler<T>::gsElasticityAssembler(gsMultiPatch<T> const & patches,
                                                gsMultiBasis<T> const & basis,
                                                gsBoundaryConditions<T> const & bconditions,
                                                const gsFunction<T> & body_force)
{
    // Always concieved as a meaningful class, now gsPde is just a container for
    // the domain, boundary conditions and the right-hand side;
    // any derived class can surve this purpuse, for example gsPoissonPde;
    // TUDO: change/remove gsPde from gsAssembler logic
    gsPiecewiseFunction<T> rightHandSides;
    rightHandSides.addPiece(body_force);
    typename gsPde<T>::Ptr pde( new gsPoissonPde<T>(patches,bconditions,rightHandSides) );
    // gsAssembler<>::initialize requires a vector of bases, one for each unknown;
    // different bases are used to compute Dirichlet DoFs;
    // but always the first basis is used for the assembly;
    // TODO: change gsAssembler logic
    m_dim = body_force.targetDim();
    for (int d = 0; d < m_dim; ++d)
        m_bases.push_back(basis);

    Base::initialize(pde, m_bases, defaultOptions());
}

template <class T>
gsOptionList gsElasticityAssembler<T>::defaultOptions()
{
    gsOptionList opt = Base::defaultOptions();
    opt.addReal("YoungsModulus","Youngs modulus of the material",200e9);
    opt.addReal("PoissonsRatio","Poisson's ratio of the material",0.33);
    opt.addReal("Density","Density of the material",8e3);
    opt.addReal("TimeFactor","Time factor for the time-dependent forces",1.);
    opt.addInt("MaterialLaw","Material law: 0 for St. Venant-Kirchhof, 1 for Neo-Hooke",0);
    return opt;
}

template <class T>
void gsElasticityAssembler<T>::refresh()
{
    GISMO_ASSERT(m_dim == m_pde_ptr->domain().parDim(), "The RHS dimension and the domain dimension don't match!");
    GISMO_ASSERT(m_dim == 2 || m_dim == 3, "Only two- and three-dimenstion domains are supported!");

    std::vector<gsDofMapper> m_dofMappers(m_dim);
    for (index_t d = 0; d < m_dim; d++)
        m_bases[d].getMapper((dirichlet::strategy)m_options.getInt("DirichletStrategy"),
                             iFace::glue,m_pde_ptr->bc(),m_dofMappers[d],d,true);

    gsVector<unsigned> dims;
    dims.setOnes(m_dim);
    m_system = gsSparseSystem<T>(m_dofMappers, dims);
}

template<class T>
void gsElasticityAssembler<T>::assemble()
{
    if ( this->numDofs() == 0 )
    {
        gsWarn << " No internal DOFs. Computed Dirichlet boundary only.\n";
        return;
    }

    m_system.reserve(m_bases[0], m_options, 1);

    for (int d = 0; d < m_dim; ++d)
        Base::computeDirichletDofs(d);

    // Compute volumetric integrals and write to the global linear system
    Base::template push<gsVisitorLinearElasticity<T> >();
    // Compute surface integrals and write to the global rhs vector
    Base::template push<gsVisitorElasticityNeumann<T> >(m_pde_ptr->bc().neumannSides());

    m_system.matrix().makeCompressed();
}

template <class T>
void gsElasticityAssembler<T>::constructCauchyStresses(const gsMultiPatch<T> & displacement,
                                                       gsPiecewiseFunction<T> & result,
                                                       stress_type::type type) const
{
    result.clear();
    if (type == stress_type::all_2D)
        GISMO_ASSERT(m_dim == 2, "Invalid stress type for a 2D problem");
    if (type == stress_type::normal_3D || type == stress_type::shear_3D)
        GISMO_ASSERT(m_dim == 3, "Invalid stress type for a 3D problem");

    T E = m_options.getReal("YoungsModulus");
    T pr = m_options.getReal("PoissonsRatio");
    T lambda = E * pr / ( ( 1. + pr ) * ( 1. - 2. * pr ) );
    T mu     = E / ( 2. * ( 1. + pr ) );

    for (index_t p = 0; p < m_pde_ptr->domain().nPatches(); ++p )
        result.addPiecePointer(new gsCauchyStressFunction<T>(displacement,p,type,lambda,mu));
}

template<class T>
void gsElasticityAssembler<T>::assemble(const gsMultiPatch<T> & deformed)
{
    GISMO_NO_IMPLEMENTATION
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
*/
}


template<class T>
void  gsElasticityAssembler<T>::reComputeDirichletDofs(gsMultiPatch<T> &deformed)
{
    GISMO_NO_IMPLEMENTATION
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
    GISMO_NO_IMPLEMENTATION
    /*GISMO_ASSERT(this->numDofs() == m_system.rhs().rows(), "Something went wrong, assemble() not called?");

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
    }*/
}

template<class T>
void  gsElasticityAssembler<T>::updateSolution(const gsMatrix<T>& solVector,
                                               gsMultiPatch<T>& result) const
{
    GISMO_NO_IMPLEMENTATION
    /*GISMO_ASSERT(this->numDofs() == m_system.rhs().rows(), "Something went wrong, assemble() not called?");

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
    }*/
}

template<class T>
void gsElasticityAssembler<T>::assembleMass()
{
    GISMO_NO_IMPLEMENTATION
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
void gsElasticityAssembler<T>::deformGeometry(const gsMatrix<T> & solVector,
                                              gsMultiPatch<T> &result)
{
    GISMO_NO_IMPLEMENTATION
    /*result.clear();
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
    }*/
}

}// namespace gismo ends
