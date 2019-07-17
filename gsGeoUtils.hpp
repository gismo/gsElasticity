/** @file gsGeoUtils.hpp

    @brief Provides isogeometric meshing and modelling routines.

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s):
        A.Shamanskiy (2016 - ...., TU Kaiserslautern)
*/

#pragma once

#include <gsElasticity/gsGeoUtils.h>

#include <gsCore/gsField.h>
#include <gsCore/gsFuncData.h>
#include <gsCore/gsFunctionExpr.h>
#include <gsNurbs/gsBSpline.h>
#include <gsNurbs/gsTensorBSplineBasis.h>
#include <gsAssembler/gsQuadRule.h>
#include <gsAssembler/gsQuadrature.h>

#include <gsElasticity/gsElasticityAssembler.h>
#include <gsElasticity/gsElasticityFunctions.h>
#include <gsElasticity/gsWriteParaviewMultiPhysics.h>

namespace gismo
{

//-----------------------------------//
//--------- Mesh Analysis -----------//
//-----------------------------------//

template <class T>
void plotGeometry(gsMultiPatch<T> const & domain, std::string fileName, index_t numSamples)
{
    std::string fileNameOnly = fileName.substr(fileName.find_last_of("/\\")+1);
    gsParaviewCollection collectionMesh(fileName + "_mesh");
    gsParaviewCollection collectionJac(fileName + "_jac");
    index_t res;

    bool plotJac = true;
    if (numSamples == 0)
        plotJac = false;

    gsPiecewiseFunction<T> dets;
    for (size_t p = 0; p < domain.nPatches(); ++p)
        dets.addPiecePointer(new gsDetFunction<T>(domain,p));
    gsField<> detField(domain,dets,true);
    std::map<std::string,const gsField<> *> fields;
    fields["Jacobian"] = &detField;
    gsWriteParaviewMultiPhysics(fields,fileName,numSamples,true);

    for (size_t p = 0; p < domain.nPatches(); ++p)
    {
        collectionMesh.addPart(fileNameOnly,p,"_mesh.vtp"); 
        if (plotJac)
            collectionJac.addPart(fileNameOnly,p,".vts");
        else
            res = system(("rm " + fileName + std::to_string(p) + ".vts").c_str());
        GISMO_ENSURE(res == 0, "Problems with deleting files\n");
    }
    res = system(("rm " + fileName + ".pvd").c_str());
    GISMO_ENSURE(res == 0, "Problems with deleting files\n");
    (void)res;

    collectionMesh.save();
    if (plotJac)
        collectionJac.save();
}

template <class T>
void plotDeformation(const gsMultiPatch<T> & initDomain, const std::vector<gsMultiPatch<T> > & displacements,
                                             std::string fileName, index_t numSamplingPoints)
{
    gsInfo << "Plotting deformed configurations...\n";

    std::string fileNameOnly = fileName.substr(fileName.find_last_of("/\\")+1);
    gsParaviewCollection collectionMesh(fileName + "_mesh");
    gsParaviewCollection collectionJac(fileName + "_jac");
    index_t res;

    gsMultiPatch<T> configuration;
    for (size_t p = 0; p < initDomain.nPatches(); ++p)
        configuration.addPatch(initDomain.patch(p).clone());

    gsPiecewiseFunction<T> dets;
    for (size_t p = 0; p < configuration.nPatches(); ++p)
        dets.addPiecePointer(new gsDetFunction<T>(configuration,p));

    bool plotJac = true;
    if (numSamplingPoints == 0)
        plotJac = false;

    gsInfo << "Step: 0/" << displacements.size() << std::endl;

    gsField<T> detField(configuration,dets,true);
    std::map<std::string,const gsField<T> *> fields;
    fields["Jacobian"] = &detField;
    gsWriteParaviewMultiPhysics(fields,fileName+std::to_string(0),numSamplingPoints == 0 ? 1 : numSamplingPoints,true);

    for (size_t p = 0; p < configuration.nPatches(); ++p)
    {
        collectionMesh.addTimestep(fileNameOnly + std::to_string(0),p,0,"_mesh.vtp");
        if (plotJac)
            collectionJac.addTimestep(fileNameOnly + std::to_string(0),p,0,".vts");
        else
        {
            res = system(("rm " + fileName + std::to_string(0) + std::to_string(p) + ".vts").c_str());
            GISMO_ENSURE(res == 0, "Problems with deleting files\n");
        }
    }
    res = system(("rm " + fileName + std::to_string(0) + ".pvd").c_str());
    GISMO_ENSURE(res == 0, "Problems with deleting files\n");

    for (unsigned s = 0; s < displacements.size(); ++s)
    {
        gsInfo << "Step: " << s+1 << "/" << displacements.size() << std::endl;

        for (size_t p = 0; p < configuration.nPatches(); ++p)
        {
            configuration.patch(p).coefs() += displacements[s].patch(p).coefs();
            if (s > 0)
               configuration.patch(p).coefs() -= displacements[s-1].patch(p).coefs();
        }

        gsWriteParaviewMultiPhysics(fields,fileName+std::to_string(s+1),numSamplingPoints == 0 ? 1 : numSamplingPoints,true);
        for (size_t p = 0; p < configuration.nPatches(); ++p)
        {
            collectionMesh.addTimestep(fileNameOnly + std::to_string(s+1),p,s+1,"_mesh.vtp");

            if (plotJac)
                collectionJac.addTimestep(fileNameOnly + std::to_string(s+1),p,s+1,".vts");
            else
            {
                res = system(("rm " + fileName + std::to_string(s+1) + std::to_string(p) + ".vts").c_str());
                GISMO_ENSURE(res == 0, "Problems with deleting files\n");
            }
        }
        res = system(("rm " + fileName + std::to_string(s+1) + ".pvd").c_str());
        GISMO_ENSURE(res == 0, "Problems with deleting files\n");
    }

    collectionMesh.save();
    if (plotJac)
        collectionJac.save();

    (void)res;
}


template <class T>
index_t checkGeometry(gsMultiPatch<T> const & domain)
{
    index_t corruptedPatch = -1;
    gsMapData<T> md;
    md.flags = NEED_DERIV;
    gsMatrix<T> points;

    for (size_t p = 0; p < domain.nPatches(); ++p)
    {
        gsVector<index_t> numNodes(domain.dim());
        for (short_t i = 0; i < domain.dim(); ++i)
            numNodes.at(i) = domain.basis(p).degree(i)+1;
        gsQuadRule<T> quRule = gsQuadrature::get<T>(1,numNodes);

        typename gsBasis<T>::domainIter domIt = domain.basis(p).makeDomainIterator(boundary::none);
        for (; domIt->good(); domIt->next())
        {
            genSamplingPoints(domIt->lowerCorner(),domIt->upperCorner(),quRule,points);
            md.points = points;
            domain.patch(p).computeMap(md);
            for (index_t q = 0; q < points.cols(); ++q)
                if (md.jacobian(q).determinant() <= 0)
                    return p;
        }
    }
    return corruptedPatch;
}

template <class T>
T geometryJacRatio(gsMultiPatch<T> const & domain)
{
    std::vector<T> maxs, mins;
    gsMapData<T> md;
    md.flags = NEED_DERIV;
    gsMatrix<T> points;

    for (size_t p = 0; p < domain.nPatches(); ++p)
    {
        gsVector<index_t> numNodes(domain.dim());
        for (short_t i = 0; i < domain.dim(); ++i)
            numNodes.at(i) = domain.basis(p).degree(i)+1;
        gsQuadRule<T> quRule = gsQuadrature::get<T>(1,numNodes);

        typename gsBasis<T>::domainIter domIt = domain.basis(p).makeDomainIterator(boundary::none);
        for (; domIt->good(); domIt->next())
        {
            genSamplingPoints(domIt->lowerCorner(),domIt->upperCorner(),quRule,points);
            md.points = points;
            domain.patch(p).computeMap(md);

            T min = md.jacobian(0).determinant();
            T max = min;
            for (index_t q = 1; q < points.cols(); ++q)
            {
                T jac = md.jacobian(q).determinant();
                if (jac > max)
                    max = jac;
                if (jac < min)
                    min = jac;
            }

            maxs.push_back(max);
            mins.push_back(min);
        }
    }

    return *(std::min_element(mins.begin(),mins.end())) / *(std::max_element(maxs.begin(),maxs.end()));
}

//--------------------------------------------------------------------//
//------------------------ Modelling ---------------------------------//
//--------------------------------------------------------------------//

template<class T>
typename gsGeometry<T>::uPtr simplifyCurve(gsGeometry<T> const & curve,
                                          index_t additionalPoints, index_t degree,
                                          index_t numSamples)
{
    GISMO_ENSURE(curve.domainDim() == 1 ,"That's not a curve.\n");
    index_t deg = degree == 0 ? curve.degree(0) : degree;
    index_t num = deg + 1 + additionalPoints;

    gsKnotVector<T> knots(0.0,1.0, num - deg - 1, deg + 1);
    gsKnotVector<T> knotVector(knots);
    gsBSplineBasis<T> basis(knotVector);

    gsMatrix<T> params(1,numSamples);
    for (index_t p = 0; p < numSamples; ++p)
        params.at(p) = 1.*p/(numSamples-1);

    gsMatrix<T> curveValues;
    curve.eval_into(params,curveValues);

    gsMatrix<T> lens(numSamples,1);
    lens.at(0) = 0.;
    for (index_t p = 1; p < numSamples; ++p)
        lens.at(p) = lens.at(p-1) + distance(curveValues,p,curveValues,p-1,true);

    gsMatrix<T> lenParams(1,numSamples);
    for (index_t p = 0; p < numSamples; ++p)
        lenParams.at(p) = lens.at(p)/lens.at(numSamples-1);

    typename gsGeometry<T>::uPtr simpleCurve = fittingDirichlet(lenParams,curveValues,basis);

    simpleCurve->eval_into(params,curveValues);
    for (index_t p = 1; p < numSamples; ++p)
        lens.at(p) = lens.at(p-1) + distance(curveValues,p,curveValues,p-1,true);

    for (index_t p = 0; p < numSamples; ++p)
        lenParams.at(p) = lens.at(p)/lens.at(numSamples-1);

    return fittingDirichlet(lenParams,curveValues,curve.basis());
}

template<class T>
T curveDistance(gsGeometry<T> const & curveA,
                gsGeometry<T> const & curveB,
                index_t numSamples)
{
    gsMatrix<T> params(1,numSamples);
    for (index_t p = 0; p < numSamples; ++p)
        params.at(p) = 1.*p/(numSamples-1);

    gsMatrix<T> pointsA, pointsB;
    curveA.eval_into(params,pointsA);
    curveB.eval_into(params,pointsB);

    T dist = 0.;
    for (int p = 0; p < numSamples; ++p)
        dist += pow(distance(pointsA,p,pointsB,p,true),2);

    return sqrt(dist/numSamples);
}

template <class T>
typename gsGeometry<T>::uPtr fittingDirichlet(gsMatrix<T> const & params,
                                              gsMatrix<T> const & points,
                                              gsBasis<T> const & basis)
{
    index_t numSamples = params.cols();
    unsigned num = basis.size();
    index_t dim = points.rows();

    gsSparseMatrix<T> A(num,num);
    gsMatrix<T> b(num,dim);
    b.setZero();

    gsMatrix<T> basisValues;
    gsMatrix<unsigned> activeBasis;

    basis.eval_into(params,basisValues);
    basis.active_into(params,activeBasis);

    for (index_t p = 1; p < numSamples; ++p)
    {
        index_t numActive = activeBasis.rows();
        for (index_t i = 0; i < numActive; ++i)
        {
            if (activeBasis(i,p) == 0)
            {
                for (index_t j = 0; j < numActive; ++j)
                {
                    if (activeBasis(j,p) != 0 && activeBasis(j,p) != num-1 )
                        b.row(activeBasis(j,p)-1) -= basisValues(i,p) * basisValues(j,p) *
                                                     points.col(0).transpose() * (params.at(p)-params.at(p-1));
                }
            }
            else if (activeBasis(i,p) == num-1)
            {
                for (index_t j = 0; j < numActive; ++j)
                {
                    if (activeBasis(j,p) != 0 && activeBasis(j,p) != num-1 )
                        b.row(activeBasis(j,p)-1) -= basisValues(i,p) * basisValues(j,p) *
                                                     points.col(numSamples-1).transpose() * (params.at(p)-params.at(p-1));
                }
            }
            else
            {
                b.row(activeBasis(i,p)-1) += basisValues(i,p) * points.col(p).transpose() * (params.at(p)-params.at(p-1));
                for (index_t j = 0; j < numActive; ++j)
                    if (activeBasis(j,p) != 0 && activeBasis(j,p) != num-1 )
                        A(activeBasis(i,p)-1, activeBasis(j,p)-1) += basisValues(i,p) * basisValues(j,p) * (params.at(p)-params.at(p-1));
            }
        }
    }

    A.makeCompressed();
    typename gsSparseSolver<T>::CGDiagonal solver(A);
    gsMatrix<T> x = solver.solve(b);

    gsMatrix<T> coefs(num,dim);
    coefs.row(0) = points.col(0).transpose();
    coefs.block(1,0,num-2,dim) = x.block(0,0,num-2,dim);
    coefs.row(num-1) = points.col(numSamples-1).transpose();

    return basis.makeGeometry(give(coefs));
}

template<class T>
typename gsGeometry<T>::uPtr genPatchInterpolation(gsGeometry<T> const & A, gsGeometry<T> const & B,
                                                   index_t deg, index_t num, bool xiDir)
{
    GISMO_ENSURE(A.parDim() == B.parDim(), "Geometries are incompatible: different parametric dimensions: " +
                                           std::to_string(A.parDim()) + " and " + std::to_string(B.parDim()) + "\n");
    short_t pDim = A.parDim();
    GISMO_ASSERT(pDim == 1 || pDim ==2, "Can only interpolate between curves or surfaces. Given geometries have parametric dimension " +
                                        std::to_string(pDim) + "\n");
    for (index_t d = 0; d < pDim; ++d)
        GISMO_ENSURE(A.degree(d) == B.degree(d), "Geometries are incompatible: different splines degrees in dimension" +
                                                 std::to_string(d) + ": " + std::to_string(A.degree(d)) +
                                                 " and " + std::to_string(B.degree(d)) + "\n");

    GISMO_ENSURE(A.targetDim() == B.targetDim(), "Geometries are incompatible: different physical dimensions: " +
                                                 std::to_string(A.targetDim()) + " and " + std::to_string(B.targetDim()) + "\n");
    short_t tDim = A.targetDim();
    GISMO_ASSERT(A.coefsSize() == B.coefsSize(), "Geometries are incompatible: different number of control points: " +
                                                 std::to_string(A.coefsSize()) + " and " + std::to_string(B.coefsSize()) + "\n");
    index_t baseNum = A.coefsSize();

    gsKnotVector<T> newKnots(0.0,1.0, num - deg - 1, deg + 1);

    gsMultiPatch<> temp;
    temp.addPatch(A.clone());

    if (pDim == 1)
    {

        gsMatrix<T> coefs(baseNum*num,tDim);
        T part = 1./(num-deg);
        T pos = 0.;
        for (index_t i = 0; i < num; ++i)
        {
            if (i < deg)
                pos += part/deg*i;
            else if (deg <= i && i < num-deg )
                pos += part;
            else
                pos += part/deg*(num-i);

            for (index_t j = 0; j < baseNum; ++j)
            {
                if (xiDir)
                    coefs.row(j*num+i) = combine(A.coefs(),B.coefs(),pos,j,j).row(0);
                else
                    coefs.row(i*baseNum+j) = combine(A.coefs(),B.coefs(),pos,j,j).row(0);
            }
        }

        if (xiDir)
        {
            gsTensorBSplineBasis<2,T> basis(newKnots,
                                            (static_cast<gsBSpline<T> &>(temp.patch(0))).knots());
            return basis.makeGeometry(give(coefs));
        }
        else
        {
            gsTensorBSplineBasis<2,T> basis((static_cast<gsBSpline<T> &>(temp.patch(0))).knots(),
                                            newKnots);
            return basis.makeGeometry(give(coefs));
        }
    }
    else
    {
        gsTensorBSplineBasis<3,T> basis((static_cast<gsTensorBSpline<2,T> &>(temp.patch(0))).knots(0),
                                        (static_cast<gsTensorBSpline<2,T> &>(temp.patch(0))).knots(1),
                                        newKnots);

        gsMatrix<T> coefs(baseNum*num,tDim);
        T part = 1./(num-deg);
        T pos = 0.;
        for (index_t i = 0; i < num; ++i)
        {
            if (i < deg)
                pos += part/deg*i;
            else if (deg <= i && i < num-deg )
                pos += part;
            else
                pos += part/deg*(num-i);

            for (index_t j = 0; j < baseNum; ++j)
                coefs.row(i*baseNum+j) = combine(A.coefs(),B.coefs(),pos,j,j).row(0);
        }

        return basis.makeGeometry(give(coefs));
    }
}

template<class T>
typename gsGeometry<T>::uPtr genPatchScaling(gsGeometry<T> const & boundary,
                                             index_t deg, index_t num,
                                             T scaling, gsVector<T> const & center)
{
    typename gsGeometry<T>::uPtr scaledBoundary = boundary.clone();
    scaledBoundary->translate(-1*center);
    scaledBoundary->scale(scaling);
    scaledBoundary->translate(center);
    return genPatchInterpolation(boundary,*scaledBoundary,deg,num);
}

template<class T>
typename gsGeometry<T>::uPtr genLine(index_t deg, index_t num,
                                     gsMatrix<T> const & A, gsMatrix<T> const & B,
                                     index_t iA, index_t iB)
{
    GISMO_ENSURE(num - deg - 1 >= 0,"Too few DoFs\n");
    GISMO_ENSURE(A.cols() == B.cols(),"Points have different dimensions\n");
    gsKnotVector<T> knots(0.0,1.0, num - deg - 1, deg + 1);
    gsBSplineBasis<T> basis(knots);

    index_t dim = A.cols();
    gsMatrix<T> coefs(num,dim);

    T part = 1./(num-deg);
    T pos = 0.;
    for (index_t i = 0; i < num; ++i)
    {
        if (i < deg)
            pos += part/deg*i;
        else if (deg <= i && i < num-deg )
            pos += part;
        else
            pos += part/deg*(num-i);

        coefs.row(i) = combine(A,B,pos,iA,iB).row(0);
    }

    return basis.makeGeometry(give(coefs));
}

template<class T>
typename gsGeometry<T>::uPtr genCircle(index_t deg, index_t num,
                                       T radius, T x0, T y0,
                                       T angle0, T arcAngle)
{
    GISMO_ENSURE(num - deg - 1 >= 0,"Too few DoFs\n");
    gsKnotVector<T> knots(0.0,1.0, num - deg - 1, deg + 1);
    gsBSplineBasis<T> basis(knots);
    return genCircle(basis,radius,x0,y0,angle0,arcAngle);
}

template<class T>
typename gsGeometry<T>::uPtr genCircle(gsBasis<T> & basis,
                                       T radius,
                                       T x0, T y0,
                                       T angle0, T arcAngle)
{
    index_t numPoints = 1000;
    gsMatrix<> params(1,numPoints);
    gsMatrix<> points(2,numPoints);
    for (index_t i = 0; i < numPoints; ++i)
    {
        params(0,i) = 1.*i/(numPoints-1);
        points(0,i) = x0 + radius*cos(angle0 + i*arcAngle/(numPoints-1));
        points(1,i) = y0 + radius*sin(angle0 + i*arcAngle/(numPoints-1));
    }

    return fittingDirichlet(params,points,basis);
}

template<class T>
typename gsGeometry<T>::uPtr genQuad(index_t xiDeg, index_t xiNum, index_t etaDeg, index_t etaNum,
                                     gsMatrix<T> const & A, gsMatrix<T> const & B,
                                     gsMatrix<T> const & C, gsMatrix<T> const & D,
                                     index_t iA, index_t iB, index_t iC, index_t iD)
{
    typename gsGeometry<T>::uPtr sideAB = genLine(xiDeg,xiNum,A,B,iA,iB);
    typename gsGeometry<T>::uPtr sideCD = genLine(xiDeg,xiNum,C,D,iC,iD);

    return genPatchInterpolation(*sideAB,*sideCD,etaDeg,etaNum);
}

template<class T>
typename gsGeometry<T>::uPtr genSphere(index_t xiDeg, index_t xiNum, index_t etaDeg, index_t etaNum,
                                       T xi0, T xi1, T eta0, T eta1)
{
    GISMO_ENSURE(xiNum - xiDeg - 1 >= 0,"Too few DoFs\n");
    GISMO_ENSURE(etaNum - etaDeg - 1 >= 0,"Too few DoFs\n");

    gsKnotVector<T> xiKnots(0.,1.,xiNum - xiDeg - 1, xiDeg + 1);
    gsKnotVector<T> etaKnots(0.,1.,etaNum - etaDeg - 1, etaDeg + 1);
    return genSphere(xiKnots,etaKnots,xi0,xi1,eta0,eta1);

}

template<class T>
typename gsGeometry<T>::uPtr genSphere(gsKnotVector<T> & xiKnots, gsKnotVector<T> & etaKnots,
                                       T xi0, T xi1, T eta0, T eta1)
{
    gsBSplineBasis<T> xiBasis(xiKnots);
    typename gsGeometry<T>::uPtr xiCircle = genCircle(xiBasis,1.,0.,0.,xi0,xi1-xi0);
    gsBSplineBasis<T> etaBasis(etaKnots);
    typename gsGeometry<T>::uPtr etaCircle = genCircle(xiBasis,1.,0.,0.,xi0,xi1-xi0);

    gsTensorBSplineBasis<2,T> basis(xiKnots,etaKnots);

    index_t xiNum = xiCircle->coefsSize();
    index_t etaNum = etaCircle->coefsSize();
    gsMatrix<T> coefs(xiNum*etaNum,3);
    for (index_t i = 0; i < xiNum; ++i)
    {
        for (index_t j = 0; j < etaNum; ++j)
        {
            coefs(j*xiNum+i,0) = xiCircle->coef(i,0)*etaCircle->coef(j,1);
            coefs(j*xiNum+i,1) = xiCircle->coef(i,1)*etaCircle->coef(j,1);
            coefs(j*xiNum+i,2) = etaCircle->coef(j,0);
        }
    }

    return basis.makeGeometry(give(coefs));
}

template<class T>
typename gsGeometry<T>::uPtr genCylinder(gsGeometry<T> const & base,
                                         index_t deg, index_t num, T height)
{
    GISMO_ENSURE(num - deg - 1 >= 0,"Too few DoFs\n");

    typename gsGeometry<T>::uPtr botBoundary = base.clone();
    typename gsGeometry<T>::uPtr topBoundary = base.clone();
    if (base.targetDim() < 3)
    {
        botBoundary->embed3d();
        topBoundary->embed3d();
    }

    topBoundary->translate(gsVector<T,3>::vec(0.,0.,height));
    return genPatchInterpolation(*botBoundary,*topBoundary,deg,num);
}

template<class T>
typename gsGeometry<T>::uPtr genScrew(gsGeometry<T> const & base,
                                      index_t deg, index_t num,
                                      T height, T pitch, T x0, T y0)
{
    GISMO_ENSURE(num - deg - 1 >= 0,"Too few DoFs\n");

    short_t pDim = base.parDim();
    GISMO_ASSERT(pDim == 1 || pDim ==2,"Wrong geometry type\n");

    gsKnotVector<> zKnots(0.0,1.0, num - deg - 1, deg + 1);
    gsBSplineBasis<> zBasis(zKnots);

    index_t numBase = base.coefsSize();

    index_t numPoints = 10000;
    gsMatrix<> params(1,numPoints);
    gsMatrix<> points(3,numPoints);
    for (index_t i = 0; i < numPoints; ++i)
    {
        params(0,i) = 1.*i/(numPoints-1);
        points(0,i) = cos(params(0,i)*2*EIGEN_PI*pitch/360);
        points(1,i) = sin(params(0,i)*2*EIGEN_PI*pitch/360);
        points(2,i) = height*params(0,i);
    }

    typename gsGeometry<T>::uPtr helix = fittingDirichlet(params,points,zBasis);

    gsMatrix<T> coefs(numBase*num,3);
    gsMultiPatch<T> temp;
    temp.addPatch(base.clone());

    T oldAngle = 0.;
    T oldRadius = 1.;
    for (index_t i = 0; i < num; ++i)
    {
        T x = helix->coef(i,0);
        T y = helix->coef(i,1);
        T radius = sqrt(x*x+y*y);
        T angle = atan2(y,x);
        temp.patch(0).translate(gsVector<T,2>::vec(-1*x0,-1*y0));
        temp.patch(0).rotate(angle-oldAngle);
        temp.patch(0).scale(radius/oldRadius);
        temp.patch(0).translate(gsVector<T,2>::vec(x0,y0));

        for (index_t j = 0; j < numBase; j++)
        {
            coefs(i*numBase + j,0) = temp.patch(0).coef(j,0);
            coefs(i*numBase + j,1) = temp.patch(0).coef(j,1);
            coefs(i*numBase + j,2) = helix->coef(i,2);
        }
        oldRadius = radius;
        oldAngle = angle;
    }

    if (pDim == 1)
    {
        gsTensorBSplineBasis<2,T> basis(static_cast<gsBSpline<T> &>(temp.patch(0)).knots(0),
                                        zKnots);
        return basis.makeGeometry(give(coefs));
    }
    else
    {
        gsTensorBSplineBasis<3,T> basis(static_cast<gsTensorBSpline<2,T> &>(temp.patch(0)).knots(0),
                                        static_cast<gsTensorBSpline<2,T> &>(temp.patch(0)).knots(1),
                                        zKnots);
        return basis.makeGeometry(give(coefs));
    }
}

//------------------------------------------------------//
//----------------- Auxiliary functions ----------------//
//------------------------------------------------------//

template<class T>
gsMatrix<T> combine(gsMatrix<T> const & A, gsMatrix<T> const & B, T x,
                    index_t iA, index_t iB, bool cols)
{
    if (cols)
    {
        GISMO_ENSURE(A.rows() == B.rows(),"Points have different dimensions\n");
        index_t dim = A.rows();
        gsMatrix<T> combination(dim,1);
        combination.col(0) = (1-x)*A.col(iA) + x*B.col(iB);
        return combination;
    }
    else
    {
        GISMO_ENSURE(A.cols() == B.cols(),"Points have different dimensions\n");
        index_t dim = A.cols();
        gsMatrix<T> combination(1,dim);
        combination.row(0) = (1-x)*A.row(iA) + x*B.row(iB);
        return combination;
    }
}

template <class T>
T distance(gsMatrix<T> const & A, index_t i, gsMatrix<T> const & B, index_t j, bool cols)
{
    T dist = 0.;

    if (cols)
    {
        GISMO_ENSURE(A.rows() == B.rows(),"Wrong matrix size\n");
        for (index_t d = 0; d < A.rows(); ++d)
            dist = sqrt(pow(dist,2)+pow(A(d,i)-B(d,j),2));
    }
    else
    {
        GISMO_ENSURE(A.cols() == B.cols(),"Wrong matrix size\n");
        for (index_t d = 0; d < A.cols(); ++d)
            dist = sqrt(pow(dist,2)+pow(A(i,d)-B(j,d),2));
    }

    return dist;
}

} // namespace ends
