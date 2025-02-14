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
#include <gsNurbs/gsTensorNurbs.h>
#include <gsAssembler/gsQuadRule.h>
#include <gsAssembler/gsQuadrature.h>
#include <gsModeling/gsCoonsPatch.h>

#include <gsElasticity/gsElasticityAssembler.h>
#include <gsElasticity/gsElasticityFunctions.h>
#include <gsElasticity/gsWriteParaviewMultiPhysics.h>
#include <gsUtils/gsMesh/gsMesh.h>
#include <gsIO/gsWriteParaview.h>

namespace gismo
{

//-----------------------------------//
//--------- Mesh Analysis -----------//
//-----------------------------------//

template <class T>
void plotGeometry(gsMultiPatch<T> const & domain, std::string fileName, index_t numSamples)
{
    std::string fileNameOnly = gsFileManager::getFilename(fileName);
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
        collectionMesh.addPart(fileNameOnly + "_mesh.vtp",0,"Mesh",p);
        if (plotJac)
            collectionJac.addPart(fileNameOnly + ".vts",0,"Solution",p);
        else
            res = system(("rm " + fileName + util::to_string(p) + ".vts").c_str());
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
void plotGeometry(const gsMultiPatch<T> & domain, std::string const & fileName,
                  gsParaviewCollection & collection, index_t step)
{
    for (size_t p = 0; p < domain.nPatches(); ++p)
    {
        gsMesh<T> mesh(domain.basis(p),8);
        domain.patch(p).evaluateMesh(mesh);
        std::string patchFileName = fileName + util::to_string(step) + "_" + util::to_string(p);
        gsWriteParaview(mesh,patchFileName,false);
        collection.addPart(gsFileManager::getFilename(patchFileName),step,"Solution",p);
    }
}

template <class T>
void plotDeformation(const gsMultiPatch<T> & initDomain, const std::vector<gsMultiPatch<T> > & displacements,
                                             std::string fileName, index_t numSamplingPoints)
{
    gsInfo << "Plotting deformed configurations...\n";

    std::string fileNameOnly = gsFileManager::getFilename(fileName);
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

    gsProgressBar bar;
    bar.display(0, displacements.size());

    gsField<T> detField(configuration,dets,true);
    std::map<std::string,const gsField<T> *> fields;
    fields["Jacobian"] = &detField;
    gsWriteParaviewMultiPhysics(fields,fileName+util::to_string(0),numSamplingPoints == 0 ? 1 : numSamplingPoints,true);

    for (size_t p = 0; p < configuration.nPatches(); ++p)
    {
        collectionMesh.addPart(fileNameOnly + util::to_string(0) + util::to_string(p) + "_mesh.vtp",0,"Mesh",p);
        if (plotJac)
            collectionJac.addPart(fileNameOnly + util::to_string(0) + util::to_string(p) + ".vts",0,"Solution",p);
        else
        {
            res = system(("rm " + fileName + util::to_string(0) + util::to_string(p) + ".vts").c_str());
            GISMO_ENSURE(res == 0, "Problems with deleting files\n");
        }
    }
    res = system(("rm " + fileName + util::to_string(0) + ".pvd").c_str());
    GISMO_ENSURE(res == 0, "Problems with deleting files\n");

    for (size_t s = 0; s < displacements.size(); ++s)
    {

        bar.display(s+1, displacements.size());

        for (size_t p = 0; p < configuration.nPatches(); ++p)
        {
            configuration.patch(p).coefs() += displacements[s].patch(p).coefs();
            if (s > 0)
               configuration.patch(p).coefs() -= displacements[s-1].patch(p).coefs();
        }

        gsWriteParaviewMultiPhysics(fields,fileName+util::to_string(s+1),numSamplingPoints == 0 ? 1 : numSamplingPoints,true);
        for (size_t p = 0; p < configuration.nPatches(); ++p)
        {
            collectionMesh.addPart(fileNameOnly + util::to_string(s+1) + util::to_string(p) +  "_mesh.vtp",s+1,"Mesh",p);
            if (plotJac)
                collectionJac.addPart(fileNameOnly + util::to_string(s+1) + util::to_string(p) + ".vts",s+1,"Solution",p);
            else
            {
                res = system(("rm " + fileName + util::to_string(s+1) + util::to_string(p) + ".vts").c_str());
                GISMO_ENSURE(res == 0, "Problems with deleting files\n");
            }
        }
        res = system(("rm " + fileName + util::to_string(s+1) + ".pvd").c_str());
        GISMO_ENSURE(res == 0, "Problems with deleting files\n");
    }

    collectionMesh.save();
    if (plotJac)
        collectionJac.save();

    (void)res;
}


template <class T>
void plotDeformation(const gsMultiPatch<T> & initDomain, const gsMultiPatch<T> & displacement,
                     std::string const & fileName, gsParaviewCollection & collection, index_t step)
{
    GISMO_ENSURE(initDomain.nPatches() == displacement.nPatches(), "Wrong number of patches! Geometry has " +
                 util::to_string(initDomain.nPatches()) + " patches. Displacement has " + util::to_string(displacement.nPatches()) + " patches.");

    gsMultiPatch<T> configuration;
    for (size_t p = 0; p < initDomain.nPatches(); ++p)
    {
        configuration.addPatch(initDomain.patch(p).clone());
        configuration.patch(p).coefs() += displacement.patch(p).coefs();
    }

    for (size_t p = 0; p < configuration.nPatches(); ++p)
    {
        gsMesh<T> mesh(configuration.basis(p),8);
        configuration.patch(p).evaluateMesh(mesh);
        std::string patchFileName = fileName + util::to_string(step) + "_" + util::to_string(p);
        gsWriteParaview(mesh,patchFileName,false);
        collection.addPart(gsFileManager::getFilename(patchFileName),step,"Solution",p);
    }
}

template <class T>
index_t checkGeometry(gsMultiPatch<T> const & domain)
{
    index_t corruptedPatch = -1;
    bool continueIt = true;
    for (size_t p = 0; p < domain.nPatches() && continueIt; ++p)
    {
#pragma omp parallel
        {
            gsMapData<T> md;
            md.flags = NEED_DERIV;
            gsMatrix<T> points;

            gsVector<index_t> numNodes(domain.dim());
            for (short_t i = 0; i < domain.dim(); ++i)
                numNodes.at(i) = domain.basis(p).degree(i)+1;
            gsQuadRule<T> quRule = gsQuadrature::get<T>(gsQuadrature::GaussLegendre,numNodes);

            // Initialize domain element iterator
            typename gsBasis<T>::domainIter domIt    = domain.basis(p).domain()->beginBdr(boundary::none);
            typename gsBasis<T>::domainIter domItEnd = domain.basis(p).domain()->endBdr(boundary::none);

#           ifdef _OPENMP
            const int tid = omp_get_thread_num();
            const int nt  = omp_get_num_threads();
            domIt += tid;
            for ( ; domIt<domItEnd; domIt+=nt )
#           else
            for (; domIt<domItEnd && continueIt; ++domIt )
#           endif
            {
                genSamplingPoints(domIt.lowerCorner(),domIt.upperCorner(),quRule,points);
                md.points = points;
                domain.patch(p).computeMap(md);
                for (index_t q = 0; q < points.cols() && continueIt; ++q)

                    if (md.jacobian(q).determinant() <= 0)
                    {
                        gsInfo << "Bad patch: " << p << "\nBad point:\n" << points.col(q) << "\nDet: " << md.jacobian(q).determinant() << std::endl;
                        corruptedPatch = p;
                        continueIt = false;
                    }
            }
        }
    }
    return corruptedPatch;
}


template <class T>
index_t checkDisplacement(gsMultiPatch<T> const & domain, gsMultiPatch<T> const & displacement)
{
    index_t corruptedPatch = -1;
    bool continueIt = true;
    for (size_t p = 0; p < domain.nPatches() && continueIt; ++p)
    {
#pragma omp parallel
        {
            gsMapData<T> mdG, mdU;
            mdG.flags = NEED_DERIV;
            mdU.flags = NEED_DERIV;
            gsMatrix<T> points;

            gsVector<index_t> numNodes(domain.dim());
            for (short_t i = 0; i < domain.dim(); ++i)
                numNodes.at(i) = displacement.basis(p).degree(i)+1;
            gsQuadRule<T> quRule = gsQuadrature::get<T>(gsQuadrature::GaussLegendre,numNodes);

            // Initialize domain element iterator
            typename gsBasis<T>::domainIter domIt    = displacement.basis(p).domain()->beginBdr(boundary::none);
            typename gsBasis<T>::domainIter domItEnd = displacement.basis(p).domain()->endBdr(boundary::none);

#           ifdef _OPENMP
            const int tid = omp_get_thread_num();
            const int nt  = omp_get_num_threads();
            domIt += tid;
            for ( ; domIt<domItEnd; domIt+=nt )
#           else
            for (; domIt<domItEnd && continueIt; ++domIt )
#           endif
            {
                genSamplingPoints(domIt.lowerCorner(),domIt.upperCorner(),quRule,points);
                mdG.points = points;
                mdU.points = points;
                domain.patch(p).computeMap(mdG);
                displacement.patch(p).computeMap(mdU);
                for (index_t q = 0; q < points.cols() && continueIt; ++q)
                {
                    gsMatrix<T> physDispJac = mdU.jacobian(q)*(mdG.jacobian(q).cramerInverse());
                    gsMatrix<T> F = gsMatrix<T>::Identity(domain.dim(),domain.dim()) + physDispJac;
                    if (F.determinant() <= 0)
                    {
                        gsInfo << "Bad patch: " << p << "\nBad point:\n" << points.col(q) << "\nDet: " << F.determinant() << std::endl;
                        corruptedPatch = p;
                        continueIt = false;
                    }
                }
            }
        }
    }
    return corruptedPatch;
}

template <class T>
T normL2(gsMultiPatch<T> const & domain, gsMultiPatch<T> const & solution)
{
    T norm = 0;
    for (size_t p = 0; p < domain.nPatches(); ++p)
    {
#pragma omp parallel
        {
            gsMapData<T> mdGeo(NEED_MEASURE);
            gsMatrix<T> values, quNodes;
            gsVector<T> quWeights;

            gsVector<index_t> numNodes(domain.dim());
            for (short_t i = 0; i < domain.dim(); ++i)
                numNodes.at(i) = solution.basis(p).degree(i)+1;
            gsQuadRule<T> quRule = gsQuadrature::get<T>(gsQuadrature::GaussLegendre,numNodes);

            // Initialize domain element iterator
            typename gsBasis<T>::domainIter domIt    = solution.basis(p).domain()->beginBdr(boundary::none);
            typename gsBasis<T>::domainIter domItEnd = solution.basis(p).domain()->endBdr(boundary::none);

#           ifdef _OPENMP
            const int tid = omp_get_thread_num();
            const int nt  = omp_get_num_threads();
            domIt += tid;
            for ( ; domIt<domItEnd; domIt+=nt )
#           else
            for (; domIt<domItEnd; ++domIt )
#           endif
            {
                quRule.mapTo( domIt.lowerCorner(), domIt.upperCorner(), quNodes, quWeights );
                mdGeo.points = quNodes;
                domain.patch(p).computeMap(mdGeo);
                solution.patch(p).eval_into(quNodes,values);
                T tempNorm = 0;
                for (index_t q = 0; q < quNodes.cols(); ++q)
                    tempNorm += mdGeo.measure(q)*quWeights.at(q)*(values.col(q).transpose()*values.col(q))(0,0);
#pragma omp critical
                norm += tempNorm;
            }
        }
    }
    return sqrt(norm);
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
        gsQuadRule<T> quRule = gsQuadrature::get<T>(gsQuadrature::GaussLegendre,numNodes);

        // Initialize domain element iterator
        typename gsBasis<T>::domainIter domIt    = domain.basis(p).domain()->beginBdr(boundary::none);
        typename gsBasis<T>::domainIter domItEnd = domain.basis(p).domain()->endBdr(boundary::none);
        for (; domIt<domItEnd; ++domIt )
        {
            genSamplingPoints(domIt.lowerCorner(),domIt.upperCorner(),quRule,points);
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

template <class T>
T displacementJacRatio(const gsMultiPatch<T> & domain,const gsMultiPatch<T> & displacement)
{
    std::vector<T> maxs, mins;
    gsMapData<T> mdG, mdU;
    mdG.flags = NEED_DERIV;
    mdU.flags = NEED_DERIV;
    gsMatrix<T> points;
    size_t dim = domain.parDim();

    for (size_t p = 0; p < displacement.nPatches(); ++p)
    {
        gsVector<index_t> numNodes(domain.dim());
        for (short_t i = 0; i < domain.dim(); ++i)
            numNodes.at(i) = displacement.basis(p).degree(i)+1;
        gsQuadRule<T> quRule = gsQuadrature::get<T>(gsQuadrature::GaussLegendre,numNodes);

        // Initialize domain element iterator
        typename gsBasis<T>::domainIter domIt    = domain.basis(p).domain()->beginBdr(boundary::none);
        typename gsBasis<T>::domainIter domItEnd = domain.basis(p).domain()->endBdr(boundary::none);
        for (; domIt<domItEnd; ++domIt )
        {
            genSamplingPoints(domIt.lowerCorner(),domIt.upperCorner(),quRule,points);
            mdG.points = points;
            mdU.points = points;
            domain.patch(p).computeMap(mdG);
            displacement.patch(p).computeMap(mdU);

            T minJ = (gsMatrix<T>::Identity(dim,dim) + mdU.jacobian(0)*(mdG.jacobian(0).cramerInverse())).determinant();
            T maxJ = minJ;
            for (int q = 1; q < points.cols(); ++q)
            {
                T J = (gsMatrix<T>::Identity(dim,dim) + mdU.jacobian(q)*(mdG.jacobian(q).cramerInverse())).determinant();
                if (J > maxJ)
                    maxJ = J;
                if (J < minJ)
                    minJ = J;
            }

            maxs.push_back(maxJ);
            mins.push_back(minJ);
        }
    }

    return *(std::min_element(mins.begin(),mins.end())) / *(std::max_element(maxs.begin(),maxs.end()));
}

template <class T>
void genSamplingPoints(const gsVector<T> & lower, const gsVector<T> & upper,
                       const gsQuadRule<T> & quRule, gsMatrix<T> & points)
{
    gsMatrix<T> quadPoints;
    gsVector<T> tempVector; // temporary argument for the gsQuadrule::mapTo function
    quRule.mapTo(lower,upper,quadPoints,tempVector);

    gsVector<unsigned> nPoints(quadPoints.rows());
    for (index_t d = 0; d < quadPoints.rows(); ++d)
        nPoints.at(d) = 2;
    gsMatrix<T> corners = gsPointGrid(lower,upper,nPoints);

    points.resize(quadPoints.rows(),quadPoints.cols()+corners.cols());
    points << quadPoints,corners;
}

template <class T>
T patchLength(const gsGeometry<T> & geo, short_t dir)
{
    GISMO_ENSURE(dir >= 0 && dir <= geo.parDim(),"Invalid parametric direction: " + util::to_string(dir)+
                                                 ", geometry dimension: " + util::to_string(geo.parDim()));
    switch (geo.parDim())
    {
        case 1: return curveLength<T>(geo);
        case 2:
        {
            switch (dir)
            {
                case 0:
                {
                    typename gsGeometry<T>::uPtr sideS = geo.boundary(boundary::south);
                    typename gsGeometry<T>::uPtr sideN = geo.boundary(boundary::north);
                    return (curveLength<T>(*sideS) + curveLength<T>(*sideN))/2;
                }
                case 1:
                {
                    typename gsGeometry<T>::uPtr sideW = geo.boundary(boundary::west);
                    typename gsGeometry<T>::uPtr sideE = geo.boundary(boundary::east);
                    return (curveLength<T>(*sideW) + curveLength<T>(*sideE))/2;
                }
                default: return 0.;
            }
        }
        case 3:
        {
            switch (dir)
            {
                case 0:
                {
                    typename gsGeometry<T>::uPtr sideB = geo.boundary(boundary::back);
                    typename gsGeometry<T>::uPtr sideF = geo.boundary(boundary::front);
                    typename gsGeometry<T>::uPtr edgeBS = sideB->boundary(boundary::south);
                    typename gsGeometry<T>::uPtr edgeBN = sideB->boundary(boundary::north);
                    typename gsGeometry<T>::uPtr edgeFS = sideF->boundary(boundary::south);
                    typename gsGeometry<T>::uPtr edgeFN = sideF->boundary(boundary::north);
                    return (curveLength<T>(*edgeBS) + curveLength<T>(*edgeBN) +
                            curveLength<T>(*edgeFS) + curveLength<T>(*edgeFN))/4;
                }
                case 1:
                {
                    typename gsGeometry<T>::uPtr sideB = geo.boundary(boundary::back);
                    typename gsGeometry<T>::uPtr sideF = geo.boundary(boundary::front);
                    typename gsGeometry<T>::uPtr edgeBW = sideB->boundary(boundary::west);
                    typename gsGeometry<T>::uPtr edgeBE = sideB->boundary(boundary::east);
                    typename gsGeometry<T>::uPtr edgeFW = sideF->boundary(boundary::west);
                    typename gsGeometry<T>::uPtr edgeFE = sideF->boundary(boundary::east);
                    return (curveLength<T>(*edgeBW) + curveLength<T>(*edgeBE) +
                            curveLength<T>(*edgeFW) + curveLength<T>(*edgeFE))/4;
                }
                case 2:
                {
                    typename gsGeometry<T>::uPtr sideN = geo.boundary(boundary::north);
                    typename gsGeometry<T>::uPtr sideS = geo.boundary(boundary::south);
                    typename gsGeometry<T>::uPtr edgeNW = sideN->boundary(boundary::west);
                    typename gsGeometry<T>::uPtr edgeNE = sideN->boundary(boundary::east);
                    typename gsGeometry<T>::uPtr edgeSW = sideS->boundary(boundary::west);
                    typename gsGeometry<T>::uPtr edgeSE = sideS->boundary(boundary::east);
                    return (curveLength<T>(*edgeNW) + curveLength<T>(*edgeNE) +
                            curveLength<T>(*edgeSW) + curveLength<T>(*edgeSE))/4;
                }
                default: return 0.;
            }
        }
        default: return 0.;
    }
}

template <class T>
T curveLength(const gsGeometry<T> & geo)
{
    GISMO_ENSURE(geo.parDim() == 1,"This is not a curve! Dim: " + util::to_string(geo.parDim()));
    T length = 0.;

    gsVector<index_t> numNodes(1);
    numNodes << geo.basis().degree(0)+1;
    gsQuadRule<T> quRule = gsQuadrature::get<T>(gsQuadrature::GaussLegendre,numNodes);
    gsMatrix<T> qPoints;
    gsVector<T> qWeights;
    gsMapData<T> md;
    md.flags = NEED_DERIV;

    // Initialize domain element iterator
    typename gsBasis<T>::domainIter domIt    = geo.basis().domain()->beginBdr(boundary::none);
    typename gsBasis<T>::domainIter domItEnd = geo.basis().domain()->endBdr(boundary::none);
    for (; domIt<domItEnd; ++domIt )
    {
        quRule.mapTo(domIt.lowerCorner(),domIt.upperCorner(),qPoints,qWeights);
        md.points = qPoints;
        geo.computeMap(md);
        for (index_t q = 0; q < qWeights.rows(); ++q)
            length += qWeights.at(q)*md.jacobian(q).norm();
    }
    return length;
}

template <class T>
gsVector<unsigned> distributePoints(const gsGeometry<T> & geo, unsigned numPoints)
{
    short_t dim = geo.parDim();
    GISMO_ENSURE(dim >= 1 && dim <= 3,"Invalid patch dimension: " + util::to_string(dim));
    gsVector<unsigned> numPointsPerDir(dim);
    gsVector<T> parLengths(dim);

    T volume = 1.;
    for (short_t d = 0; d < dim; ++d)
    {
        parLengths.at(d) = patchLength<T>(geo,d);
        volume *= parLengths.at(d);
    }

    T unit = pow(1.0*numPoints/volume,1./geo.parDim());
    for (short_t d = 0; d < geo.parDim(); ++d)
        numPointsPerDir.at(d) = math::ceil(unit*parLengths.at(d)) > 1 ? math::ceil(unit*parLengths.at(d)) : 2;
    return numPointsPerDir;
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
    index_t num = basis.size();
    index_t dim = points.rows();

    gsSparseMatrix<T> A(num,num);
    gsMatrix<T> b(num,dim);
    b.setZero();

    gsMatrix<T> basisValues;
    gsMatrix<index_t> activeBasis;

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
                                           util::to_string(A.parDim()) + " and " + util::to_string(B.parDim()) + "\n");
    short_t pDim = A.parDim();
    GISMO_ASSERT(pDim == 1 || pDim ==2, "Can only interpolate between curves or surfaces. Given geometries have parametric dimension " +
                                        util::to_string(pDim) + "\n");
    for (index_t d = 0; d < pDim; ++d)
        GISMO_ENSURE(A.degree(d) == B.degree(d), "Geometries are incompatible: different splines degrees in dimension" +
                                                 util::to_string(d) + ": " + util::to_string(A.degree(d)) +
                                                 " and " + util::to_string(B.degree(d)) + "\n");

    GISMO_ENSURE(A.targetDim() == B.targetDim(), "Geometries are incompatible: different physical dimensions: " +
                                                 util::to_string(A.targetDim()) + " and " + util::to_string(B.targetDim()) + "\n");
    short_t tDim = A.targetDim();
    GISMO_ASSERT(A.coefsSize() == B.coefsSize(), "Geometries are incompatible: different number of control points: " +
                                                 util::to_string(A.coefsSize()) + " and " + util::to_string(B.coefsSize()) + "\n");
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

    gsKnotVector<T> knotsCoarse(0.0,1.0, 0, 2);
    index_t dim = A.cols();
    gsMatrix<T> coefs(2,dim);
    coefs.row(0) = A.row(iA);
    coefs.row(1) = B.row(iB);
    // gen simplest possible Bspline line
    typename gsBSpline<T>::uPtr line1 = typename gsBSpline<T>::uPtr(new gsBSpline<T>(knotsCoarse,coefs));
    // refine it
    for (index_t i = 0; i < deg - 1; ++i)
        line1->degreeElevate();
    gsKnotVector<T> knotsFine(0.0,1.0, num - deg - 1, deg + 1);
    for (index_t i = 0; i < num - deg - 1; ++i)
        line1->insertKnot(knotsFine.at(i+deg+1));

    return line1;
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
    GISMO_UNUSED(eta0);
    GISMO_UNUSED(eta1);
    gsBSplineBasis<T> xiBasis(xiKnots);
    typename gsGeometry<T>::uPtr xiCircle = genCircle(xiBasis,(T)1.,(T)0.,(T)0.,xi0,xi1-xi0);
    gsBSplineBasis<T> etaBasis(etaKnots);
    typename gsGeometry<T>::uPtr etaCircle = genCircle(xiBasis,(T)1.,(T)0.,(T)0.,xi0,xi1-xi0);

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

template<class T>
typename gsGeometry<T>::uPtr genSpring(gsGeometry<T> const & crossSection,
                                       T R, T springPitch,
                                       index_t N, bool nurbs)
{
    GISMO_ENSURE(crossSection.parDim() == 2, "Wrong dimension of the cross-section: " + util::to_string(crossSection.parDim()));

    // cross-section knot vectors
    gsKnotVector<T> xKnots = nurbs ? static_cast<const gsTensorNurbsBasis<2,T> &>(crossSection.basis()).knots(0)
                                   : static_cast<const gsTensorBSplineBasis<2,T> &>(crossSection.basis()).knots(0);
    gsKnotVector<T> yKnots = nurbs ? static_cast<const gsTensorNurbsBasis<2,T> &>(crossSection.basis()).knots(1)
                                   : static_cast<const gsTensorBSplineBasis<2,T> &>(crossSection.basis()).knots(1);
    // knot vector for the helical guiding line
    std::vector<T> knots;
    for (index_t i = 0; i < 3; ++i)
        knots.push_back(0.);
    for (index_t i = 0; i < N-1; ++i)
        for (index_t j = 0; j < 2; ++j)
            knots.push_back(1./N*(i+1));
    for (index_t i = 0; i < 3; ++i)
        knots.push_back(1.);
    gsKnotVector<T> zKnots(knots,2);

    // NURBS weights
    index_t numP = crossSection.coefs().rows();
    gsMatrix<T> crossWeights = nurbs ? static_cast<const gsTensorNurbsBasis<2,T> &>(crossSection.basis()).weights()
                                     : gsMatrix<T>::Ones(numP,1);
    gsMatrix<T> weights(numP*(2*N+1),1);
    for (index_t i = 0; i < 2*N+1; ++i)
        weights.middleRows(i*numP,numP) = crossWeights*(i%2==0 ? 1 : 1./sqrt(2));

    // control points of the cross-section
    gsMatrix<T> crossCoefs(numP,3);
    crossCoefs.setZero();
    crossCoefs.col(0) = crossSection.coefs().col(0);
    crossCoefs.col(2) = crossSection.coefs().col(1);
    for (index_t i = 0; i < numP; ++i)
        crossCoefs(i,0) += R;

    // control points of the spring
    gsMatrix<T> coefs(numP*(2*N+1),3);
    for (index_t i = 0; i < 2*N+1; ++i)
    {
        gsMatrix<T> rotation = gsEigen::AngleAxis<T>(-1.0*EIGEN_PI/4*i,gsEigen::Vector3f(0.,0.,1.)).toRotationMatrix().transpose();
        gsMatrix<T> scaling = gsEigen::DiagonalMatrix<T,3>((i%2 == 0 ? 1. : sqrt(2)),1.,1.);
        coefs.middleRows(i*numP,numP) = crossCoefs * scaling  * rotation;
        for (index_t j = 0; j < numP; ++j)
            coefs(i*numP + j,2) += i*springPitch/8;
    }

    return typename gsGeometry<T>::uPtr(new gsTensorNurbs<3,T>(xKnots,yKnots,zKnots,coefs,weights));
}

template<class T>
void genMuscleMP(gsGeometry<T> const & muscleSurface, gsMultiPatch<T> & result)
{
    index_t innerDim = 2;

    gsMultiPatch<T> temp;
    temp.addPatch(muscleSurface.clone());
    static_cast<gsTensorBSpline<2,T> &>(temp.patch(0)).insertKnot(
                static_cast<const gsTensorBSplineBasis<2,T> &>(muscleSurface.basis()).knots(1).at(4),1);
    static_cast<gsTensorBSpline<2,T> &>(temp.patch(0)).insertKnot(
                static_cast<const gsTensorBSplineBasis<2,T> &>(muscleSurface.basis()).knots(1).at(4),1);
    static_cast<gsTensorBSpline<2,T> &>(temp.patch(0)).insertKnot(0.5,1);
    static_cast<gsTensorBSpline<2,T> &>(temp.patch(0)).insertKnot(0.5,1);
    static_cast<gsTensorBSpline<2,T> &>(temp.patch(0)).insertKnot(
                static_cast<const gsTensorBSplineBasis<2,T> &>(muscleSurface.basis()).knots(1).at(7),1);
    static_cast<gsTensorBSpline<2,T> &>(temp.patch(0)).insertKnot(
                static_cast<const gsTensorBSplineBasis<2,T> &>(muscleSurface.basis()).knots(1).at(7),1);

    std::vector<T> knotsEtaShort;
    for (index_t i = 0; i < 4; ++i)
        knotsEtaShort.push_back(0.);
    for (index_t i = 0; i < 4; ++i)
        knotsEtaShort.push_back(1.);
    gsKnotVector<T> knotVectorEtaShort(knotsEtaShort,3);

    std::vector<T> knotsEtaLong;
    for (index_t i = 0; i < 4; ++i)
        knotsEtaLong.push_back(0.);
    knotsEtaLong.push_back(0.5);
    for (index_t i = 0; i < 4; ++i)
        knotsEtaLong.push_back(1.);
    gsKnotVector<T> knotVectorEtaLong(knotsEtaLong,3);

    std::vector<T> knotsInner;
    for (index_t i = 0; i < innerDim+1; ++i)
        knotsInner.push_back(0.);
    for (index_t i = 0; i < innerDim+1; ++i)
        knotsInner.push_back(1.);
    gsKnotVector<T> knotVectorInner(knotsInner,innerDim);

    gsTensorBSplineBasis<3,T> basis1(static_cast<const gsTensorBSplineBasis<2,T> &>(temp.basis(0)).knots(0),
                                     knotVectorEtaShort,knotVectorInner);
    gsTensorBSplineBasis<3,T> basis2(static_cast<const gsTensorBSplineBasis<2,T> &>(temp.basis(0)).knots(0),
                                     knotVectorEtaLong,knotVectorInner);
    gsTensorBSplineBasis<3,T> basis3(static_cast<const gsTensorBSplineBasis<2,T> &>(temp.basis(0)).knots(0),
                                     knotVectorEtaShort,knotVectorInner);
    gsTensorBSplineBasis<3,T> basis4(static_cast<const gsTensorBSplineBasis<2,T> &>(temp.basis(0)).knots(0),
                                     knotVectorEtaLong,knotVectorInner);
    gsTensorBSplineBasis<3,T> basis5(static_cast<const gsTensorBSplineBasis<2,T> &>(temp.basis(0)).knots(0),
                                     knotVectorEtaShort,knotVectorEtaLong);

    gsMatrix<T> coefs1(7*4*(innerDim+1),3);
    gsMatrix<T> coefs2(7*5*(innerDim+1),3);
    gsMatrix<T> coefs3(7*4*(innerDim+1),3);
    gsMatrix<T> coefs4(7*5*(innerDim+1),3);
    gsMatrix<T> coefs5(7*4*5,3);
    for (index_t i = 0; i < 7; ++i)
    {   // for each cross-section
        // compute four corners of the interior hub-patch
        gsMatrix<T> cornerA = combine(temp.patch(0).coefs(),temp.patch(0).coefs(),(T)1./3.,0*7+i,7*7+i); // right 1/3 left 1/3
        gsMatrix<T> cornerC = combine(temp.patch(0).coefs(),temp.patch(0).coefs(),(T)2./3.,0*7+i,7*7+i); // right 3/5 left 2/3
        gsMatrix<T> cornerB = combine(temp.patch(0).coefs(),temp.patch(0).coefs(),(T)1./3.,3*7+i,10*7+i); // right 2/5 left 1/3
        gsMatrix<T> cornerD = combine(temp.patch(0).coefs(),temp.patch(0).coefs(),(T)2./3.,3*7+i,10*7+i); // right 2/3 left 2/3

        if (i == 0) // only for the left
        {
            cornerB = combine(cornerB,cornerC,(T)0.5);
            cornerD = combine(cornerD,cornerC,(T)0.5);
        }

        // gen patch1
        gsMultiPatch<T> patch1Curves;
        patch1Curves.addPatch(genLine(innerDim,innerDim+1,temp.patch(0).coefs(),cornerA,0*7+i));
        patch1Curves.addPatch(genLine(innerDim,innerDim+1,temp.patch(0).coefs(),cornerB,3*7+i));
        gsMatrix<T> bdryABcoefs(4,3);
        for (index_t j = 0; j < 4; ++j)
            bdryABcoefs.row(j) = temp.patch(0).coefs().row(j*7+i);
        patch1Curves.addPatch(typename gsGeometry<T>::uPtr(new gsBSpline<T>(knotVectorEtaShort,bdryABcoefs)));
        patch1Curves.addPatch(genLine(3,4,cornerA,cornerB));
        gsCoonsPatch<T> coonsPatch1(patch1Curves);
        coonsPatch1.compute();

        for (index_t j = 0; j < 4; ++j)
            for (index_t k = 0; k < innerDim+1; ++k)
                coefs1.row(k*4*7+j*7+(6-i)) = coonsPatch1.result().coefs().row(k*4+j);

        // gen patch2
        gsMultiPatch<T> patch2Curves;
        patch2Curves.addPatch(genLine(innerDim,innerDim+1,temp.patch(0).coefs(),cornerB,3*7+i));
        patch2Curves.addPatch(genLine(innerDim,innerDim+1,temp.patch(0).coefs(),cornerC,7*7+i));
        gsMatrix<T> bdryBCcoefs(5,3);
        for (index_t j = 0; j < 5; ++j)
            bdryBCcoefs.row(j) = temp.patch(0).coefs().row((j+3)*7+i);
        patch2Curves.addPatch(typename gsGeometry<T>::uPtr(new gsBSpline<T>(knotVectorEtaLong,bdryBCcoefs)));
        patch2Curves.addPatch(genLine(3,5,cornerB,cornerC));
        gsCoonsPatch<T> coonsPatch2(patch2Curves);
        coonsPatch2.compute();

        for (index_t j = 0; j < 5; ++j)
            for (index_t k = 0; k < innerDim+1; ++k)
                coefs2.row(k*5*7+j*7+(6-i)) = coonsPatch2.result().coefs().row(k*5+j);

        // gen patch3
        gsMultiPatch<T> patch3Curves;
        patch3Curves.addPatch(genLine(innerDim,innerDim+1,temp.patch(0).coefs(),cornerC,7*7+i));
        patch3Curves.addPatch(genLine(innerDim,innerDim+1,temp.patch(0).coefs(),cornerD,10*7+i));
        gsMatrix<T> bdryCDcoefs(4,3);
        for (index_t j = 0; j < 4; ++j)
            bdryCDcoefs.row(j) = temp.patch(0).coefs().row((j+7)*7+i);
        patch3Curves.addPatch(typename gsGeometry<T>::uPtr(new gsBSpline<T>(knotVectorEtaShort,bdryCDcoefs)));
        patch3Curves.addPatch(genLine(3,4,cornerC,cornerD));
        gsCoonsPatch<T> coonsPatch3(patch3Curves);
        coonsPatch3.compute();

        for (index_t j = 0; j < 4; ++j)
            for (index_t k = 0; k < innerDim+1; ++k)
                coefs3.row(k*4*7+j*7+(6-i)) = coonsPatch3.result().coefs().row(k*4+j);

        // gen patch4
        gsMultiPatch<T> patch4Curves;
        patch4Curves.addPatch(genLine(innerDim,innerDim+1,temp.patch(0).coefs(),cornerD,10*7+i));
        patch4Curves.addPatch(genLine(innerDim,innerDim+1,temp.patch(0).coefs(),cornerA,14*7+i));
        gsMatrix<T> bdryDAcoefs(5,3);
        for (index_t j = 0; j < 5; ++j)
            bdryDAcoefs.row(j) = temp.patch(0).coefs().row((j+10)*7+i);
        patch4Curves.addPatch(typename gsGeometry<T>::uPtr(new gsBSpline<T>(knotVectorEtaLong,bdryDAcoefs)));
        patch4Curves.addPatch(genLine(3,5,cornerD,cornerA));
        gsCoonsPatch<T> coonsPatch4(patch4Curves);
        coonsPatch4.compute();

        for (index_t j = 0; j < 5; ++j)
            for (index_t k = 0; k < innerDim+1; ++k)
                coefs4.row(k*5*7+j*7+(6-i)) = coonsPatch4.result().coefs().row(k*5+j);

        // gen patch5
        gsMultiPatch<T> patch5Curves;
        patch5Curves.addPatch(genLine(3,5,cornerA,cornerD));
        patch5Curves.addPatch(genLine(3,5,cornerB,cornerC));
        patch5Curves.addPatch(genLine(3,4,cornerA,cornerB));
        patch5Curves.addPatch(genLine(3,4,cornerD,cornerC));

        gsCoonsPatch<T> coonsPatch5(patch5Curves);
        coonsPatch5.compute();

        for (index_t j = 0; j < 4; ++j)
            for (index_t k = 0; k < 5; ++k)
                coefs5.row(k*4*7+j*7+(6-i)) = coonsPatch5.result().coefs().row(k*4+j);

    }

    result.addPatch(basis1.makeGeometry(coefs1));
    result.addPatch(basis2.makeGeometry(coefs2));
    result.addPatch(basis3.makeGeometry(coefs3));
    result.addPatch(basis4.makeGeometry(coefs4));
    result.addPatch(basis5.makeGeometry(coefs5));

    result.computeTopology();
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
