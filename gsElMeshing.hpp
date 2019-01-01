/** @file gsElMeshing.hpp

    @brief Provides isogeometric meshing and modelling routines.

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s):
        A.Shamanskiy (2016 - ...., TU Kaiserslautern)
*/

#pragma once

#include <gsElasticity/gsElMeshing.h>

#include <gsNurbs/gsTensorBSplineBasis.h>
#include <gsModeling/gsFitting.h>

#include <gsElasticity/gsElasticityAssembler.h>
#include <gsElasticity/gsElasticityNewton.h>
#include <gsIO/gsWriteParaview.h>
#include <gsCore/gsFunctionExpr.h>

#include <gsElasticity/gsElasticityFunctions.h>
#include <gsCore/gsField.h>
#include <gsNurbs/gsBSpline.h>

#include <gsCore/gsLinearAlgebra.h>
#include <iostream>
#include <fstream>


namespace gismo
{

//-----------------------------------//
//--------- Mesh deformation --------//
//-----------------------------------//



template <class T>
void computeDeformationInc(std::vector<std::vector<gsMatrix<T> > > & result,
            gsMultiPatch<T> const & domain, gsBoundaryConditions<T> const & bdry,
            int numSteps, T poissonRatio)
{
    int numP = domain.nPatches();
    gsMultiPatch<T> geo;
    for (int p = 0; p < numP; ++p)
        geo.addPatch(domain.patch(p));

    result.resize(numSteps);
    for (int s = 0; s < numSteps; ++s)
        result[s].resize(numP);

    int pDim = domain.parDim();

    gsBoundaryConditions<T> bcInfo;
    for (auto it = domain.bBegin(); it != domain.bEnd(); ++it)
        for (int d = 0; d < pDim; ++d)
            bcInfo.addCondition(it->patch,it->side(),condition_type::dirichlet,0,d);

    std::map<std::pair<int,int>,gsMatrix<T> > deformCoefs;
    for (auto it = bdry.dirichletBegin(); it != bdry.dirichletEnd(); ++it)
        deformCoefs[std::pair<int,int>(it->patch(),it->side())] = (static_cast<gsGeometry<T> &>(*(it->function())).coefs() -

                                                                   domain.patch(it->patch()).boundary(it->side())->coefs()) / numSteps;
    std::vector<std::string> zeros;
    for (int d = 0; d < pDim; ++d)
        zeros.push_back("0.");
    gsFunctionExpr<T> g(zeros,pDim);

    for (int s = 0; s < numSteps; ++s)
    {
        gsInfo << "Incremental step " << s+1 << "/" << numSteps << "..." << std::endl;

        gsMultiBasis<T> basis(geo);
        gsElasticityAssembler<T> assembler(geo,basis,bcInfo,g);
        assembler.options().setReal("PoissonsRatio",poissonRatio);
        for (auto it = bdry.dirichletBegin(); it != bdry.dirichletEnd(); ++it)
            assembler.setDirichletDofs(it->patch(),it->side(),deformCoefs.at(std::pair<int,int>(it->patch(),it->side())));

        assembler.assemble();
        gsSparseSolver<>::LU solver(assembler.matrix());
        gsVector<> solVector = solver.solve(assembler.rhs());

        assembler.deformGeometry(solVector,geo);
        geo.computeTopology();

        gsMultiPatch<T> displacement;
        assembler.constructSolution(solVector,displacement);

        for (int p = 0; p < numP; ++p)
            result[s][p] = displacement.patch(p).coefs();
    }
}

template <class T>
void deformNonLinearly(gsMultiPatch<T> & result, gsMultiPatch<T> const & domain,
                       std::vector<std::vector<gsMatrix<T> > > const & deformation,
                       gsMatrix<T> const & sumSolutionVector,
                       gsBoundaryConditions<T> const & bdry, T poissonRatio,
                       T tolerance, int maxNumIterations)
{/*
    result.clear();
    for (int p = 0; p < domain.nPatches(); ++p)
    {
        result.addPatch(domain.patch(p).clone());
        for (unsigned s = 0; s < deformation.size(); ++s)
            result.patch(p).coefs() += deformation[s][p];
    }

    gsMultiBasis<T> basis(domain);

    gsBoundaryConditions<T> bcInfo;
    int pDim = domain.parDim();
    for (auto it = domain.bBegin(); it != domain.bEnd(); ++it)
        for (int d = 0; d < pDim; ++d)
            bcInfo.addCondition(it->patch,it->side(),condition_type::dirichlet,0,d);

    std::map<std::pair<int,int>,gsMatrix<T> > deformCoefs;
    for (auto it = bdry.dirichletBegin(); it != bdry.dirichletEnd(); ++it)
        deformCoefs[std::pair<int,int>(it->patch(),it->side())] = (static_cast<gsGeometry<T> &>(*(it->function())).coefs() -
                                                                   domain.patch(it->patch()).boundary(it->side())->coefs());


    gsConstantFunction<T> g(0.,0.,2);
    gsElasticityAssembler<T> assembler(domain,basis,200.,poissonRatio,1.,bcInfo,g);
    for (auto it = bdry.dirichletBegin(); it != bdry.dirichletEnd(); ++it)
        assembler.setDirichletDoFs(deformCoefs.at(std::pair<int,int>(it->patch(),it->side())),it->patch(),it->side());

    assembler.set_MaterialLaw(0);

    gsElasticityNewton<T> newtonSolver(assembler,result);
    newtonSolver.setMaxIterations(maxNumIterations);
    newtonSolver.setTolerance(tolerance);
    newtonSolver.setSolution(sumSolutionVector);

    newtonSolver.solveGiven();
    result.clear();
    for (int p = 0; p < domain.nPatches(); ++p)
        result.addPatch(newtonSolver.solution().patch(p).clone());*/

}

template <class T>
void plotDeformation(std::vector<std::vector<gsMatrix<T> > > const & deformation,
                     gsMultiPatch<T> const & domain, std::string fileName,
                     bool plotJac, int numSamples)
{
    std::string fileNameOnly = fileName.substr(fileName.find_last_of("/\\")+1);
    gsParaviewCollection collectionMesh(fileName + "_mesh");
    gsParaviewCollection collectionJac(fileName + "_jac");

    gsMultiPatch<T> configuration;
    for (int p = 0; p < domain.nPatches(); ++p)
        configuration.addPatch(domain.patch(p));

    gsPiecewiseFunction<T> dets;
    for (int p = 0; p < configuration.nPatches(); ++p)
        dets.addPiecePointer(new gsDetFunction<T>(configuration,p));

    int samples = numSamples;
    if (!plotJac)
        samples = 1;

    gsField<> detField(configuration,dets,true);
    std::map<std::string,const gsField<> *> fields;
    fields["Jacobian"] = &detField;
    gsWriteParaviewMultiPhysics(fields,fileName+std::to_string(0),samples,true);

    gsInfo << "Plotting initial configuration...\n";
    for (int p = 0; p < domain.nPatches(); ++p)
    {
        collectionMesh.addTimestep(fileNameOnly + std::to_string(0),p,0,"_mesh.vtp");
        int res = system(("rm " + fileName + std::to_string(0) + ".pvd").c_str());
        if (plotJac)
            collectionJac.addTimestep(fileNameOnly + std::to_string(0),p,0,".vts");
        else
            res = system(("rm " + fileName + std::to_string(0) + std::to_string(p) + ".vts").c_str());
        GISMO_ASSERT(res == 0, "Problems with deleting files\n");
    }

    for (unsigned s = 0; s < deformation.size(); ++s)
    {
        gsInfo << "Plotting configuration " << s+1 << "/" << deformation.size() << "..." << std::endl;

        for (int p = 0; p < domain.nPatches(); ++p)
            configuration.patch(p).coefs() += deformation[s][p];

        gsWriteParaviewMultiPhysics(fields,fileName+std::to_string(s+1),samples,true);
        for (int p = 0; p < domain.nPatches(); ++p)
        {
            collectionMesh.addTimestep(fileNameOnly + std::to_string(s+1),p,s+1,"_mesh.vtp");
            int res = system(("rm " + fileName + std::to_string(s+1) + ".pvd").c_str());
            if (plotJac)
                collectionJac.addTimestep(fileNameOnly + std::to_string(s+1),p,s+1,".vts");
            else
                res = system(("rm " + fileName + std::to_string(s+1) + std::to_string(p) + ".vts").c_str());

            GISMO_ASSERT(res == 0, "Problems with deleting files\n");
        }
    }

    collectionMesh.save();
    if (plotJac)
        collectionJac.save();
}

template <class T>
void applyDeformation(std::vector<std::vector<gsMatrix<T> > > const & deformation,
                      gsMultiPatch<T> const & initDomain, gsMultiPatch<T> & domain)
{
    domain.clear();
    for (int p = 0; p < initDomain.nPatches(); ++p)
    {
        domain.addPatch(initDomain.patch(p).clone());
        for (unsigned s = 0; s < deformation.size(); ++s)
            domain.patch(p).coefs() += deformation[s][p];
    }
}

template <class T>
T measureMinMaxJ(gsMultiPatch<T> const & domain, int numSamples)
{
    /*std::vector<T> maxs, mins;

    for (int p = 0; p < domain.nPatches(); ++p)
    {
        gsMatrix<> ab = domain.patch(p).support();
        gsVector<> a = ab.col(0);
        gsVector<> b = ab.col(1);
        gsVector<unsigned> np = uniformSampleCount(a,b,numSamples);
        gsMatrix<> pts = gsPointGrid(a,b,np);

        typename gsGeometryEvaluator<T>::uPtr geoEval(getEvaluator(NEED_MEASURE,domain.patch(p)));
        geoEval->evaluateAt(pts);
        T min = abs(geoEval->jacDet(0));
        T max = geoEval->jacDet(0);
        for (int i = 1; i < pts.cols(); ++i)
        {
            if (geoEval->jacDet(i) > max)
                max = geoEval->jacDet(i);
            if (abs(geoEval->jacDet(i)) < min)
                min = abs(geoEval->jacDet(i));
        }

        maxs.push_back(max);
        mins.push_back(min);
    }

    return *(std::min_element(mins.begin(),mins.end())) / *(std::max_element(maxs.begin(),maxs.end()));*/
    return 0.;
}

template <class T>
void analyzeDeformation(std::vector<std::vector<gsMatrix<T> > > const & deformation,
                        gsMultiPatch<T> const & domain, int measPerStep,
                        std::string fileName, int numSamples)
{/*
    int steps = deformation.size();
    std::ofstream res;
    res.open(fileName + "_" + std::to_string(steps) + ".res");
    std::ofstream resFine;
    resFine.open(fileName + "_" + std::to_string(steps) + ".frs");

    T initialMeasure = measure(domain);
    res << 1. << std::endl;
    resFine << 1. << std::endl;

    gsMultiPatch<T> configuration;
    for (int p = 0; p < domain.nPatches(); ++p)
        configuration.addPatch(domain.patch(p));

    T newMeasure = 1.;
    for (int s = 0; s < steps; ++s)
    {
        gsInfo << "Analyzing configuration " << s+1 << "/" << deformation.size() << "..." << std::endl;
        for (int i = 0; i < measPerStep; ++i)
        {
            for (int p = 0; p < domain.nPatches(); ++p)
                configuration.patch(p).coefs() += deformation[s][p];

            newMeasure = measure(configuration);
            resFine << newMeasure/initialMeasure << std::endl;
        }
        res << newMeasure/initialMeasure << std::endl;
    }

    res.close();
    resFine.close();*/
}

//--------------------------------------------------------------------//
//------------------------ Modelling ---------------------------------//
//--------------------------------------------------------------------//

template<class T>
typename gsGeometry<T>::uPtr simplifyCurve(gsGeometry<T> const & curve,
                                          int additionalPoints, int numSamples)
{
    GISMO_ASSERT(curve.domainDim() == 1 ,"That's not a curve.\n");
    int deg = curve.degree(0);
    int num = deg + 1 + additionalPoints;
    GISMO_ASSERT(num >= deg+1 ,"Coarse basis size is too small.\n");
    GISMO_ASSERT(unsigned(num) <= curve.coefsSize() ,"Coarse basis if finer than the original.\n");

    gsKnotVector<T> knots(0.0,1.0, num - deg - 1, deg + 1);
    gsKnotVector<T> knotVector(knots);
    gsBSplineBasis<T> basis(knotVector);

    gsMatrix<T> params(1,numSamples);
    for (int p = 0; p < numSamples; ++p)
        params.at(p) = 1.*p/(numSamples-1);

    gsMatrix<T> curveValues;
    curve.eval_into(params,curveValues);

    gsMatrix<T> lens(numSamples,1);
    lens.at(0) = 0.;
    for (int p = 1; p < numSamples; ++p)
        lens.at(p) = lens.at(p-1) + distance(curveValues,p,curveValues,p-1,true);

    gsMatrix<T> lenParams(1,numSamples);
    for (int p = 0; p < numSamples; ++p)
        lenParams.at(p) = lens.at(p)/lens.at(numSamples-1);

    typename gsGeometry<T>::uPtr simpleCurve = fittingDirichlet(lenParams,curveValues,basis);

    simpleCurve->eval_into(params,curveValues);
    for (int p = 1; p < numSamples; ++p)
        lens.at(p) = lens.at(p-1) + distance(curveValues,p,curveValues,p-1,true);

    for (int p = 0; p < numSamples; ++p)
        lenParams.at(p) = lens.at(p)/lens.at(numSamples-1);

    return fittingDirichlet(lenParams,curveValues,curve.basis());
}

template <class T>
typename gsGeometry<T>::uPtr fittingDirichlet(gsMatrix<T> const & params,
                                              gsMatrix<T> const & points,
                                              gsBasis<T> const & basis)
{
    int numSamples = params.cols();
    unsigned num = basis.size();
    int dim = points.rows();

    gsSparseMatrix<T> A(num,num);
    gsMatrix<T> b(num,dim);
    b.setZero();

    gsMatrix<T> basisValues;
    gsMatrix<unsigned> activeBasis;

    basis.eval_into(params,basisValues);
    basis.active_into(params,activeBasis);

    for (int p = 1; p < numSamples; ++p)
    {
        int numActive = activeBasis.rows();
        for (int i = 0; i < numActive; ++i)
        {
            if (activeBasis(i,p) == 0)
            {
                for (int j = 0; j < numActive; ++j)
                {
                    if (activeBasis(j,p) != 0 && activeBasis(j,p) != num-1 )
                        b.row(activeBasis(j,p)-1) -= basisValues(i,p) * basisValues(j,p) *
                                                     points.col(0).transpose() * (params.at(p)-params.at(p-1));
                }
            }
            else if (activeBasis(i,p) == num-1)
            {
                for (int j = 0; j < numActive; ++j)
                {
                    if (activeBasis(j,p) != 0 && activeBasis(j,p) != num-1 )
                        b.row(activeBasis(j,p)-1) -= basisValues(i,p) * basisValues(j,p) *
                                                     points.col(numSamples-1).transpose() * (params.at(p)-params.at(p-1));
                }
            }
            else
            {
                b.row(activeBasis(i,p)-1) += basisValues(i,p) * points.col(p).transpose() * (params.at(p)-params.at(p-1));
                for (int j = 0; j < numActive; ++j)
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
                                                   int deg, int num, bool xiDir)
{
    GISMO_ASSERT(A.parDim() == B.parDim(), "Geometries are incompatible: different parametric dimensions: " +
                                           std::to_string(A.parDim()) + " and " + std::to_string(B.parDim()) + "\n");
    int pDim = A.parDim();
    GISMO_ASSERT(pDim == 1 || pDim ==2, "Can only interpolate between curves or surfaces. Given geometries have parametric dimension " +
                                        std::to_string(pDim) + "\n");
    for (int d = 0; d < pDim; ++d)
        GISMO_ASSERT(A.degree(d) == B.degree(d), "Geometries are incompatible: different splines degrees in dimension" +
                                                 std::to_string(d) + ": " + std::to_string(A.degree(d)) +
                                                 " and " + std::to_string(B.degree(d)) + "\n");

    GISMO_ASSERT(A.targetDim() == B.targetDim(), "Geometries are incompatible: different physical dimensions: " +
                                                 std::to_string(A.targetDim()) + " and " + std::to_string(B.targetDim()) + "\n");
    int tDim = A.targetDim();
    GISMO_ASSERT(A.coefsSize() == B.coefsSize(), "Geometries are incompatible: different number of control points: " +
                                                 std::to_string(A.coefsSize()) + " and " + std::to_string(B.coefsSize()) + "\n");
    int baseNum = A.coefsSize();

    gsKnotVector<T> newKnots(0.0,1.0, num - deg - 1, deg + 1);

    gsMultiPatch<> temp;
    temp.addPatch(A.clone());

    if (pDim == 1)
    {

        gsMatrix<T> coefs(baseNum*num,tDim);
        T part = 1./(num-deg);
        T pos = 0.;
        for (int i = 0; i < num; ++i)
        {
            if (i < deg)
                pos += part/deg*i;
            else if (deg <= i && i < num-deg )
                pos += part;
            else
                pos += part/deg*(num-i);

            for (int j = 0; j < baseNum; ++j)
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
        for (int i = 0; i < num; ++i)
        {
            if (i < deg)
                pos += part/deg*i;
            else if (deg <= i && i < num-deg )
                pos += part;
            else
                pos += part/deg*(num-i);

            for (int j = 0; j < baseNum; ++j)
                coefs.row(i*baseNum+j) = combine(A.coefs(),B.coefs(),pos,j,j).row(0);
        }

        return basis.makeGeometry(give(coefs));
    }
}

template<class T>
typename gsGeometry<T>::uPtr genPatchScaling(gsGeometry<T> const & boundary,
                                             int deg, int num,
                                             T scaling, gsVector<T> const & center)
{
    typename gsGeometry<T>::uPtr scaledBoundary = boundary.clone();
    scaledBoundary->translate(-1*center);
    scaledBoundary->scale(scaling);
    scaledBoundary->translate(center);
    return genPatchInterpolation(boundary,*scaledBoundary,deg,num);
}

template<class T>
typename gsGeometry<T>::uPtr genLine(int deg, int num,
                                     gsMatrix<T> const & A, gsMatrix<T> const & B,
                                     int iA, int iB)
{
    GISMO_ASSERT(num - deg - 1 >= 0,"Too few DoFs\n");
    GISMO_ASSERT(A.cols() == B.cols(),"Points have different dimensions\n");
    gsKnotVector<T> knots(0.0,1.0, num - deg - 1, deg + 1);
    gsBSplineBasis<T> basis(knots);

    int dim = A.cols();
    gsMatrix<T> coefs(num,dim);

    T part = 1./(num-deg);
    T pos = 0.;
    for (int i = 0; i < num; ++i)
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
typename gsGeometry<T>::uPtr genCircle(int deg, int num,
                                       T radius, T x0, T y0,
                                       T angle0, T arcAngle)
{
    GISMO_ASSERT(num - deg - 1 >= 0,"Too few DoFs\n");
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
    int numPoints = 1000;
    gsMatrix<> params(1,numPoints);
    gsMatrix<> points(2,numPoints);
    for (int i = 0; i < numPoints; ++i)
    {
        params(0,i) = 1.*i/(numPoints-1);
        points(0,i) = x0 + radius*cos(angle0 + i*arcAngle/(numPoints-1));
        points(1,i) = y0 + radius*sin(angle0 + i*arcAngle/(numPoints-1));
    }

    return fittingDirichlet(params,points,basis);
}

template<class T>
typename gsGeometry<T>::uPtr genQuad(int xiDeg, int xiNum, int etaDeg, int etaNum,
                                     gsMatrix<T> const & A, gsMatrix<T> const & B,
                                     gsMatrix<T> const & C, gsMatrix<T> const & D,
                                     int iA, int iB, int iC, int iD)
{
    typename gsGeometry<T>::uPtr sideAB = genLine(xiDeg,xiNum,A,B,iA,iB);
    typename gsGeometry<T>::uPtr sideCD = genLine(xiDeg,xiNum,C,D,iC,iD);

    return genPatchInterpolation(*sideAB,*sideCD,etaDeg,etaNum);
}

template<class T>
typename gsGeometry<T>::uPtr genSphere(int xiDeg, int xiNum, int etaDeg, int etaNum,
                                       T xi0, T xi1, T eta0, T eta1)
{
    GISMO_ASSERT(xiNum - xiDeg - 1 >= 0,"Too few DoFs\n");
    GISMO_ASSERT(etaNum - etaDeg - 1 >= 0,"Too few DoFs\n");

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

    int xiNum = xiCircle->coefsSize();
    int etaNum = etaCircle->coefsSize();
    gsMatrix<T> coefs(xiNum*etaNum,3);
    for (int i = 0; i < xiNum; ++i)
    {
        for (int j = 0; j < etaNum; ++j)
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
                                         int deg, int num, T height)
{
    GISMO_ASSERT(num - deg - 1 >= 0,"Too few DoFs\n");

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
                                      int deg, int num,
                                      T height, T pitch, T x0, T y0)
{
    GISMO_ASSERT(num - deg - 1 >= 0,"Too few DoFs\n");

    int pDim = base.parDim();
    GISMO_ASSERT(pDim == 1 || pDim ==2,"Wrong geometry type\n");

    gsKnotVector<> zKnots(0.0,1.0, num - deg - 1, deg + 1);
    gsBSplineBasis<> zBasis(zKnots);

    int  numBase = base.coefsSize();

    int numPoints = 10000;
    gsMatrix<> params(1,numPoints);
    gsMatrix<> points(3,numPoints);
    for (int i = 0; i < numPoints; ++i)
    {
        params(0,i) = 1.*i/(numPoints-1);
        points(0,i) = cos(params(0,i)*2*M_PI*pitch/360);
        points(1,i) = sin(params(0,i)*2*M_PI*pitch/360);
        points(2,i) = height*params(0,i);
    }

    typename gsGeometry<T>::uPtr helix = fittingDirichlet(params,points,zBasis);

    gsMatrix<T> coefs(numBase*num,3);
    gsMultiPatch<T> temp;
    temp.addPatch(base.clone());

    T oldAngle = 0.;
    T oldRadius = 1.;
    for (int i = 0; i < num; ++i)
    {
        T x = helix->coef(i,0);
        T y = helix->coef(i,1);
        T radius = sqrt(x*x+y*y);
        T angle = atan2(y,x);
        temp.patch(0).translate(gsVector<T,2>::vec(-1*x0,-1*y0));
        temp.patch(0).rotate(angle-oldAngle);
        temp.patch(0).scale(radius/oldRadius);
        temp.patch(0).translate(gsVector<T,2>::vec(x0,y0));

        for (int j = 0; j < numBase; j++)
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
                    int iA, int iB, bool cols)
{
    if (cols)
    {
        GISMO_ASSERT(A.rows() == B.rows(),"Points have different dimensions\n");
        int dim = A.rows();
        gsMatrix<T> combination(dim,1);
        combination.col(0) = (1-x)*A.col(iA) + x*B.col(iB);
        return combination;
    }
    else
    {
        GISMO_ASSERT(A.cols() == B.cols(),"Points have different dimensions\n");
        int dim = A.cols();
        gsMatrix<T> combination(1,dim);
        combination.row(0) = (1-x)*A.row(iA) + x*B.row(iB);
        return combination;
    }
}

template <class T>
T distance(gsMatrix<T> const & A, int i, gsMatrix<T> const & B, int j, bool cols)
{
    T dist = 0.;

    if (cols)
    {
        GISMO_ASSERT(A.rows() == B.rows(),"Wrong matrix size\n");
        for (int d = 0; d < A.rows(); ++d)
            dist = sqrt(pow(dist,2)+pow(A(d,i)-B(d,j),2));
    }
    else
    {
        GISMO_ASSERT(A.cols() == B.cols(),"Wrong matrix size\n");
        for (int d = 0; d < A.cols(); ++d)
            dist = sqrt(pow(dist,2)+pow(A(i,d)-B(j,d),2));
    }

    return dist;
}

} // namespace ends
