/** @file gsCouplingUtilities.hpp

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

   Author(s): A. Shamanskiy
*/

#include <gsThermoElasticity/gsCouplingUtilities.h>

namespace gismo
{

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
gsFluxFunctionExternal<T>::gsFluxFunctionExternal(gsMatrix<T> const & points, gsMatrix<T> const & normGrads, T conductivity)
    : m_points(points),
      m_grads(normGrads),
      m_cond(conductivity)
{

}

template <class T>
void gsFluxFunctionExternal<T>::eval_into(gsMatrix<T> const & u, gsMatrix<T> & result) const
{
    int numP = u.cols();

    result.resize(1,numP);
    T ratio;
    int index;

    gsMatrix<T> point(1,2);

    for (int i = 0; i < numP; ++i)
    {
        point(0,0) = u(0,i);
        point(0,1) = u(1,i);
        gsProjectPiecewiseLin(point,m_points,ratio,index);
        result.at(i) = m_cond*(m_grads.at(index)*(1.-ratio) + m_grads.at(index+1)*ratio);
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
    gsMatrix<T> grads;
    m_geo.patch(m_patch).deriv_into(u,grads);

    for (int i = 0; i < u.cols(); ++i)
    {
        result(0,i) = geoEval->jacDet(i);
    }

}

template <class T>
void gsExtractNormalGradients(gsMultiPatch<T> const & sourceGeo, gsMultiPatch<T> const & sourceSol,
                              std::vector<int> patches, std::vector<boundary::side> sides,
                              gsMatrix<T> & points, gsMatrix<T> & normGrads)
{
    points.resize(0,2);
    normGrads.resize(0,1);

    for (unsigned p = 0; p < patches.size(); ++p)
    {
        gsMatrix<T> coefs = sourceGeo.patch(patches[p]).boundary(sides[p])->coefs();
        int numPoints = coefs.rows();
        gsFluxFunction<T> fluxFunc(patches[p],sides[p],sourceGeo,sourceSol);
        gsMatrix<real_t> fluxes;
        fluxFunc.eval_into(coefs.transpose(),fluxes);

        if (p==0)
        {
            points.conservativeResize(numPoints,2);
            normGrads.conservativeResize(numPoints,1);
            points.block(0,0,numPoints,2) = coefs.block(0,0,numPoints,2);
            normGrads.block(0,0,numPoints,1) = fluxes.transpose().block(0,0,numPoints,1);
        }
        else
        {
            if (abs(coefs(0,0)- points(points.cols()-1,0)) >1e-9 || abs(coefs(0,1) - points(points.cols()-1,1))>1e-9)
            {
                gsMatrix<T> newCoefs = coefs;
                gsMatrix<T> newFluxes = fluxes;
                for (int i = 0; i < numPoints; ++i)
                {
                    coefs(i,0) = newCoefs(numPoints-i-1,0);
                    coefs(i,1) = newCoefs(numPoints-i-1,1);

                    fluxes(0,i) = newFluxes(0,numPoints-i-1);
                }
            }

            normGrads(normGrads.rows()-1,0) = (normGrads(normGrads.rows()-1,0) + fluxes.at(0))/2.;
            points.conservativeResize(points.rows()+numPoints-1,2);
            normGrads.conservativeResize(normGrads.rows()+numPoints-1,1);
            points.block(points.rows()-numPoints+1,0,numPoints-1,2) = coefs.block(1,0,numPoints-1,2);
            normGrads.block(normGrads.rows()-numPoints+1,0,numPoints-1,1) = fluxes.transpose().block(1,0,numPoints-1,1);
        }
    }
}

template <class T>
void gsProjectPiecewiseLin(gsMatrix<T> const & point, gsMatrix<T> const & set, T & ratio, int & index)
{
    T dist;
    index = 0;
    gsMatrix<T> A(1,2);
    A.at(0) = set(0,0);
    A.at(1) = set(0,1);
    gsMatrix<T> B(1,2);
    B.at(0) = set(1,0);
    B.at(1) = set(1,1);

    gsDistPointSegment(A,B,point,ratio,dist);

    for (int i = 1; i < set.rows()-1; ++i)
    {
        T newDist;
        T newRatio;
        A.at(0) = set(i,0);
        A.at(1) = set(i,1);
        B.at(0) = set(i+1,0);
        B.at(1) = set(i+1,1);

        gsDistPointSegment(A,B,point,newRatio,newDist);
        if (newDist < dist)
        {
            dist = newDist;
            ratio = newRatio;
            index = i;
        }
    }
}

template <class T>
void gsDistPointSegment(gsMatrix<T> const & A, gsMatrix<T> const & B, gsMatrix<T> const & C,
                        T & dist, T & ratio)
{
    T ABxAC = (B.at(0)-A.at(0))*(C.at(0)-A.at(0)) + (B.at(1)-A.at(1))*(C.at(1)-A.at(1));
    T BAxBC = (A.at(0)-B.at(0))*(C.at(0)-B.at(0)) + (A.at(1)-B.at(1))*(C.at(1)-B.at(1));
    if (ABxAC < 0)
    {
        ratio = 0.;
        dist = sqrt(pow(C.at(0)-A.at(0),2) + pow(C.at(1)-A.at(1),2));
    }
    else if (BAxBC < 0)
    {
        ratio = 1.;
        dist = sqrt(pow(C.at(0)-B.at(0),2) + pow(C.at(1)-B.at(1),2));
    }
    else
    {
        T AC = sqrt(pow(C.at(0)-A.at(0),2) + pow(C.at(1)-A.at(1),2));
        T AB = sqrt(pow(B.at(0)-A.at(0),2) + pow(B.at(1)-A.at(1),2));
        T cosA = ABxAC / AC / AB;
        T AO = AC*cosA;

        ratio = AO/AB;
        dist = sqrt(pow(AC,2) - pow(AO,2));
    }
}

} // namespace ends
