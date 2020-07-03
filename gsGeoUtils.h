/** @file gsGeoUtils.h

    @brief Provides isogeometric meshing and modelling routines.

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s):
        A.Shamanskiy (2016 - ...., TU Kaiserslautern)
*/

#pragma once

#include <gsCore/gsMultiPatch.h>

namespace gismo
{

class gsParaviewCollection;

//-----------------------------------//
//--------- Mesh Analysis -----------//
//-----------------------------------//

/// @brief Plots the mesh and the jacobian (if <numSamples> > 0) to Paraview
template <class T>
void plotGeometry(gsMultiPatch<T> const & domain, std::string fileName, index_t numSamples);

/// plot an isogeometric mesh and add to collection
template <class T>
void plotGeometry(const gsMultiPatch<T> & domain, std::string const & fileName,
                  gsParaviewCollection & collection, index_t step);

/// use all saved displacement fields to plot the all intermediate deformed configurations of the computational domain;
/// always plots the deformed isoparametric mesh; plots the Jacobian determinant of the deformed configuration if *numSamplingPoints* > 0
template <class T>
void plotDeformation(const gsMultiPatch<T> & initDomain, const std::vector<gsMultiPatch<T> > & displacements,
                     std::string fileName, index_t numSamplingPoints = 10000);

/// plot a deformed isogeometric mesh and add it to a Paraview collection
template <class T>
void plotDeformation(const gsMultiPatch<T> & initDomain, const gsMultiPatch<T> & displacement,
                     std::string const & fileName, gsParaviewCollection & collection, index_t step);

/// @brief Checks whether configuration is bijective, i.e. det(Jac(geo)) > 0;
/// returns -1 if yes or the number of the first invalid patch;
/// samples the Jacobian elementwise at the quadrature points and the corners
template <class T>
index_t checkGeometry(gsMultiPatch<T> const & domain);

/// @brief Checks whether the deformed configuration is bijective, i.e. det(Jac(geo+disp)) > 0;
/// returns -1 if yes or the number of the first invalid patch;
/// samples the Jacobian elementwise at the quadrature points and the corners
template <class T>
index_t checkDisplacement(gsMultiPatch<T> const & domain, gsMultiPatch<T> const & displacement);

/// @ Compute norm of the isogeometric solution
template <class T>
T normL2(gsMultiPatch<T> const & domain, gsMultiPatch<T> const & solution);

/// @brief Returns min(Jacobian determinant) divided by max(Jacobian determinant);
/// samples the Jacobian elementwise at the quadrature points and the corners
template <class T>
T geometryJacRatio(gsMultiPatch<T> const & domain);

/// @brief Returns min(Jacobian determinant) divided by max(Jacobian determinant) for geo+disp
/// samples the Jacobian elementwise at the quadrature points and the corners
template <class T>
T displacementJacRatio(gsMultiPatch<T> const & domain, gsMultiPatch<T> const & displacement);

/// @brief Generates a matrix of sampling points for a given parametric element;
/// includes quadrature points for the element as well as the corner points
template <class T>
void genSamplingPoints(const gsVector<T> & lower, const gsVector<T> & upper,
                       const gsQuadRule<T> & quRule, gsMatrix<T> & points);

/// @brief compute length of a patch in a given parametric direction as a mean of all boundary edges corresponding to this direction
template <class T>
T patchLength(const gsGeometry<T> & geo, short_t dir = 0);

/// @brief compute curve length
template <class T>
T curveLength(const gsGeometry<T> & geo);

/// @brief distributes sampling points according to the length of the patch in each parametric direction
template <class T>
gsVector<unsigned> distributePoints(const gsGeometry<T> & geo, unsigned numPoints);

//-----------------------------------//
//----------- Modelling--------------//
//-----------------------------------//

/// @brief generates a simplified curve by fitting with the coarsest basis of the same degree;
/// then reparametrizes it using the basis of the original curve;
/// <additionalPoints> increases the number of degrees of freedom;
/// <numSamples> is a number of sampling points for reparametrization
template<class T>
typename gsGeometry<T>::uPtr simplifyCurve(gsGeometry<T> const & curve,
                                          index_t additionalPoints = 0, index_t degree = 0,
                                          index_t numSamples = 1000);

/// @brief returns a distance in L2 sense between two curves parametrized from 0 to 1
template<class T>
T curveDistance(gsGeometry<T> const & curveA,
                gsGeometry<T> const & curveB,
                index_t numSamples = 1000);

/// @brief fits a given parametrized point cloud with a curve using a given basis;
/// the resulting curve interpolates the first and the last points
template <class T>
typename gsGeometry<T>::uPtr fittingDirichlet(gsMatrix<T> const & params,
                                              gsMatrix<T> const & points,
                                              gsBasis<T> const & basis);

/// @grief generates a tensor product B-spline patch by interpolating between the two given B-spline patches
/// of a dimension one lower;
/// in 2D case, interpolates in eta/south-north direction by default;
/// in 3D case, always interpolates in the third direction/front-back
template<class T>
typename gsGeometry<T>::uPtr genPatchInterpolation(gsGeometry<T> const & A, gsGeometry<T> const & B,
                                                   index_t deg, index_t num, bool xiDir = false);

/// @brief generates a tensor product B-spline    bdry   south | front
/// patch by scaling a given geometry object      \  /     |   |   |
/// towards a given center point;                 (x,y)  north |  back
/// oppositely lying bdry is generated by scaling the original boundary with <scaling> coeff
template<class T>
typename gsGeometry<T>::uPtr genPatchScaling(gsGeometry<T> const & boundary,
                                             index_t deg, index_t num,
                                             T scaling, gsVector<T> const & center);

/// @breif generates a uniformly parametrized straight line with a B-spline basis;
/// end points are given as ROWS of matrices
template<class T>
typename gsGeometry<T>::uPtr genLine(index_t deg, index_t num,
                                     gsMatrix<T> const & A, gsMatrix<T> const & B,
                                     index_t iA = 0, index_t iB = 0);

/// @breif generates a uniformly parametrized circular arc with a B-spline basis;
/// generates a closed circle if <arcAngle> is 2pi
template<class T>
typename gsGeometry<T>::uPtr genCircle(index_t deg, index_t num,
                                       T radius = 1., T x0 = 0., T y0 = 0.,
                                       T angle0 = 0., T arcAngle = 2*M_PI);

/// @breif generates a uniformly parametrized circular arc with a given basis;
/// generates a closed circle if <arcAngle> is 2pi
template<class T>
typename gsGeometry<T>::uPtr genCircle(gsBasis<T> & basis,
                                       T radius = 1., T x0 = 0., T y0 = 0.,
                                       T angle0 = 0., T arcAngle = 2*M_PI);

/// @brief generates a quad patch given by its four   C---D
/// corners with the following orientation;           |   |
/// the points are given as ROWS of matrices          A---B
template<class T>
typename gsGeometry<T>::uPtr genQuad(index_t xiDeg, index_t xiNum, index_t etaDeg, index_t etaNum,
                                     gsMatrix<T> const & A, gsMatrix<T> const & B,
                                     gsMatrix<T> const & C, gsMatrix<T> const & D,
                                     index_t iA = 0, index_t iB = 0, index_t iC = 0, index_t iD = 0);

/// @brief generates a tensor product B-spline spherical patch with radius 1 and center at 0
/// given the degrees and number of control points in two parametric dimensions
template<class T>
typename gsGeometry<T>::uPtr genSphere(index_t xiDeg, index_t xiNum, index_t etaDeg, index_t etaNum,
                                       T xi0 = 0., T xi1 = 2*M_PI,
                                       T eta0 = -M_PI/2, T eta1 = M_PI/2);

/// @brief generates a tensor product B-spline spherical patch with radius 1 and center at 0
/// given knot vectors for two parametric dimensions
template<class T>
typename gsGeometry<T>::uPtr genSphere(gsKnotVector<T> & xiKnots, gsKnotVector<T> & etaKnots,
                                       T xi0 = 0., T xi1 = 2*M_PI,
                                       T eta0 = -M_PI/2, T eta1 = M_PI/2);

/// @brief generates a 3D tensor product B-spline cylindrical patch
template<class T>
typename gsGeometry<T>::uPtr genCylinder(gsGeometry<T> const & base,
                                         index_t deg, index_t num, T height);

/// @brief generates a 3D tensor product B-spline screw-like patch
template<class T>
typename gsGeometry<T>::uPtr genScrew(gsGeometry<T> const & base,
                                      index_t deg, index_t num,
                                      T height, T pitch, T x0 = 0., T y0 = 0.);

/// @brief generates a 3D NURBS spring using provided geometry as a cross-section
template<class T>
typename gsGeometry<T>::uPtr genSpring(gsGeometry<T> const & crossSection,
                                       T springRadius = 6.0, T springPitch = 2.60258,
                                       index_t numQuarterSegments = 12, bool nurbs = false);

/// @brief This is more of a script than a function. I use it to generate a multi-parametrization for the biceps model given its surface
template<class T>
void genMuscleMP(gsGeometry<T> const & muscleSurface, gsMultiPatch<T> & result);

//----------------------------------------//
//----------- Auxiliary functions --------//
//----------------------------------------//

/// @brief compute a convex combintation (1-x)a+xb
template<class T>
inline T combine(T a, T b, T x) { return a*(1-x)+b*x; }

/// @brief compute a convex combination of two points given as ROWS
/// of matrices numbered <iA> and <iB>;
/// set <cols> to <true> to give points as COLUMNS
template<class T>
gsMatrix<T> combine(gsMatrix<T> const & A, gsMatrix<T> const & B, T x,
                    index_t iA = 0, index_t iB = 0, bool cols = false);

/// @brief compute a distance between the point number <i> in the set <A>
/// and the point number <j> in the set <B>;
/// by default, points are given as ROWS of matrices;
/// set <cols> to <true> to give points as COLUMNS
template <class T>
T distance(gsMatrix<T> const & A, gsMatrix<T> const & B,
           index_t i = 0, index_t j = 0, bool cols = false);

} // namespace ends

