/** @file gsElMeshing.h

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

//-----------------------------------//
//--------- Mesh deformation --------//
//-----------------------------------//

/// @brief Computes a deformation of a given initial domain using linear elasticity with incremental Dirichlet BC
/// such that the domain's boundary coincides with a given set of boundary curves;
/// returns a set of displacement control points for every incremental steps and every patch
template <class T>
void computeDeformation(std::vector<std::vector<gsMatrix<T> > > & deformation,
                        gsMultiPatch<T> const & initDomain, gsBoundaryConditions<T> const & bdryCurves,
                        index_t numSteps = 3, T poissonRatio = 0.49);

/// @brief Plots steps of the computed incremental deformation; outputs a deformed mesh for every incremental step;
/// set <numSamples> greater than 0 to plot the Jacobian of the geometry mapping using <numSamples> points
template <class T>
void plotDeformation(std::vector<std::vector<gsMatrix<T> > > const & deformation,
                     gsMultiPatch<T> const & intiDomain, std::string fileName,
                     index_t numSamples = 0);

/// @brief Deformes the initial domain using a previously computed incremental <deformation>
template <class T>
void applyDeformation(std::vector<std::vector<gsMatrix<T> > > const & deformation,
                      gsMultiPatch<T> const & initDomain, gsMultiPatch<T> & domain);

/// @brief Construct an isogeometric displacement field
template <class T>
void constructDisplacement(std::vector<std::vector<gsMatrix<T> > > const & deformation,
                           gsMultiPatch<T> const & initDomain, gsMultiPatch<T> & displacement);

/// @brief Computes a deformation of a given initial domain using nonlinear elasticity and a given initial guess;
/// such that the domain's boundary coincides with a given set of boundary curves;
/// returns a deformed/parametrized domain
template <class T>
void computeDeformationNonlin(gsMultiPatch<T> & domain, gsMultiPatch<T> const & initDomain,
                              gsBoundaryConditions<T> const & bdryCurves, gsMultiPatch<T> const & initGuess,
                              index_t materialLaw = 0, T poissonRatio = 0.49,
                              T tolerance = 1e-12, index_t maxNumIterations = 50);

/// @brief Checks whether configuration is bijective, i.e. det(Jac) > 0;
/// return -1 if yes or a number of the first invalid patch
template <class T>
index_t checkGeometry(gsMultiPatch<T> const & domain);

/// @brief Plots the mesh and the jacobian (if <numSamples> > 0) to Paraview
template <class T>
void plotGeometry(gsMultiPatch<T> const & domain, std::string fileName, index_t numSamples);

/// @brief Returns max(Jacobian determinant) divided by min(Jacobian determinant)
template <class T>
T measureMinMaxJ(gsMultiPatch<T> const & domain);

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
