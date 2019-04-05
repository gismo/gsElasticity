#include <gsCore/gsTemplateTools.h>

#include <gsElasticity/gsElMeshing.h>
#include <gsElasticity/gsElMeshing.hpp>

namespace gismo
{

//-----------------------------------//
//--------- Mesh deformation --------//
//-----------------------------------//

TEMPLATE_INST index_t computeMeshDeformation(std::vector<gsMultiPatch<real_t> > & displacements, gsMultiPatch<real_t> const & initDomain,
                                             gsBoundaryConditions<real_t> const & bdryCurves, real_t poissonRatio,
                                             index_t numSteps, index_t maxAdapt, real_t qualityRatio,
                                             bool finalize, real_t tolerance, index_t maxNumIterations);

TEMPLATE_INST void plotDeformation(std::vector<gsMultiPatch<real_t> > & displacements, gsMultiPatch<real_t> const & initDomain,
                                   std::string fileName, index_t numSamples);

TEMPLATE_INST void plotGeometry(gsMultiPatch<real_t> const & domain, std::string fileName, index_t numSamples);

TEMPLATE_INST index_t checkGeometry(gsMultiPatch<real_t> const & domain);

TEMPLATE_INST real_t geometryJacRatio(gsMultiPatch<real_t> const & domain);

//-----------------------------------//
//----------- Modelling -------------//
//-----------------------------------//

TEMPLATE_INST gsGeometry<real_t>::uPtr simplifyCurve(gsGeometry<real_t> const & curve,
                                                     index_t additionalPoints, index_t degree,
                                                     index_t numSamples);

TEMPLATE_INST real_t curveDistance(gsGeometry<real_t> const & curveA,
                                   gsGeometry<real_t> const & curveB,
                                   index_t numSamples);

TEMPLATE_INST gsGeometry<real_t>::uPtr fittingDirichlet(gsMatrix<real_t> const & params,
                                                        gsMatrix<real_t> const & points,
                                                        gsBasis<real_t> const & basis);

TEMPLATE_INST gsGeometry<real_t>::uPtr genPatchInterpolation(gsGeometry<real_t> const & A, gsGeometry<real_t> const & B,
                                                             index_t deg, index_t num, bool xiDir);

TEMPLATE_INST gsGeometry<real_t>::uPtr genPatchScaling(gsGeometry<real_t> const & boundary,
                                                                index_t deg, index_t num,
                                                                real_t scaling, gsVector<real_t> const & center);

TEMPLATE_INST gsGeometry<real_t>::uPtr genLine(index_t deg, index_t num,
                                               gsMatrix<real_t> const & A, gsMatrix<real_t> const & B,
                                               index_t iA, index_t iB);

TEMPLATE_INST gsGeometry<real_t>::uPtr genCircle(index_t deg, index_t num,
                                                 real_t radius, real_t x0, real_t y0,
                                                 real_t angle0, real_t arcAngle);

TEMPLATE_INST gsGeometry<real_t>::uPtr genCircle(gsBasis<real_t> & basis,
                                                 real_t radius, real_t x0, real_t y0,
                                                 real_t angle0, real_t arcAngle);

TEMPLATE_INST gsGeometry<real_t>::uPtr genQuad(index_t xiDeg, index_t xiNum, index_t etaDeg, index_t etaNum,
                                               gsMatrix<real_t> const & A, gsMatrix<real_t> const & B,
                                               gsMatrix<real_t> const & C, gsMatrix<real_t> const & D,
                                               index_t iA, index_t iB, index_t iC, index_t iD);

TEMPLATE_INST gsGeometry<real_t>::uPtr genSphere(index_t xiDeg, index_t xiNum, index_t etaDeg, index_t etaNum,
                                                 real_t xi0, real_t xi1,
                                                 real_t eta0, real_t eta1);

TEMPLATE_INST gsGeometry<real_t>::uPtr genSphere(gsKnotVector<real_t> & xiKnots, gsKnotVector<real_t> & etaKnots,
                                                 real_t xi0, real_t xi1,
                                                 real_t eta0, real_t eta1);

TEMPLATE_INST gsGeometry<real_t>::uPtr genCylinder(gsGeometry<real_t> const & base,
                                                   index_t deg, index_t num, real_t height);

TEMPLATE_INST gsGeometry<real_t>::uPtr genScrew(gsGeometry<real_t> const & base,
                                                index_t deg, index_t num, real_t height, real_t pitch,
                                                real_t x0, real_t y0);

//----------------------------------------//
//----------- Auxiliary functions --------//
//----------------------------------------//

TEMPLATE_INST real_t combine(real_t a, real_t b, real_t x);

TEMPLATE_INST gsMatrix<real_t> combine(gsMatrix<real_t> const & A, gsMatrix<real_t> const & B, real_t x,
                                       index_t iA, index_t iB, bool cols);

TEMPLATE_INST real_t distance(gsMatrix<real_t> const & A, index_t i, gsMatrix<real_t> const & B, index_t j, bool cols);

}
