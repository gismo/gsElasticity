#include <gsCore/gsTemplateTools.h>

#include <gsElasticity/gsGeoUtils.h>
#include <gsElasticity/gsGeoUtils.hpp>

namespace gismo
{

//-----------------------------------//
//--------- Mesh Analysis -----------//
//-----------------------------------//

TEMPLATE_INST void plotGeometry(gsMultiPatch<real_t> const & domain, std::string fileName, index_t numSamples);

TEMPLATE_INST void plotGeometry(const gsMultiPatch<real_t> & domain,std::string const & fileName,
                                gsParaviewCollection & collection, index_t step);

TEMPLATE_INST void plotDeformation(const gsMultiPatch<real_t> & initDomain, const std::vector<gsMultiPatch<real_t> > & displacements,
                     std::string fileName, index_t numSamplingPoints);

TEMPLATE_INST void plotDeformation(const gsMultiPatch<real_t> & initDomain, const gsMultiPatch<real_t> & displacement,
                                   std::string const & fileName, gsParaviewCollection & collection, index_t step);

TEMPLATE_INST index_t checkGeometry(gsMultiPatch<real_t> const & domain);

TEMPLATE_INST index_t checkDisplacement(gsMultiPatch<real_t> const & domain, gsMultiPatch<real_t> const & displacement);

TEMPLATE_INST real_t normL2(gsMultiPatch<real_t> const & domain, gsMultiPatch<real_t> const & solution);


TEMPLATE_INST real_t geometryJacRatio(gsMultiPatch<real_t> const & domain);

TEMPLATE_INST real_t displacementJacRatio(gsMultiPatch<real_t> const & domain, gsMultiPatch<real_t> const & displacement);

TEMPLATE_INST void genSamplingPoints(const gsVector<real_t> & lower, const gsVector<real_t> & upper,
                                     const gsQuadRule<real_t> & quRule, gsMatrix<real_t> & points);

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

TEMPLATE_INST gsGeometry<real_t>::uPtr genSpring(real_t springRadius, real_t springPitch,
                                                 real_t wireRadius, index_t numQuarterSegments);

//----------------------------------------//
//----------- Auxiliary functions --------//
//----------------------------------------//

TEMPLATE_INST real_t combine(real_t a, real_t b, real_t x);

TEMPLATE_INST gsMatrix<real_t> combine(gsMatrix<real_t> const & A, gsMatrix<real_t> const & B, real_t x,
                                       index_t iA, index_t iB, bool cols);

TEMPLATE_INST real_t distance(gsMatrix<real_t> const & A, index_t i, gsMatrix<real_t> const & B, index_t j, bool cols);

}
