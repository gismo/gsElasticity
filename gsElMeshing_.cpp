#include <gsCore/gsTemplateTools.h>

#include <gsElasticity/gsElMeshing.h>
#include <gsElasticity/gsElMeshing.hpp>

namespace gismo
{

//-----------------------------------//
//--------- Mesh deformation --------//
//-----------------------------------//

TEMPLATE_INST void computeDeformation(std::vector<std::vector<gsMatrix<real_t> > > & deformation,
                                      gsMultiPatch<real_t> const & initDomain, gsBoundaryConditions<real_t> const & bdryCurves,
                                      index_t numSteps = 3, real_t poissonRatio = 0.49);

TEMPLATE_INST void plotDeformation(std::vector<std::vector<gsMatrix<real_t> > > const & deformation,
                                   gsMultiPatch<real_t> const & initDomain, std::string fileName,
                                   index_t numSamples = 0);

TEMPLATE_INST void applyDeformation(std::vector<std::vector<gsMatrix<real_t> > > const & deformation,
                                    gsMultiPatch<real_t> const & initDomain, gsMultiPatch<real_t> & domain);

TEMPLATE_INST void constructDisplacement(std::vector<std::vector<gsMatrix<real_t> > > const & deformation,
                                         gsMultiPatch<real_t> const & initDomain, gsMultiPatch<real_t> & displacement);

TEMPLATE_INST void computeDeformationNonlin(gsMultiPatch<real_t> & domain, gsMultiPatch<real_t> const & initDomain,
                                            gsBoundaryConditions<real_t> const & bdryCurves, gsMultiPatch<real_t> const & initGuess,
                                            index_t materialLaw = 0, real_t poissonRatio = 0.49,
                                            real_t tolerance = 1e-12, index_t maxNumIterations = 50);

TEMPLATE_INST real_t measureMinMaxJ(gsMultiPatch<real_t> const & domain, index_t numSamples = 10000);

TEMPLATE_INST void analyzeDeformation(std::vector<std::vector<gsMatrix<real_t> > > const & deformation,
                                      gsMultiPatch<real_t> const & domain, index_t measPerStep,
                                      std::string fileName, index_t numSamples = 10000);

//-----------------------------------//
//----------- Modelling -------------//
//-----------------------------------//

TEMPLATE_INST gsGeometry<real_t>::uPtr simplifyCurve(gsGeometry<real_t> const & curve,
                                                     index_t additionalPoints = 0, index_t numSamples = 1000);

TEMPLATE_INST gsGeometry<real_t>::uPtr fittingDirichlet(gsMatrix<real_t> const & params,
                                                        gsMatrix<real_t> const & points,
                                                        gsBasis<real_t> const & basis);

TEMPLATE_INST gsGeometry<real_t>::uPtr genPatchInterpolation(gsGeometry<real_t> const & A, gsGeometry<real_t> const & B,
                                                             index_t deg, index_t num, bool xiDir = false);

TEMPLATE_INST gsGeometry<real_t>::uPtr genPatchScaling(gsGeometry<real_t> const & boundary,
                                                                index_t deg, index_t num,
                                                                real_t scaling, gsVector<real_t> const & center);

TEMPLATE_INST gsGeometry<real_t>::uPtr genLine(index_t deg, index_t num,
                                               gsMatrix<real_t> const & A, gsMatrix<real_t> const & B,
                                               index_t iA = 0, index_t iB = 0);

TEMPLATE_INST gsGeometry<real_t>::uPtr genCircle(index_t deg, index_t num,
                                                 real_t radius = 1., real_t x0 = 0., real_t y0 = 0.,
                                                 real_t angle0 = 0., real_t arcAngle = 2*M_PI);

TEMPLATE_INST gsGeometry<real_t>::uPtr genCircle(gsBasis<real_t> & basis,
                                                 real_t radius = 1., real_t x0 = 0., real_t y0 = 0.,
                                                 real_t angle0 = 0., real_t arcAngle = 2*M_PI);

TEMPLATE_INST gsGeometry<real_t>::uPtr genQuad(index_t xiDeg, index_t xiNum, index_t etaDeg, index_t etaNum,
                                               gsMatrix<real_t> const & A, gsMatrix<real_t> const & B,
                                               gsMatrix<real_t> const & C, gsMatrix<real_t> const & D,
                                               index_t iA = 0, index_t iB = 0, index_t iC = 0, index_t iD = 0);

TEMPLATE_INST gsGeometry<real_t>::uPtr genSphere(index_t xiDeg, index_t xiNum, index_t etaDeg, index_t etaNum,
                                                 real_t xi0 = 0., real_t xi1 = 2*M_PI,
                                                 real_t eta0 = -M_PI/2, real_t eta1 = M_PI/2);

TEMPLATE_INST gsGeometry<real_t>::uPtr genSphere(gsKnotVector<real_t> & xiKnots, gsKnotVector<real_t> & etaKnots,
                                                 real_t xi0 = 0., real_t xi1 = 2*M_PI,
                                                 real_t eta0 = -M_PI/2, real_t eta1 = M_PI/2);

TEMPLATE_INST gsGeometry<real_t>::uPtr genCylinder(gsGeometry<real_t> const & base,
                                                   index_t deg, index_t num, real_t height);

TEMPLATE_INST gsGeometry<real_t>::uPtr genScrew(gsGeometry<real_t> const & base,
                                                index_t deg, index_t num, real_t height, real_t pitch,
                                                real_t x0 = 0., real_t y0 = 0.);

//----------------------------------------//
//----------- Auxiliary functions --------//
//----------------------------------------//

TEMPLATE_INST real_t combine(real_t a, real_t b, real_t x);

TEMPLATE_INST gsMatrix<real_t> combine(gsMatrix<real_t> const & A, gsMatrix<real_t> const & B, real_t x,
                                       index_t iA = 0, index_t iB = 0, bool cols = false);

TEMPLATE_INST real_t distance(gsMatrix<real_t> const & A, index_t i, gsMatrix<real_t> const & B, index_t j, bool cols = false );

}
