#include <gsCore/gsTemplateTools.h>

#include <gsElasticity/gsElMeshing.h>
#include <gsElasticity/gsElMeshing.hpp>

namespace gismo
{

//-----------------------------------//
//--------- Mesh deformation --------//
//-----------------------------------//

TEMPLATE_INST void computeDeformationInc(std::vector<std::vector<gsMatrix<real_t> > > & result,
                          gsMultiPatch<real_t> const & domain, gsBoundaryConditions<real_t> const & bdry,
                          int numSteps = 3, real_t poissonRatio = 0.49);

TEMPLATE_INST void deformNonLinearly(gsMultiPatch<real_t> & result, gsMultiPatch<real_t> const & domain,
                                     std::vector<std::vector<gsMatrix<real_t> > > const & deformation,
                                     gsMatrix<real_t> const & sumSolutionVector,
                                     gsBoundaryConditions<real_t> const & bdry, real_t poissonRatio = 0.49,
                                     real_t tolerance = 1e-12, int maxNumIterations = 20);

TEMPLATE_INST void plotDeformation(std::vector<std::vector<gsMatrix<real_t> > > const & deformation,
                                   gsMultiPatch<real_t> const & domain, std::string fileName,
                                   bool plotJac = false, int numSamples = 1);

TEMPLATE_INST void applyDeformation(std::vector<std::vector<gsMatrix<real_t> > > const & deformation,
                                    gsMultiPatch<real_t> const & initDomain, gsMultiPatch<real_t> & domain);

TEMPLATE_INST real_t measureMinMaxJ(gsMultiPatch<real_t> const & domain, int numSamples = 10000);

TEMPLATE_INST void analyzeDeformation(std::vector<std::vector<gsMatrix<real_t> > > const & deformation,
                                      gsMultiPatch<real_t> const & domain, int measPerStep,
                                      std::string fileName, int numSamples = 10000);

//-----------------------------------//
//----------- Modelling -------------//
//-----------------------------------//

TEMPLATE_INST gsGeometry<real_t>::uPtr simplifyCurve(gsGeometry<real_t> const & curve,
                                                     int additionalPoints = 0, int numSamples = 1000);

TEMPLATE_INST gsGeometry<real_t>::uPtr fittingDirichlet(gsMatrix<real_t> const & params,
                                                        gsMatrix<real_t> const & points,
                                                        gsBasis<real_t> const & basis);

TEMPLATE_INST gsGeometry<real_t>::uPtr genPatchInterpolation(gsGeometry<real_t> const & A, gsGeometry<real_t> const & B,
                                                             int deg, int num, bool xiDir = false);

TEMPLATE_INST gsGeometry<real_t>::uPtr genPatchScaling(gsGeometry<real_t> const & boundary,
                                                                int deg, int num,
                                                                real_t scaling, gsVector<real_t> const & center);

TEMPLATE_INST gsGeometry<real_t>::uPtr genLine(int deg, int num,
                                               gsMatrix<real_t> const & A, gsMatrix<real_t> const & B,
                                               int iA = 0, int iB = 0);

TEMPLATE_INST gsGeometry<real_t>::uPtr genCircle(int deg, int num,
                                                 real_t radius = 1., real_t x0 = 0., real_t y0 = 0.,
                                                 real_t angle0 = 0., real_t arcAngle = 2*M_PI);

TEMPLATE_INST gsGeometry<real_t>::uPtr genCircle(gsBasis<real_t> & basis,
                                                 real_t radius = 1., real_t x0 = 0., real_t y0 = 0.,
                                                 real_t angle0 = 0., real_t arcAngle = 2*M_PI);

TEMPLATE_INST gsGeometry<real_t>::uPtr genQuad(int xiDeg, int xiNum, int etaDeg, int etaNum,
                                               gsMatrix<real_t> const & A, gsMatrix<real_t> const & B,
                                               gsMatrix<real_t> const & C, gsMatrix<real_t> const & D,
                                               int iA = 0, int iB = 0, int iC = 0, int iD = 0);

TEMPLATE_INST gsGeometry<real_t>::uPtr genSphere(int xiDeg, int xiNum, int etaDeg, int etaNum,
                                                 real_t xi0 = 0., real_t xi1 = 2*M_PI,
                                                 real_t eta0 = -M_PI/2, real_t eta1 = M_PI/2);

TEMPLATE_INST gsGeometry<real_t>::uPtr genSphere(gsKnotVector<real_t> & xiKnots, gsKnotVector<real_t> & etaKnots,
                                                 real_t xi0 = 0., real_t xi1 = 2*M_PI,
                                                 real_t eta0 = -M_PI/2, real_t eta1 = M_PI/2);

TEMPLATE_INST gsGeometry<real_t>::uPtr genCylinder(gsGeometry<real_t> const & base,
                                                   int deg, int num, real_t height);

TEMPLATE_INST gsGeometry<real_t>::uPtr genScrew(gsGeometry<real_t> const & base,
                                                int deg, int num, real_t height, real_t pitch,
                                                real_t x0 = 0., real_t y0 = 0.);

//----------------------------------------//
//----------- Auxiliary functions --------//
//----------------------------------------//

TEMPLATE_INST real_t combine(real_t a, real_t b, real_t x);

TEMPLATE_INST gsMatrix<real_t> combine(gsMatrix<real_t> const & A, gsMatrix<real_t> const & B, real_t x,
                                       int iA = 0, int iB = 0, bool cols = false);

TEMPLATE_INST real_t distance(gsMatrix<real_t> const & A, int i, gsMatrix<real_t> const & B, int j, bool cols = false );

}
