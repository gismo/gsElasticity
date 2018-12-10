#include <gsCore/gsTemplateTools.h>
#include <gsThermoElasticity/gsCouplingUtilities.h>
#include <gsThermoElasticity/gsCouplingUtilities.hpp>
#include <gsThermoElasticity/gsElPoissonAssembler.h>
#include <gsThermoElasticity/gsElPoissonAssembler.hpp>
#include <gsThermoElasticity/gsElThermoAssembler.h>
#include <gsThermoElasticity/gsElThermoAssembler.hpp>

namespace gismo
{
    CLASS_TEMPLATE_INST gsElPoissonAssembler<real_t>;
    CLASS_TEMPLATE_INST gsElThermoAssembler<real_t>;
    CLASS_TEMPLATE_INST gsFluxFunction<real_t>;
    CLASS_TEMPLATE_INST gsFluxFunctionExternal<real_t>;
    CLASS_TEMPLATE_INST gsGradFunction<real_t>;
    CLASS_TEMPLATE_INST gsPhysGradFunction<real_t>;
    CLASS_TEMPLATE_INST gsDetFunction<real_t>;

    TEMPLATE_INST void gsExtractNormalGradients(gsMultiPatch<real_t> const & sourceGeo, gsMultiPatch<real_t> const & sourceSol,
                                                std::vector<int> patches, std::vector<boundary::side> sides,
                                                gsMatrix<real_t> & points, gsMatrix<real_t> & normGrads);

    TEMPLATE_INST void gsProjectPiecewiseLin(gsMatrix<real_t> const & point, gsMatrix<real_t> const & set, real_t & ratio, int & index);

    TEMPLATE_INST void gsDistPointSegment(gsMatrix<real_t> const & A, gsMatrix<real_t> const & B, gsMatrix<real_t> const & C,
                                          real_t & dist, real_t & ratio);
}

