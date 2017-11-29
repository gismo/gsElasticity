#include <gsCore/gsTemplateTools.h>
#include <gsThermoElasticity/gsElPoissonAssembler.h>
#include <gsThermoElasticity/gsElPoissonAssembler.hpp>
#include <gsThermoElasticity/gsElThermoAssembler.h>
#include <gsThermoElasticity/gsElThermoAssembler.hpp>

namespace gismo
{
    CLASS_TEMPLATE_INST gsElPoissonAssembler<real_t>;
    CLASS_TEMPLATE_INST gsElThermoAssembler<real_t>;
}

