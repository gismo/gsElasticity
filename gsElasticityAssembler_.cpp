
#include <gsCore/gsTemplateTools.h>

#include <gsElasticity/gsElasticityFunctions.h>
#include <gsElasticity/gsElasticityFunctions.hpp>

#include <gsElasticity/gsElasticityAssembler.h>
#include <gsElasticity/gsElasticityAssembler.hpp>

namespace gismo
{
    CLASS_TEMPLATE_INST gsCauchyStressFunction<real_t>;
    CLASS_TEMPLATE_INST gsElasticityAssembler<real_t>;
}
