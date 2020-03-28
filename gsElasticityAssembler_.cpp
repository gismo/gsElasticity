
#include <gsCore/gsTemplateTools.h>

#include <gsElasticity/gsElasticityFunctions.h>
#include <gsElasticity/gsElasticityFunctions.hpp>

#include <gsElasticity/gsElasticityAssembler.h>
#include <gsElasticity/gsElasticityAssembler.hpp>

namespace gismo
{
    CLASS_TEMPLATE_INST gsCauchyStressFunction<real_t>;
    CLASS_TEMPLATE_INST gsDetFunction<real_t>;
    CLASS_TEMPLATE_INST gsStiffFunction<real_t>;
    CLASS_TEMPLATE_INST gsFsiLoad<real_t>;
    CLASS_TEMPLATE_INST gsElasticityAssembler<real_t>;
}
