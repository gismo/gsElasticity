
#include <gsCore/gsTemplateTools.h>

#include <gsElasticity/gsElasticityAssembler.h>
#include <gsElasticity/gsElasticityAssembler.hpp>

#include <gsElasticity/gsElasticityMassAssembler.h>
#include <gsElasticity/gsElasticityMassAssembler.hpp>

#include <gsElasticity/gsElasticityNewton.h>

namespace gismo
{
    CLASS_TEMPLATE_INST gsElasticityAssembler<real_t>;
	CLASS_TEMPLATE_INST gsElasticityMassAssembler<real_t>;
}
