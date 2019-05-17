#include <gsCore/gsTemplateTools.h>

#include <gsElasticity/gsElasticityNewton.h>
#include <gsElasticity/gsElasticityNewton.hpp>

#include <gsElasticity/gsElasticityNewtonDeLuxe.h>
#include <gsElasticity/gsElasticityNewtonDeLuxe.hpp>

namespace gismo
{
    CLASS_TEMPLATE_INST gsElasticityNewton<real_t>;
    CLASS_TEMPLATE_INST gsElasticityNewtonDeLuxe<real_t>;
}
