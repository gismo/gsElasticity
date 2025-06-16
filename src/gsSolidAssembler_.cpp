#include <gsCore/gsTemplateTools.h>

#include <gsElasticity/gsSolidAssembler.h>
#include <gsElasticity/gsSolidAssembler.hpp>

#include <gsElasticity/gsLinearMaterial.h>
#include <gsElasticity/gsLinearDegradedMaterial.h>
#include <gsElasticity/gsNeoHookeLogMaterial.h>

namespace gismo
{

  CLASS_TEMPLATE_INST gsSolidAssemblerBase<real_t>;

  CLASS_TEMPLATE_INST gsSolidAssembler<2,real_t,gsLinearMaterial<real_t>>;
  CLASS_TEMPLATE_INST gsSolidAssembler<3,real_t,gsLinearMaterial<real_t>>;

  CLASS_TEMPLATE_INST gsSolidAssembler<2,real_t,gsLinearDegradedMaterial<real_t>>;
  CLASS_TEMPLATE_INST gsSolidAssembler<3,real_t,gsLinearDegradedMaterial<real_t>>;

}

