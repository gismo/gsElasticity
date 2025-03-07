#include <gsCore/gsTemplateTools.h>

#include <gsElasticity/gsLinearMaterial.h>
#include <gsElasticity/gsLinearDegradedMaterial.h>
#include <gsElasticity/gsNeoHookeLogMaterial.h>
#include <gsElasticity/gsNeoHookeQuadMaterial.h>

namespace gismo
{
  CLASS_TEMPLATE_INST gsLinearMaterial<real_t>;
  CLASS_TEMPLATE_INST gsLinearDegradedMaterial<real_t>;
  CLASS_TEMPLATE_INST gsNeoHookeLogMaterial<real_t>;
  CLASS_TEMPLATE_INST gsNeoHookeQuadMaterial<real_t>;

}

