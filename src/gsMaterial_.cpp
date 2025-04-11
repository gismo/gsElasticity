#include <gsCore/gsTemplateTools.h>

#include <gsElasticity/src/gsLinearMaterial.h>
#include <gsElasticity/src/gsNeoHookeLogMaterial.h>
#include <gsElasticity/src/gsNeoHookeQuadMaterial.h>

namespace gismo
{
  CLASS_TEMPLATE_INST gsLinearMaterial<real_t>;
  CLASS_TEMPLATE_INST gsNeoHookeLogMaterial<real_t>;
  CLASS_TEMPLATE_INST gsNeoHookeQuadMaterial<real_t>;

}

