#include <gsCore/gsTemplateTools.h>

#include <gsElasticity/gsPhaseFieldAssembler.h>
#include <gsElasticity/gsPhaseFieldAssembler.hpp>

namespace gismo
{

  CLASS_TEMPLATE_INST gsPhaseFieldAssembler<real_t,PForder::Second,PFmode::AT2>;

}

