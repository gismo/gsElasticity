#include <gsCore/gsTemplateTools.h>

#include <gsElasticity/gsPhaseFieldAssembler.h>
#include <gsElasticity/gsPhaseFieldAssembler.hpp>

namespace gismo
{

  CLASS_TEMPLATE_INST gsPhaseFieldAssemblerBase<real_t>;
  CLASS_TEMPLATE_INST gsPhaseFieldAssembler<real_t,PForder::Second,PFmode::AT2>;
  CLASS_TEMPLATE_INST gsPhaseFieldAssembler<real_t,PForder::Fourth,PFmode::AT2>;

}

