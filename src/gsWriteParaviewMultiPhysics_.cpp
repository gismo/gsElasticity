
#include <gsCore/gsTemplateTools.h>

#include <gsElasticity/src/gsWriteParaviewMultiPhysics.h>
#include <gsElasticity/src/gsWriteParaviewMultiPhysics.hpp>

namespace gismo
{
TEMPLATE_INST
void gsWriteParaviewMultiPhysics(std::map<std::string, const gsField<real_t>* > fields, std::string const & fn,
                     unsigned npts, bool mesh, bool cnet);

TEMPLATE_INST
void gsWriteParaviewMultiPhysicsTimeStep(std::map<std::string, const gsField<real_t> *> fields, std::string const & fn,
                                         gsParaviewCollection & collection, int time, unsigned npts);

TEMPLATE_INST
void gsWriteParaviewMultiPhysicsSinglePatch(std::map<std::string,const gsField<real_t>* > fields,
                                const unsigned patchNum,
                                std::string const & fn,
                                unsigned npts);

TEMPLATE_INST
void gsWriteParaviewMultiTPgrid(gsMatrix<real_t> const& points,
                                std::map<std::string, gsMatrix<real_t> >& data,
                                const gsVector<index_t> & np,
                                std::string const & fn);
}
