#include "casm/app/EnumeratorHandler_impl.hh"
#include "casm/clex/ScelEnum.hh"
#include "casm/clex/ConfigEnumAllOccupations.hh"
#include "casm/clex/SuperConfigEnum.hh"
#include "casm/kinetics/DiffusionTransformationEnum.hh"
#include "casm/kinetics/DiffTransConfigEnumOccPerturbations.hh"

namespace CASM {

  EnumeratorHandler::EnumeratorHandler(const ProjectSettings &set) :
    m_set(&set),
    m_enumerator(make_standard_enumerator_map()) {

    load_enumerator_plugins(
      *m_set,
      std::inserter(*m_enumerator, m_enumerator->end()),
      std::inserter(m_lib, m_lib.end()));
  }

  namespace {
    typedef std::insert_iterator<std::map<std::string, std::shared_ptr<RuntimeLibrary> > > runtimelib_it_type;
    typedef std::insert_iterator<EnumeratorMap> enum_it_type;
  }

  template std::pair<enum_it_type, runtimelib_it_type> load_enumerator_plugins(
    const ProjectSettings &set,
    enum_it_type enum_it,
    runtimelib_it_type lib_it);

}
