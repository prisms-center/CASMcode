#include "casm/app/EnumeratorHandler_impl.hh"
#include "casm/app/enum/standard_enumerator_interfaces.hh"

namespace CASM {

EnumeratorHandler::EnumeratorHandler(ProjectSettings const &set)
    : m_set(&set), m_enumerator(make_standard_enumerator_interfaces()) {
  load_enumerator_plugins(*m_set, std::back_inserter(m_enumerator),
                          std::inserter(m_lib, m_lib.end()));
}

EnumeratorHandler::~EnumeratorHandler() {
  // order of deletion matters
  m_enumerator.clear();
  m_lib.clear();
}

namespace {
typedef std::insert_iterator<
    std::map<std::string, std::shared_ptr<RuntimeLibrary> > >
    runtimelib_it_type;
typedef std::insert_iterator<EnumInterfaceVector> enum_it_type;
}  // namespace

template std::pair<enum_it_type, runtimelib_it_type> load_enumerator_plugins(
    ProjectSettings const &set, enum_it_type enum_it,
    runtimelib_it_type lib_it);

}  // namespace CASM
