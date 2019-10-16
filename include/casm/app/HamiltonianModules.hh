#ifndef CASM_HamiltonianModules
#define CASM_HamiltonianModules

#include <map>
#include <memory>

#include "casm/misc/cloneable_ptr.hh"

namespace CASM {
  template<typename T>
  class ParsingDictionary;

  class RuntimeLibrary;
  class ProjectSettings;

  namespace DoFType {
    class BasicTraits;
  }

  namespace xtal {
    namespace SpeciesAttribute_impl {
      class BasicTraits;
    }
  }


  class HamiltonianModules {
  public:

    using DoFDictionary = ParsingDictionary<DoFType::BasicTraits>;
    using SpeciesAttributeDictionary = ParsingDictionary<CASM::xtal::SpeciesAttribute_impl::BasicTraits>;

    HamiltonianModules(ProjectSettings const *set = nullptr);

    ~HamiltonianModules();

    DoFDictionary &dof_dict();

    DoFDictionary const &dof_dict()const;

    SpeciesAttributeDictionary &species_attribute_dict();

    SpeciesAttributeDictionary const &species_attribute_dict()const;

    //std::unique_ptr<HamiltonianModules > clone() const;

  private:

    notstd::cloneable_ptr<DoFDictionary >  m_dof_dict;

    notstd::cloneable_ptr<SpeciesAttributeDictionary>  m_species_attribute_dict;

    std::map<std::string, std::shared_ptr<RuntimeLibrary> > m_dof_lib;

    std::map<std::string, std::shared_ptr<RuntimeLibrary> > m_species_attribute_lib;

  };

  /// \brief Load DoF plugins from a CASM project
  template<typename DoFDictInserter, typename RuntimeLibInserter>
  std::pair<DoFDictInserter, RuntimeLibInserter>
  load_dof_plugins(
    const ProjectSettings &set,
    DoFDictInserter dict_it,
    RuntimeLibInserter lib_it);

  /// \brief Load SpeciesAttribute plugins from a CASM project
  template<typename AttributeDictInserter, typename RuntimeLibInserter>
  std::pair<AttributeDictInserter, RuntimeLibInserter>
  load_species_attribute_plugins(
    const ProjectSettings &set,
    AttributeDictInserter dict_it,
    RuntimeLibInserter lib_it);

}

#endif
