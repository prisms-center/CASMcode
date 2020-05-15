#ifndef CASM_HamiltonianModules
#define CASM_HamiltonianModules

#include <map>
#include <memory>
#include "casm/basis_set/DoFTraits.hh"
#include "casm/crystallography/AnisoValTraits.hh"
#include "casm/misc/cloneable_ptr.hh"
#include "casm/misc/ParsingDictionary.hh"

namespace CASM {
  template<typename T>
  class ParsingDictionary;

  class RuntimeLibrary;
  class ProjectSettings;

  class SymRepBuilderInterface;
  class AnisoValTraits;


  namespace DoFType {
    class Traits;
  }

  class HamiltonianModules: public notstd::Cloneable {
    CLONEABLE_NEEDS_DESTRUCTOR_DEF(HamiltonianModules)
  public:


    using DoFDictionary = ParsingDictionary<DoFType::Traits>;
    using SymRepBuilderDictionary = ParsingDictionary<SymRepBuilderInterface>;
    using AnisoValDictionary = ParsingDictionary<AnisoValTraits>;

    HamiltonianModules(ProjectSettings const *set = nullptr);

    DoFDictionary &dof_dict();

    DoFDictionary const &dof_dict()const;

    AnisoValDictionary &aniso_val_dict();

    AnisoValDictionary const &aniso_val_dict()const;

    SymRepBuilderDictionary &symrep_builder_dict();

    SymRepBuilderDictionary const &symrep_builder_dict() const;

  private:

    notstd::cloneable_ptr<DoFDictionary >  m_dof_dict;

    notstd::cloneable_ptr<AnisoValDictionary>  m_aniso_val_dict;

    notstd::cloneable_ptr<SymRepBuilderDictionary>  m_symrep_builder_dict;

    std::map<std::string, std::shared_ptr<RuntimeLibrary> > m_dof_lib;

    std::map<std::string, std::shared_ptr<RuntimeLibrary> > m_symrep_builder_lib;

  };

  template<>
  HamiltonianModules::AnisoValDictionary make_parsing_dictionary<AnisoValTraits>();

  template<>
  HamiltonianModules::SymRepBuilderDictionary make_parsing_dictionary<SymRepBuilderInterface>();



  /// \brief Load DoF plugins from a CASM project
  template<typename DoFDictInserter, typename RuntimeLibInserter>
  std::pair<DoFDictInserter, RuntimeLibInserter>
  load_dof_plugins(
    const ProjectSettings &set,
    DoFDictInserter dict_it,
    RuntimeLibInserter lib_it);

  /// \brief Load SymRepBuilder plugins from a CASM project
  template<typename SymRepBuilderDictInserter, typename RuntimeLibInserter>
  std::pair<SymRepBuilderDictInserter, RuntimeLibInserter>
  load_symrep_builder_plugins(
    const ProjectSettings &set,
    SymRepBuilderDictInserter dict_it,
    RuntimeLibInserter lib_it);

}

#endif
