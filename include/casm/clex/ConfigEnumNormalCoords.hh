#ifndef CASM_ConfigEnumNormalCoords
#define CASM_ConfigEnumNormalCoords

#include "casm/symmetry/PermuteIterator.hh"
#include "casm/enumerator/InputEnumerator.hh"
#include "casm/container/Counter.hh"
#include "casm/clex/Configuration.hh"

extern "C" {
  CASM::EnumInterfaceBase *make_ConfigEnumNormalCoords_interface();
}

namespace CASM {


  /// Enumerate strained Configurations
  ///
  /// \ingroup ConfigEnumGroup
  ///
  class ConfigEnumNormalCoords : public InputEnumeratorBase<Configuration> {

    // -- Required members -------------------

  public:

    ConfigEnumNormalCoords(const ConfigEnumInput &_init,
                           DoFKey const &_dof_key,
                           Eigen::Ref<const Eigen::MatrixXd> const &_coords);


    std::string name() const override {
      return enumerator_name;

    }

    static const std::string enumerator_name;

    static std::string interface_help();

    static int run(PrimClex const &primclex,
                   jsonParser const &kwargs,
                   Completer::EnumOption const &enum_opt,
                   EnumeratorMap const *interface_map);

    static int run(PrimClex const &_primclex,
                   ConfigEnumInput const &_in_config,
                   DoFKey const &_dof,
                   std::vector<std::string> const &_filter_expr,
                   bool dry_run);
  private:
    void _set_dof();

    /// Implements increment over all strain states
    void increment() override;


    // -- Unique -------------------
    Configuration m_current;

    DoFKey m_dof_key;

    std::set<Index> m_sites;

    Eigen::MatrixXd m_coords;

    // counts over normal coords
    //Index m_i_coord;

    //set of non-equivalent transformation matrices matrices that, along with m_counter define irreducible space
    //std::vector<Eigen::MatrixXd> m_trans_mats;
    //PermuteIterator m_perm_begin, m_perm_end;

  };

}

#endif
