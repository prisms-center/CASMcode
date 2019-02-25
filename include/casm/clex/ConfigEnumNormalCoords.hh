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

    ConfigEnumNormalCoords(ConfigEnumInput const &_init,
                           DoFKey const &_dof,
                           Eigen::Ref<const Eigen::MatrixXd> const &_axes,
                           Eigen::Ref<const Eigen::VectorXd> const &min_val,
                           Eigen::Ref<const Eigen::VectorXd> const &max_val,
                           Eigen::Ref<const Eigen::VectorXd> const &inc_val,
                           Index _min_nonzero,
                           Index _max_nonzero);



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
                   Eigen::Ref<const Eigen::MatrixXd> const &_axes,
                   Eigen::Ref<const Eigen::VectorXd> const &min_val,
                   Eigen::Ref<const Eigen::VectorXd> const &max_val,
                   Eigen::Ref<const Eigen::VectorXd> const &inc_val,
                   bool sym_axes,
                   Index _min_nonzero,
                   Index _max_nonzero,
                   std::vector<std::string> const &_filter_expr,
                   bool dry_run);

    /// Implements increment over all strain states
    void increment() override;

  private:
    void _set_dof();

    void _initialize_count();

    bool _increment_combo();

    bool _check_sparsity() const;

    bool _check_current() const;

    // -- Unique -------------------
    notstd::cloneable_ptr<Configuration> m_current;

    ConfigDoF::LocalDoFContainerType *m_dof_vals;

    DoFKey m_dof_key;

    Index m_min_nonzero;

    Index m_max_nonzero;

    Eigen::MatrixXd m_axes;

    Eigen::VectorXd m_min;

    Eigen::VectorXd m_max;

    Eigen::VectorXd m_inc;

    bool m_unit_length;

    std::vector<Index> m_sites;

    std::vector<Index> m_dof_dims;

    bool m_subset_mode;

    bool m_combo_mode;

    Index m_combo_index;

    std::vector<Index> m_combo;

    EigenCounter<Eigen::VectorXd> m_counter;

    // counts over normal coords
    //Index m_i_coord;

    //set of non-equivalent transformation matrices matrices that, along with m_counter define irreducible space
    //std::vector<Eigen::MatrixXd> m_trans_mats;
    //PermuteIterator m_perm_begin, m_perm_end;

  };

}

#endif
