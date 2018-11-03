#ifndef CASM_ConfigEnumStrain
#define CASM_ConfigEnumStrain

#include "casm/symmetry/PermuteIterator.hh"
#include "casm/strain/StrainConverter.hh"
#include "casm/enumerator/InputEnumerator.hh"
#include "casm/container/Counter.hh"
#include "casm/clex/Configuration.hh"

extern "C" {
  CASM::EnumInterfaceBase *make_ConfigEnumStrain_interface();
}

namespace CASM {

  namespace SymRepTools {
    class SubWedge;
  }

  /// Enumerate strained Configurations
  ///
  /// \ingroup ConfigEnumGroup
  ///
  class ConfigEnumStrain : public InputEnumeratorBase<Configuration> {

    // -- Required members -------------------

  public:

    ConfigEnumStrain(const Configuration &_init,
                     const std::vector<SymRepTools::SubWedge> &_wedges,
                     Eigen::VectorXd min_val,
                     Eigen::VectorXd max_val,
                     Eigen::VectorXd inc_val,
                     DoFKey const &_strain_key,
                     bool auto_range,
                     bool trim_corners);


    std::string name() const override {
      return enumerator_name;
    }

    static const std::string enumerator_name;

    static std::string interface_help();

    static int run(PrimClex const &primclex,
                   jsonParser const &kwargs,
                   Completer::EnumOption const &enum_opt);

    static int run(PrimClex const &_primclex,
                   Configuration const &_config,
                   Eigen::Ref<const Eigen::MatrixXd> const &_axes,
                   Eigen::Ref<const Eigen::VectorXd> const &min_val,
                   Eigen::Ref<const Eigen::VectorXd> const &max_val,
                   Eigen::Ref<const Eigen::VectorXd> const &inc_val,
                   bool sym_axes,
                   bool auto_range,
                   bool trim_corners,
                   std::vector<std::string> const &_filter_expr,
                   bool dry_run);
  private:

    /// Implements increment over all strain states
    void increment() override;


    // -- Unique -------------------
    DoFKey m_strain_key;

    bool m_trim_corners;

    Configuration m_current;

    // counts over strain grid
    EigenCounter<Eigen::VectorXd> m_counter;

    // counts over transformation matrices
    Index m_equiv_ind;

    std::vector<SymRepTools::SubWedge> m_wedges;

    //set of non-equivalent transformation matrices matrices that, along with m_counter define irreducible space
    //std::vector<Eigen::MatrixXd> m_trans_mats;
    PermuteIterator m_perm_begin, m_perm_end;
    Eigen::MatrixXd m_shape_factor;

    const PermuteIterator &_perm_begin() {
      return m_perm_begin;
    }
    const PermuteIterator &_perm_end() {
      return m_perm_end;
    }

  };

}

#endif
