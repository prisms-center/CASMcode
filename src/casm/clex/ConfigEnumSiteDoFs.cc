#include "casm/basis_set/DoFTraits.hh"
#include "casm/clex/ConfigEnumSiteDoFs.hh"
#include "casm/enumerator/ConfigEnumInput.hh"
#include "casm/symmetry/SupercellSymInfo_impl.hh"


namespace CASM {

  /// See `ConfigEnumSiteDoFsParams` for method and parameter details
  ConfigEnumSiteDoFs::ConfigEnumSiteDoFs(
    ConfigEnumInput const &_in_config,
    ConfigEnumSiteDoFsParams const &params):
    ConfigEnumSiteDoFs(
      _in_config,
      params.dof,
      params.axes,
      params.min_val,
      params.max_val,
      params.inc_val,
      params.min_nonzero,
      params.max_nonzero) {}


  /// See `ConfigEnumSiteDoFsParams` for method and parameter details
  ConfigEnumSiteDoFs::ConfigEnumSiteDoFs(ConfigEnumInput const &_init,
                                         DoFKey const &_dof,
                                         Eigen::Ref<const Eigen::MatrixXd> const &_axes,
                                         Eigen::Ref<const Eigen::VectorXd> const &min_val,
                                         Eigen::Ref<const Eigen::VectorXd> const &max_val,
                                         Eigen::Ref<const Eigen::VectorXd> const &inc_val,
                                         Index _min_nonzero,
                                         Index _max_nonzero) :

    m_dof_key(_dof),
    m_min_nonzero(_min_nonzero),
    m_max_nonzero(_max_nonzero),
    m_axes(_axes),
    m_min(min_val),
    m_max(max_val),
    m_inc(inc_val),
    m_unit_length(DoF::BasicTraits(_dof).unit_length()),
    m_sites(_init.sites().begin(), _init.sites().end()),
    m_subset_mode(false),
    //m_combo_mode(false),
    m_combo_index(0) {

    Configuration const &configuration = _init.configuration();

    if(configuration.size() != m_sites.size())
      m_subset_mode = true;

    if(m_axes.cols() == 0) {
      this->_invalidate();
    }
    else {
      m_current = notstd::make_cloneable<Configuration>(configuration);
      reset_properties(*m_current);

      m_dof_vals = &(m_current->configdof().local_dof(m_dof_key));
      auto const &dof_info = m_dof_vals->info();
      for(Index l : m_sites)
        m_dof_dims.push_back(dof_info[m_current->sublat(l)].dim());

      if(m_max_nonzero > m_axes.cols() / 3) {
        m_combo.resize(m_max_nonzero - 1);
        m_combo_index = m_max_nonzero;
      }
      if(_increment_combo())
        this->_initialize(&(*m_current));
      _set_dof();
      m_current->set_source(this->source(step()));
    }
  }

  std::string ConfigEnumSiteDoFs::name() const {
    return enumerator_name;
  }

  const std::string ConfigEnumSiteDoFs::enumerator_name = "ConfigEnumSiteDoFs";

  bool ConfigEnumSiteDoFs::_increment_combo() {
    Index k = max<Index>(m_combo.size(), m_min_nonzero);
    bool invalid = true;
    //std::cout << "COMBO INCREMENT: " << m_combo << "  to  ";
    while(invalid) {
      if(m_combo_index >= nchoosek(m_axes.cols(), k)) {
        ++k;
        m_combo_index = 0;
      }
      if(k > m_max_nonzero)
        break;
      invalid = false;
      m_combo = index_to_kcombination(m_combo_index, k);
      Eigen::VectorXd
      vmin(m_combo.size()),
           vmax(m_combo.size()),
           vinc(m_combo.size());
      for(Index i = 0; i < m_combo.size(); ++i) {
        Index j = m_combo.size() - 1 - i;
        vmin[i] = m_min[m_combo[j]];
        vmax[i] = m_max[m_combo[j]];
        vinc[i] = m_inc[m_combo[j]];
        if(almost_zero(vmin[i]) && almost_zero(vmax[i]))
          invalid = true;
      }
      m_combo_index++;
      m_counter = EigenCounter<Eigen::VectorXd>(vmin, vmax, vinc);
    }

    //std::cout << m_combo << "\n";

    //std::cout << "counter.valid() = " << m_counter.valid() << "; values: " << m_counter().transpose() << "\n";
    return !invalid;
  }


  void ConfigEnumSiteDoFs::_set_dof() {
    Eigen::MatrixXd vals = m_current->configdof().local_dof(m_dof_key).values();
    Eigen::VectorXd pert_vals(Eigen::VectorXd::Zero(m_axes.rows()));

    for(Index i = 0; i < m_combo.size(); ++i)
      pert_vals += m_counter[i] * m_axes.col(m_combo[i]);
    Index l = 0;

    for(Index i = 0; i < m_sites.size(); ++i) {
      for(Index j = 0; j < m_dof_dims[i]; ++j, ++l) {
        vals(j, m_sites[i]) = pert_vals(l);
      }
      if(m_unit_length) {
        double tnorm = vals.col(m_sites[i]).norm();
        if(!almost_zero(tnorm))
          vals.col(m_sites[i]) /= tnorm;
      }
    }
    m_current->configdof().set_local_dof(m_dof_key, vals);
  }

  /// Implements _increment over all occupations
  void ConfigEnumSiteDoFs::increment() {
    //std::cout << "m_combo before increment: " << m_combo << "\n";
    bool is_valid_config = false;
    while(!is_valid_config && (++m_counter || _increment_combo())) {
      //std::cout << "Counting point: " << m_counter() << "\n";
      if(_check_sparsity()) {
        //std::cout << "Sparsity is good!\n";
        _set_dof();
        is_valid_config = _check_current();
      }
      //std::cout << "is_valid? " << (is_valid_config? "yes\n" : "no\n");
    }
    //std::cout << "m_combo after increment: " << m_combo << "\n";
    //std::cout << "OUTSIDE!! is_valid? " << (is_valid_config? "yes\n" : "no\n");
    if(is_valid_config) {
      this->_increment_step();
      m_current->set_source(this->source(step()));
    }
    else {
      this->_invalidate();
    }
  }

  bool ConfigEnumSiteDoFs::_check_sparsity() const {
    if(m_min_nonzero == 0 && m_max_nonzero == m_axes.cols())
      return true;

    Index nonzero = 0;
    for(Index i = 0; i < m_counter().size(); ++i) {
      nonzero += almost_zero(m_counter[i]) ? 0 : 1;
    }
    return (m_min_nonzero <= nonzero && nonzero <= m_max_nonzero);
  }

  /// Returns true if current() is primitive and canonical
  bool ConfigEnumSiteDoFs::_check_current() const {
    return true;
    //return current().is_primitive() && _check_sparsity() && (m_subset_mode || current().is_canonical());
  }


}
