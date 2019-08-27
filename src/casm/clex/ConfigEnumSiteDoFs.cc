#include "casm/clex/ConfigEnumSiteDoFs.hh"
#include "casm/crystallography/Structure.hh"
#include "casm/clex/PrimClex.hh"
#include "casm/symmetry/SupercellSymInfo_impl.hh"
#include "casm/enumerator/Enumerator_impl.hh"
#include "casm/database/ConfigDatabase.hh"
#include "casm/clex/ScelEnum.hh"
//#include "casm/clex/Supercell.hh"
#include "casm/clex/ConfigIsEquivalent.hh"
//#include "casm/misc/CASM_math.hh"
//#include "casm/misc/CASM_Eigen_math.hh"
//#include "casm/misc/algorithm.hh"

extern "C" {
  CASM::EnumInterfaceBase *make_ConfigEnumSiteDoFs_interface() {
    return new CASM::EnumInterface<CASM::ConfigEnumSiteDoFs>();
  }
}

namespace CASM {

  struct MakeConfigInvariantSubgroup {

    MakeConfigInvariantSubgroup() {}

    template<typename PermuteOutputIterator>
    PermuteOutputIterator operator()(const Configuration &config, PermuteIterator begin, PermuteIterator end, PermuteOutputIterator result) {
      ConfigIsEquivalent f(config, config.crystallography_tol());
      return std::copy_if(begin, end, result, f);
    }

  };

  const std::string ConfigEnumSiteDoFs::enumerator_name = "ConfigEnumSiteDoFs";

  std::string ConfigEnumSiteDoFs::interface_help() {
    return
      "ConfigEnumSiteDoFs: \n\n"

      "  confignames: Array of strings (optional) \n"
      "    Names of configurations to be used as reference states. Normal coordinates are enum-\n"
      "    erated after zeroing the DoF values at selected sites of the specified configurations\n"
      "    and calculating the resulting symmetry of the selected sites.\n"
      "    Ex: \"confignames\" : [\"SCEL1_1_1_1_0_0_0/1\",\"SCEL2_2_1_1_0_0_0/3\"]\n\n"

      "  scelnames: Array of strings (optional) \n"
      "    Names of supercells used as reference states. Normal coordinates are enumerated starting\n"
      "    from the fully zeroed configuration of the specified supercells.\n"
      "    Ex: \"scelnames\" : [\"SCEL1_1_1_1_0_0_0\",\"SCEL2_2_1_1_0_0_0\"]\n\n"

      "  dof: string (required) \n"
      "    Name of site degree of freecom  for which normal coordinates are to be generated.\n"
      "    Must be one of the degrees of freedom under consideration in the current project,\n"
      "    as determined by prim.json\n\n"

      "  sublats: array of integers (optional, default none) \n"
      "    Restricts normal coordinate determination to specified sublattices. Each sublat-\n"
      "    tice index specifies the correspondings basis site in prim.json, indexed from 0.\n"
      "    Ex: \"sublats\" : [0,2]\n\n"

      "  sites: array of 4-entry integer arrays (optional, default none) \n"
      "    Restricts normal coordinate determination to specified sites. Sites are specified\n"
      "    in [b,i,j,k] convention, where 'b' is sublattice index and [i,j,k] spedifies line-\n"
      "    ar combinations of primitive-cell lattice vectors.\n"
      "    Ex: \"sites\" : [[0,0,0,0],\n"
      "                   [2,0,0,0]]\n\n"

      "  axes: matrix of doubles (optional, default=identity matrix) \n"
      "    Coordinate axes of dof grid. Rows of matrix specify unrolled elements of\n"
      "    MxN DoF matrix, where M is dimension of on-site DoF and N is number of sites\n"
      "    in specified configuration or supercell. DoF values are unrolled in column-\n"
      "    major order, such that values from a particular site are listed contiguously.\n"
      "    'axes' matrix be rank deficient.\n"
      "    Ex: \"axes\" : [[1, 1, 1, 1, 1, 1],\n"
      "                  [1,-1, 0,-1, 1, 0],\n"
      "                  [1,-1, 0, 1,-1, 0]]\n\n"

      "  sym_axes: bool (optional, default=false)\n"
      "    If true, overrides \"axes\" field and instead constructs symmetry-adapted grid axes\n"
      "    as the symmetry-adapted DoF order parameters of 'config'. Run with option \n"
      "    \"dry_run\": true to obtain analysis report including the symmetry-adapted axes\n"
      "    without adding enumerated configurations to project.\n\n"

      "  min: number, or array of numbers (optional, default = [0,...,0]) \n"
      "    Minimum, starting value of grid counter\n"
      "    Dimension must be equal to number of rows of \"axes\"\n"
      "    Ex: \"min\" : [-0.05, -0.1, -0.1]\n\n"

      "  max: number, or array of numbers (required) \n"
      "    Maximum, final value of grid counter\n"
      "    Dimension must be equal to number of rows of \"axes\"\n"
      "    Ex: \"max\" : [0.05, 0.1, 0.1]\n\n"

      "  increment: number, or array of numbers (required) \n"
      "    Amount by which to increment counter elements\n"
      "    Dimension must be equal to number of rows of \"axes\"\n"
      "    Ex: \"increment\" : [0.01, 0.01, 0.01]\n\n"

      "  min_nonzero: integer (optional, default = 0) \n"
      "    Minimum number of coordinate amplitudes that are allowed\n"
      "    to be nonzero. Must be less than or equal to number of rows of \"axes\".\n\n"

      "  max_nonzero: integer (optional, default = axes.rows()) \n"
      "    Maximum number of coordinate amplitudes that are allowed\n"
      "    to be nonzero. Must be less than or equal to number of rows of \"axes\".\n\n"

      "  filter: string (optional, default=None)\n"
      "    A query command to use to filter which Configurations are kept.          \n\n"

      "  dry_run: bool (optional, default=false)\n"
      "    Perform dry run.\n\n"

      "  Examples:\n"
      "    To enumerate all DoF perturbations of a particular configuration:\n"
      "      casm enum --method ConfigEnumSiteDoFs -i \n"
      "      '{ \n"
      "        \"config\": \"SCEL4_1_4_1_0_0_0/3\",\n"
      "        \"analysis\": true,\n"
      "        } \n"
      "      }' \n\n";
  }

  int ConfigEnumSiteDoFs::run(
    PrimClex const &primclex,
    jsonParser const &_kwargs,
    Completer::EnumOption const &enum_opt,
    EnumeratorMap const *interface_map) {

    std::vector<ConfigEnumInput> in_configs = make_enumerator_input_configs(primclex, _kwargs, enum_opt, interface_map);
    std::vector<std::string> filter_expr = make_enumerator_filter_expr(_kwargs, enum_opt);

    DoFKey dof;
    if(!in_configs.size()) {
      return 1;
    }

    bool sym_axes(false);

    Index nsites = in_configs[0].sites().size();
    std::vector<Index> dof_dims;
    Index tot_dim(0);
    for(auto const &_in : in_configs) {
      if(_in.sites().size() != nsites) {
        throw std::runtime_error("Starting configurations or supercells passed to ConfigEnumSiteDoFs must all have the same number of selected sites!\n");
      }
    }

    Index max_nonzero(-1), min_nonzero(0);
    Eigen::MatrixXd axes;
    Eigen::VectorXd min, max, inc;
    try {
      if(!_kwargs.contains("dof")) {
        throw std::runtime_error("Field \"dof\" is required.\n");
      }
      from_json(dof, _kwargs["dof"]);
      DoF::traits(dof);

      auto const &dof_info = in_configs[0].configdof().local_dof(dof).info();
      for(Index l : in_configs[0].sites()) {
        dof_dims.push_back(dof_info[in_configs[0].config().sublat(l)].dim());
        tot_dim += dof_dims.back();
      }

      if(_kwargs.contains("sym_axes")) {
        _kwargs["sym_axes"].get(sym_axes);
      }

      if(!_kwargs.contains("axes")) {
        axes = Eigen::MatrixXd::Identity(tot_dim, tot_dim);
      }
      else {
        axes = _kwargs["axes"].get<Eigen::MatrixXd>().transpose();
        if(axes.rows() != tot_dim) {
          throw std::runtime_error("Number of columns of \"axes\" must be equal to dimensionality of selected variable space ("
                                   + std::to_string(tot_dim) + "). Size as parsed: " + std::to_string(axes.rows()));
        }
        if(axes.cols() > tot_dim) {
          throw std::runtime_error("Number of coordinate axes (i.e., number of rows of field \"axes\") must be less than or equal to dimensionality of selected variable space ("
                                   + std::to_string(tot_dim) + "). Number of axes parsed: " + std::to_string(axes.cols()));
        }
      }

      //min
      if(_kwargs.contains("min")) {

        if(_kwargs["min"].is_number()) {
          min = Eigen::VectorXd::Constant(axes.cols(), _kwargs["min"].get<double>());
        }
        else {
          _kwargs["min"].get(min);
          if(min.size() != axes.cols()) {
            throw std::runtime_error("Array field \"min\" must have dimension equal to number of coordinate axes!");
          }
        }
      }
      else {
        min = Eigen::VectorXd::Constant(axes.cols(), 0);
      }

      //max
      if(!_kwargs.contains("max")) {
        throw std::runtime_error("Field \"max\" is required.\n");
      }
      if(_kwargs["max"].is_number()) {
        max = Eigen::VectorXd::Constant(axes.cols(), _kwargs["max"].get<double>());
      }
      else {
        _kwargs["max"].get(max);
        if(max.size() != axes.cols()) {
          throw std::runtime_error("Array field \"max\" must have dimension equal to number of coordinate axes!");
        }
      }

      //inc
      if(!_kwargs.contains("increment")) {
        throw std::runtime_error("Field \"increment\" is required.\n");
      }
      if(_kwargs["increment"].is_number()) {
        inc = Eigen::VectorXd::Constant(axes.cols(), _kwargs["increment"].get<double>());
      }
      else {
        _kwargs["increment"].get(inc);
        if(inc.size() != axes.cols()) {
          throw std::runtime_error("Array field \"increment\" must have dimension equal to number of coordinate axes!");
        }

      }

      _kwargs.get_if(min_nonzero, "min_nonzero");

      _kwargs.get_else(max_nonzero, "max_nonzero", axes.cols());

    }
    catch(std::exception &e) {
      throw std::runtime_error(std::string("Error parsing JSON arguments for ConfigEnumSiteDoFs: ") + e.what());
    }

    for(ConfigEnumInput const &config : in_configs) {
      Index result = run(primclex,
                         config,
                         dof,
                         axes,
                         min,
                         max,
                         inc,
                         sym_axes,
                         min_nonzero,
                         max_nonzero,
                         filter_expr,
                         CASM::dry_run(_kwargs, enum_opt));
      if(result)
        return result;
    }

    return 0;
  }

  int ConfigEnumSiteDoFs::run(PrimClex const &_primclex,
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
                              bool dry_run) {
    Configuration tconfig = _in_config.config();

    if(_in_config.sites().size() == 0) {
      tconfig.configdof().local_dof(_dof).values().setZero();
    }
    else {
      for(Index s : _in_config.sites()) {
        tconfig.configdof().local_dof(_dof).site_value(s).setZero();
      }
    }

    ConfigEnumInput config(tconfig, _in_config.sites());
    Eigen::MatrixXd axes = _axes;
    //PRINT INFO TO LOG:
    Log &log = _primclex.log();
    Eigen::IOFormat tformat(4, 0, 8, " ", "\n", "    ", "", "", "");
    if(sym_axes) {
      log << "Option \"sym_axes\" selected. Preparing to construct symmetry-adapted axes. This may take several minutes...\n\n";
      std::pair<Eigen::MatrixXd, std::vector<Index>> normcoords = collective_dof_normal_coords_and_irrep_dims(config.sites().begin(),
                                                                  config.sites().end(),
                                                                  config.supercell().sym_info(),
                                                                  _dof,
                                                                  config.group(),
                                                                  _axes);
      axes = normcoords.first.transpose();
      log << "ConfigEnumSiteDoFs summary for DoF '" << _dof << "':\n";

      //std::cout << "Axes:\n" << axes.transpose().format(tformat) << "\n";
      log << "Enumeration will be performed using symmetry-adapted normal coordinates as axes.\n"
          << "Normal coordinates partition DoF space into " << normcoords.second.size() << " subspaces.\n"
          << "Normal coordinates are:\n";
      Index l = 0;
      for(Index d = 0; d < normcoords.second.size(); ++d) {
        log << "Axes for irreducible representation " << (d + 1) << "\n  --------------\n";
        Index dim = normcoords.second[d];
        for(Index i = 0; i < dim; ++i, ++l) {
          log <<  axes.col(l).transpose().format(tformat) << "\n";
        }
      }
      if(axes.cols() != _axes.cols()) {
        throw std::runtime_error("In ConfigEnumSiteDoFs, symmetry-adapted axes do not have same dimension as provided axes. "
                                 "Please ensure that provided axes completely span one or more of subspaces listed above.");
      }
    }
    else {
      log << "Enumeration will be performed using user-specified axes. Axes are:\n";
      log << axes.transpose().format(tformat) << "\n";
    }
    log << "\nEnumeration will be performed from starting vector: \n" << min_val.transpose().format(tformat) << "\n";
    log << "\nEnumeration will be performed to final vector: \n" << max_val.transpose().format(tformat) << "\n";
    log << "\nEnumeration will be performed using grid increment vector: \n" << inc_val.transpose().format(tformat) << "\n";
    log << "---------\nEnd of summary.\n\n";


    auto constructor = [&](const ConfigEnumInput & _config) {
      return notstd::make_unique<ConfigEnumSiteDoFs>(config,
                                                     _dof,
                                                     axes,
                                                     min_val,
                                                     max_val,
                                                     inc_val,
                                                     _min_nonzero,
                                                     _max_nonzero);
    };

    int returncode = insert_configs(enumerator_name,
                                    _primclex,
                                    config.supercell(),
                                    constructor,
                                    _filter_expr,
                                    false,
                                    dry_run);

    return returncode;

  }

  ConfigEnumSiteDoFs::ConfigEnumSiteDoFs(ConfigEnumInput const &_init,
                                         DoFKey const &_dof,
                                         Eigen::Ref<const Eigen::MatrixXd> const &_axes,
                                         Eigen::Ref<const Eigen::VectorXd> const &min_val,
                                         Eigen::Ref<const Eigen::VectorXd> const &max_val,
                                         Eigen::Ref<const Eigen::VectorXd> const &inc_val,
                                         Index _min_nonzero,
                                         Index _max_nonzero) :

    //m_current(_init.config()),
    m_dof_key(_dof),
    m_min_nonzero(_min_nonzero),
    m_max_nonzero(_max_nonzero),
    m_axes(_axes),
    m_min(min_val),
    m_max(max_val),
    m_inc(inc_val),
    m_unit_length(DoFType::traits(_dof).unit_length()),
    m_sites(_init.sites().begin(), _init.sites().end()),
    m_subset_mode(false),
    //m_combo_mode(false),
    m_combo_index(0) {



    if(_init.config().size() != _init.sites().size())
      m_subset_mode = true;



    if(m_axes.cols() == 0) {
      this->_invalidate();
    }
    else {
      m_current = notstd::make_cloneable<Configuration>(_init.config());
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

