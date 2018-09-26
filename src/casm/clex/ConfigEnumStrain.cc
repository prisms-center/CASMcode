#include "casm/clex/ConfigEnumStrain.hh"
#include <algorithm>
#include "casm/crystallography/Structure.hh"
#include "casm/symmetry/SymRepTools.hh"
#include "casm/clex/PrimClex.hh"
#include "casm/database/ConfigDatabase.hh"
#include "casm/clex/Supercell.hh"
#include "casm/clex/ConfigIsEquivalent.hh"
#include "casm/clex/PrimClex.hh"
#include "casm/misc/CASM_math.hh"
#include "casm/misc/CASM_Eigen_math.hh"
#include "casm/misc/algorithm.hh"

extern "C" {
  CASM::EnumInterfaceBase *make_ConfigEnumStrain_interface() {
    return new CASM::EnumInterface<CASM::ConfigEnumStrain>();
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

  const std::string ConfigEnumStrain::enumerator_name = "ConfigEnumStrain";

  std::string ConfigEnumStrain::interface_help() {
    return
      "ConfigEnumAllStrain: \n\n"

      "  config: string (required) \n"
      "    Name of configuration used as strain reference\n"
      "    Must be single configuration. Strains are enumerated as perturbations to\n"
      "    specified configuration. Ex: \"config\" : \"SCEL2_2_1_1_0_0_0\"\n\n"

      "  min: array of doubles (optional, default = [0,...,0]) \n"
      "    Minimum, starting value of grid counter\n"
      "    Dimension must be either 6, or equal to number of rows of \"axes\"\n"
      "    Ex: \"min\" : [-0.1, -0.1, -0.1, -0.01, -0.01, -0.01]\n\n"

      "  max: array of doubles (required) \n"
      "    Maximum, final value of grid counter\n"
      "    Dimension must be either 6, or equal to number of rows of \"axes\"\n"
      "    Ex: \"min\" : [0.1, 0.1, 0.1, 0.01, 0.01, 0.01]\n\n"

      "  increment: array of doubles (required) \n"
      "    Amount by which to increment counter elements\n"
      "    Dimension must be either 6, or equal to number of rows of \"axes\"\n"
      "    Ex: \"increment\" : [0.01, 0.01, 0.01, 0.01, 0.01, 0.01]\n\n"

      "  axes: matrix of doubles (optional, default=6x6 identity matrix) \n"
      "    Coordinate axes of strain grid. Rows of matrix specify linear combinations of\n"
      "    strain basis specified in prim.json. May be rank deficient.\n"
      "    Ex: \"axes\" : [[1,0,0,0,0,0],\n"
      "                    [0,1,0,0,0,0],\n"
      "                    [0,0,1,0,0,0]]\n\n"

      "  sym_axes: bool (optional, default=false)\n"
      "    If true, overrides \"axes\" field and instead constructs symmetry-adapted grid axes\n"
      "    as the symmetry-adapted strain order parameters of 'config'. Running with option \n"
      "    \"analysis\": true will display information including the symmetry-adapted axes.\n\n"

      "  trim_corners: bool (optional, default=true) Exclude extreme grid points, if true\n"
      "    If true, any grid points outside the largest ellipsoid inscribed within the extrema\n"
      "    of the grid will be discarded\n\n"

      "  config: Name of configuration to which strain enumerations will be applied\n"
      "    Must be single configuration. Ex: \"config\" : \"SCEL2_2_1_1_0_0_0\"\n\n"

      "  filter: string (optional, default=None)\n"
      "    A query command to use to filter which Configurations are kept.          \n\n"

      "  dry_run: bool (optional, default=false)\n"
      "    Perform dry run.\n\n"

      "  analysis: bool (optional, default=false)\n"
      "    Print symmetry analysis info for chosen configuration.\n\n"

      "  Examples:\n"
      "    To enumerate all strain perturbations of a particular configuration:\n"
      "      casm enum --method ConfigEnumStrain -i \n"
      "      '{ \n"
      "        \"config\": \"SCEL4_1_4_1_0_0_0/3\",\n"
      "        \"analysis\": true,\n"
      "        } \n"
      "      }' \n\n";
  }

  int ConfigEnumStrain::run(
    const PrimClex &primclex,
    const jsonParser &_kwargs,
    const Completer::EnumOption &enum_opt) {
    bool trim_corners = _kwargs.get_if_else<bool>("trim_corners", true);
    bool sym_axes = _kwargs.get_if_else<bool>("sym_axes", false);
    bool analysis = _kwargs.get_if_else<bool>("analysis", false);
    Eigen::MatrixXd axes;
    Eigen::VectorXd min_val, max_val, inc_val;

    bool auto_range = false;

    if(!analysis) {
      if(!(_kwargs.contains("min") && _kwargs.contains("max") && _kwargs.contains("increment")))
        throw std::runtime_error("JSON options for enumeration method 'ConfigEnumStrain' must specify BOTH \"max\", and \"increment\"");
      from_json(max_val, _kwargs["max"]);
      from_json(inc_val, _kwargs["increment"]);
      if(_kwargs.contains("min")) {
        from_json(min_val, _kwargs["min"]);
      }
      else {
        auto_range = sym_axes;
        min_val = 0 * max_val;
      }
    }

    if(!_kwargs.contains("config") || !_kwargs["config"].is_string())
      throw std::runtime_error("JSON options for enumeration method 'ConfigEnumStrain' must include exactly one starting configuration, specified by field \"config\"");

    std::vector<std::string> filter_expr = make_enumerator_filter_expr(_kwargs, enum_opt);

    std::string configname = _kwargs["config"].get<std::string>();
    auto it = primclex.const_db<Configuration>().find(configname);
    if(it == primclex.const_db<Configuration>().end())
      throw std::runtime_error("In ConfigEnumStrain::run(), specified configuration " + configname + " does not exist in database.\n");

    return run(primclex,
               *it,
               axes,
               min_val,
               max_val,
               inc_val,
               sym_axes,
               auto_range,
               trim_corners,
               analysis,
               filter_expr,
               CASM::dry_run(_kwargs, enum_opt));
  }

  int ConfigEnumStrain::run(PrimClex const &_primclex,
                            Configuration const &_config,
                            Eigen::Ref<const Eigen::MatrixXd> const &_axes,
                            Eigen::Ref<const Eigen::VectorXd> const &min_val,
                            Eigen::Ref<const Eigen::VectorXd> const &max_val,
                            Eigen::Ref<const Eigen::VectorXd> const &inc_val,
                            bool sym_axes,
                            bool auto_range,
                            bool trim_corners,
                            bool analysis,
                            std::vector<std::string> const &_filter_expr,
                            bool dry_run) {

    std::vector<SymRepTools::SubWedge> wedges;
    std::vector<int> dims;
    SymGroup pg = make_sym_group(_config.point_group());
    if(!analysis && ! sym_axes)
      wedges.push_back(SymRepTools::SubWedge({SymRepTools::IrrepWedge(_axes, std::vector<Index>(_axes.cols(), 1))}));
    else {
      DoFSet const *dof_ptr = get_strain_dof(_primclex.prim());
      if(!dof_ptr)
        throw std::runtime_error("Cannot enumerate strains for project in which strain has not been specified as a degree of freedom.");
      wedges = SymRepTools::symrep_subwedges(pg, dof_ptr->symrep_ID());
    }

    if(analysis) {
      Log &log = _primclex.log();
      //PRINT INFO TO LOG:

      return 0;
    }
    auto constructor = [&](const Supercell & scel) {
      return notstd::make_unique<ConfigEnumStrain>(_config,
                                                   wedges,
                                                   min_val,
                                                   max_val,
                                                   inc_val,
                                                   auto_range,
                                                   trim_corners);
    };

    int returncode = insert_unique_canon_configs(enumerator_name,
                                                 _primclex,
                                                 _config.supercell(),
                                                 constructor,
                                                 _filter_expr,
                                                 dry_run);



  }

  ConfigEnumStrain::ConfigEnumStrain(const Configuration &_init,
                                     const std::vector<SymRepTools::SubWedge> &_wedges,
                                     Eigen::VectorXd min_val,
                                     Eigen::VectorXd max_val,
                                     Eigen::VectorXd inc_val,
                                     bool auto_range,
                                     bool trim_corners) :
    m_trim_corners(trim_corners),
    m_current(_init),
    m_equiv_ind(0),
    //m_perm_begin(_scel.permute_begin()),
    //m_perm_end(_scel.permute_end()),
    m_shape_factor(Eigen::MatrixXd::Identity(min_val.size(), min_val.size())) {

    //Condition range arrays and build shape_factor matrix
    Index nc = 0;
    if(_wedges.size() == 0)
      return;
    Index nsub = _wedges[0].irrep_wedges().size();
    for(Index s = 0; s < nsub; s++) {
      double abs_mag = 0.;
      for(Index i = 0; i < _wedges[0].irrep_wedges()[s].axes.cols(); i++)
        abs_mag = max(abs_mag, max(abs(max_val(nc)), abs(min_val(nc))));

      for(Index i = 0; i < _wedges[0].irrep_wedges()[s].axes.cols(); i++, nc++) {
        if(_wedges[0].irrep_wedges()[s].mult[i] == 1 && auto_range) {
          min_val(nc) = -max_val(nc);
        }

        if(abs_mag > TOL)
          m_shape_factor(nc, nc) /= abs_mag * abs_mag;

        if(inc_val(nc) == 0. || (max_val(nc) - min_val(nc)) / inc_val(nc) > 1e4) {
          max_val(nc) = min_val(nc) + TOL;
          inc_val(nc) = 10 * TOL;
        }
        //else {
        //max_val(nc) = absmags[s] + TOL;
        //}
      }
    }


    //m_shape_factor = m_strain_calc.sop_transf_mat().transpose() * m_shape_factor * m_strain_calc.sop_transf_mat();
    //Find first valid config
    m_counter = EigenCounter<Eigen::VectorXd>(min_val, max_val, inc_val);

    m_counter.reset();
    while(m_counter.valid() && (trim_corners && double(m_counter().transpose()*m_trans_mats[m_equiv_ind].transpose()*m_shape_factor * m_trans_mats[m_equiv_ind]*m_counter()) > 1.0 + TOL)) {
      ++m_counter;
    }

    reset_properties(m_current);
    this->_initialize(&m_current);

    if(!m_counter.valid()) {
      this->_invalidate();
    }
    m_current.set_source(this->source(step()));
  }

  // Implements _increment
  void ConfigEnumStrain::increment() {
    //bool is_valid_config(false);
    //std::cout << "Incrementing...\n";

    while(++m_counter && (m_trim_corners && double(m_counter().transpose()*m_trans_mats[m_equiv_ind].transpose()*m_shape_factor * m_trans_mats[m_equiv_ind]*m_counter()) > 1.0 + TOL)) {
      //just burning throught the count
    }

    // move to next part of wedge if necessary
    if(!m_counter.valid() && m_equiv_ind + 1 < m_trans_mats.size()) {
      m_counter.reset();
      ++m_equiv_ind;
      //std::cout << "INCREMENTING m_equiv_ind to " << m_equiv_ind << "\n";
    }

    while(m_counter &&
          (m_trim_corners && double(m_counter().transpose()*m_trans_mats[m_equiv_ind].transpose()*m_shape_factor * m_trans_mats[m_equiv_ind]*m_counter()) > 1.0 + TOL)) {
      //just burning throught the count
      ++m_counter;
    }

    if(m_counter.valid()) {
      throw std::runtime_error("UPDATE STRAIN INSERTION");
      //m_current.set_deformation(m_strain_calc.unrolled_strain_metric_to_F(m_trans_mats[m_equiv_ind] * m_counter()));
      std::cout << "Counter is " << m_counter().transpose() << "\n\n";
      //std::cout << "strain vector is \n" << m_trans_mats[m_equiv_ind]*m_counter() << "\n\n";
      //std::cout << "DEFORMATION IS\n" << m_current.deformation() << "\n\n";
      //is_valid_config = current().is_canonical(_perm_begin(), _perm_end());
      //std::cout << "counter() is: " << m_counter() << ";  is_valid_config: " << is_valid_config
      //<< ";  is_valid_counter: " << m_counter.valid() << "\n";
      _increment_step();
    }
    else {
      //std::cout << "REACHED END OF THE LINE!\n";
      _invalidate();
    }
    m_current.set_source(this->source(step()));
    //std::cout << "--FINISHED SEARCH " << _step()<< "--\n";
    return;
  }

}

