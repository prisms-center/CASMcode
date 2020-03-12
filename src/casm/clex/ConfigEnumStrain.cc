#include "casm/clex/ConfigEnumStrain.hh"
#include <algorithm>
#include "casm/crystallography/Structure.hh"
#include "casm/symmetry/SymRepTools.hh"
#include "casm/clex/PrimClex.hh"
#include "casm/enumerator/Enumerator_impl.hh"
#include "casm/database/ConfigDatabase.hh"
#include "casm/clex/ScelEnum.hh"
#include "casm/clex/Supercell.hh"
#include "casm/clex/ConfigIsEquivalent.hh"
#include "casm/misc/CASM_math.hh"
#include "casm/misc/CASM_Eigen_math.hh"
#include "casm/misc/algorithm.hh"
#include "casm/crystallography/BasicStructure_impl.hh"

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
      "ConfigEnumStrain: \n\n"

      "  confignames: Array of strings (optional) \n"
      "    Names of configurations used as initial states. Strains are enumerated\n"
      "    as perturbations to specified configurations. \n"
      "    Ex: \"confignames\" : [\"SCEL2_2_1_1_0_0_0/3\", \"SCEL4_2_2_1_0_0_0/10\"]\n\n"

      "  scelnames: Array of strings (optional) \n"
      "    Names of supercells used as initial states. Strains are enumerated\n"
      "    as perturbations of the fully zeroed configuration of the specified"
      "    supercells.\n"
      "    Ex: \"scelnames\" : [\"SCEL2_2_1_1_0_0_0\",\"SCEL4_2_3_1_0_0_0\"]\n\n"

      "  filter: string (optional, default=None)\n"
      "    A query command to use to filter which Configurations are kept.          \n\n"

      "  dry_run: bool (optional, default=false)\n"
      "    Perform dry run and print symmetry analysis info for chosen configuration.\n\n"

      "  supercells: ScelEnum JSON settings (default='{\"existing_only\"=true}')\n"
      "    Indicate supercells to use as initial states of enumeration in terms of size\n"
      "    and unit cell via a JSON object conforming to the format of 'ScelEnum' JSON\n"
      "    settings. \"scelnames\" will override \"supercells\", but if neither is specified\n"
      "    all existing supercells are used by default. See 'ScelEnum' description for details.\n\n"

      "  min: array of doubles (optional, default = [0,...,0]) \n"
      "    Minimum, starting value of grid counter\n"
      "    Dimension must be either 6, or equal to number of rows of \"axes\"\n"
      "    Ex: \"min\" : [-0.1, -0.1, -0.1, -0.01, -0.01, -0.01]\n\n"

      "  max: array of doubles (required) \n"
      "    Maximum, final value of grid counter\n"
      "    Dimension must be either 6, or equal to number of rows of \"axes\"\n"
      "    Ex: \"max\" : [0.1, 0.1, 0.1, 0.01, 0.01, 0.01]\n\n"

      "  increment: array of doubles (required) \n"
      "    Amount by which to increment counter elements\n"
      "    Dimension must be either 6, or equal to number of rows of \"axes\"\n"
      "    Ex: \"increment\" : [0.01, 0.01, 0.01, 0.01, 0.01, 0.01]\n\n"

      "  axes: matrix of doubles (optional, default=identity matrix) \n"
      "    Coordinate axes of strain grid. Rows of matrix specify linear combinations of\n"
      "    strain basis specified in prim.json. May be rank deficient.\n"
      "    Ex: \"axes\" : [[1,0,0,0,0,0],\n"
      "                    [0,1,0,0,0,0],\n"
      "                    [0,0,1,0,0,0]]\n\n"

      "  sym_axes: bool (optional, default=false)\n"
      "    If true, overrides \"axes\" field and instead constructs symmetry-adapted grid axes\n"
      "    as the symmetry-adapted strain order parameters of 'config'. Run with option \n"
      "    \"dry_run\": true to obtain analysis report including the symmetry-adapted axes\n"
      "    without adding enumerated configurations to project.\n\n"

      "  trim_corners: bool (optional, default=true) Exclude extreme grid points, if true\n"
      "    If true, any grid points outside the largest ellipsoid inscribed within the extrema\n"
      "    of the grid will be discarded\n\n"

      "  Examples:\n"
      "    To enumerate all strain perturbations of a particular configuration:\n"
      "      casm enum --method ConfigEnumStrain -i \n"
      "      '{ \n"
      "        \"config\" : \"SCEL4_1_4_1_0_0_0/3\",\n"
      "        \"increment\" : [0.01, 0.01, 0.01, 0., 0., 0.],"
      "        \"max\" : [0.05, 0.05, 0,05, 0., 0., 0.]"
      "        } \n"
      "      }' \n\n";
  }

  int ConfigEnumStrain::run(
    const PrimClex &primclex,
    const jsonParser &_kwargs,
    const Completer::EnumOption &enum_opt,
    EnumeratorMap const *interface_map) {
    //bool analysis = _kwargs.get_if_else<bool>("analysis", false);

    DoFKey strain_dof_key;
    std::vector<DoFKey> tdof_types = global_dof_types(primclex.prim());
    Index istrain = find_index_if(tdof_types,
    [](DoFKey const & other) {
      return other.find("strain") != std::string::npos;
    });
    if(istrain == tdof_types.size())
      throw std::runtime_error("Cannot enumerate strains for project in which strain has not been specified as a degree of freedom.");
    strain_dof_key = tdof_types[istrain];
    Index dim = primclex.prim().global_dof(strain_dof_key).dim();

    Eigen::MatrixXd axes;
    Eigen::VectorXd min_val, max_val, inc_val;
    std::vector<std::string> confignames;

    bool auto_range = false;
    bool trim_corners = true;
    bool sym_axes = false;
    try {
      _kwargs.get_if(trim_corners, "trim_corners");
      _kwargs.get_if(sym_axes, "sym_axes");
      if(!(_kwargs.contains("axes"))) {
        axes.setIdentity(dim, dim);
      }
      else {
        axes = _kwargs["axes"].get<Eigen::MatrixXd>().transpose();
        if(axes.rows() != dim || axes.cols() > dim) {
          throw std::runtime_error("In JSON input for ConfigEnumStrain, field \"axes\" must be a matrix having exactly " + std::to_string(dim) +
                                   " columns and no more than " + std::to_string(dim) + " rows.");
        }
      }

      //min
      if(_kwargs.contains("min")) {

        if(_kwargs["min"].is_number()) {
          min_val = Eigen::VectorXd::Constant(axes.cols(), _kwargs["min"].get<double>());
        }
        else {
          _kwargs["min"].get(min_val);
          if(min_val.size() != axes.cols()) {
            throw std::runtime_error("Array field \"min\" must have dimension equal to number of coordinate axes!");
          }
        }
      }
      else {
        auto_range = sym_axes;
        min_val = Eigen::VectorXd::Constant(axes.cols(), 0);
      }

      //max
      if(!_kwargs.contains("max")) {
        throw std::runtime_error("Field \"max\" is required.\n");
      }
      if(_kwargs["max"].is_number()) {
        max_val = Eigen::VectorXd::Constant(axes.cols(), _kwargs["max"].get<double>());
      }
      else {
        _kwargs["max"].get(max_val);
        if(max_val.size() != axes.cols()) {
          throw std::runtime_error("Array field \"max\" must have dimension equal to number of coordinate axes!");
        }
      }

      //inc
      if(!_kwargs.contains("increment")) {
        throw std::runtime_error("Field \"increment\" is required.\n");
      }
      if(_kwargs["increment"].is_number()) {
        inc_val = Eigen::VectorXd::Constant(axes.cols(), _kwargs["increment"].get<double>());
      }
      else {
        _kwargs["increment"].get(inc_val);
        if(inc_val.size() != axes.cols()) {
          throw std::runtime_error("Array field \"increment\" must have dimension equal to number of coordinate axes!");
        }

      }

    }
    catch(std::exception &e) {
      throw std::runtime_error(std::string("Error parsing JSON arguments for ConfigStrain:") + e.what());
    }

    std::vector<ConfigEnumInput> in_configs = make_enumerator_input_configs(primclex, _kwargs, enum_opt, interface_map);

    std::vector<std::string> filter_expr = make_enumerator_filter_expr(_kwargs, enum_opt);

    for(ConfigEnumInput const &config : in_configs) {
      Index result = run(primclex,
                         config,
                         axes,
                         min_val,
                         max_val,
                         inc_val,
                         sym_axes,
                         auto_range,
                         trim_corners,
                         //analysis,
                         filter_expr,
                         CASM::dry_run(_kwargs, enum_opt));
      if(result)
        return result;
    }
    return 0;

  }

  int ConfigEnumStrain::run(PrimClex const &_primclex,
                            ConfigEnumInput const &_config,
                            Eigen::Ref<const Eigen::MatrixXd> const &_axes,
                            Eigen::Ref<const Eigen::VectorXd> const &min_val,
                            Eigen::Ref<const Eigen::VectorXd> const &max_val,
                            Eigen::Ref<const Eigen::VectorXd> const &inc_val,
                            bool sym_axes,
                            bool auto_range,
                            bool trim_corners,
                            //bool analysis,
                            std::vector<std::string> const &_filter_expr,
                            bool dry_run) {
    std::vector<SymRepTools::SubWedge> wedges;
    std::vector<int> dims;
    SymGroup pg = make_point_group(_config.group(), _config.supercell().sym_info().supercell_lattice());
    DoFKey strain_dof_key;
    std::vector<DoFKey> tdof_types = global_dof_types(_primclex.prim());
    Index istrain = find_index_if(tdof_types,
    [](DoFKey const & other) {
      return other.find("strain") != std::string::npos;
    });
    if(istrain == tdof_types.size())
      throw std::runtime_error("Cannot enumerate strains for project in which strain has not been specified as a degree of freedom.");
    strain_dof_key = tdof_types[istrain];
    if(!sym_axes)
      wedges.push_back(SymRepTools::SubWedge({SymRepTools::IrrepWedge(_axes, std::vector<Index>(_axes.cols(), 1))}));
    else
      wedges = SymRepTools::symrep_subwedges(pg, _primclex.prim().global_dof(strain_dof_key).symrep_ID()).first;

    //PRINT INFO TO LOG:
    Log &log = _primclex.log();
    log << "Strain enumeration summary:\n";
    if(sym_axes) {
      log << "  Strains will be enumerated within the " << wedges.size() << " symmetrically distinct sub-wedges\n";
    }
    else {
      log << "  Strains will be enumerated using the user-specified grid\n";
    }
    Index l = 1;
    Eigen::IOFormat tformat(4, 0, 8, " ", "\n", "    ", "", "", "");
    for(SymRepTools::SubWedge const &wedge : wedges) {
      if(sym_axes)
        log << "  Sub-wedge #" << l++ << ": \n";
      for(auto const &irr : wedge.irrep_wedges()) {
        log <<  irr.axes.transpose().format(tformat) << "\n   --------------\n";

      }
      log << "\n";
    }
    log << "End of summary.\n\n";

    auto constructor = [&](const ConfigEnumInput & in_config) {
      return notstd::make_unique<ConfigEnumStrain>(in_config,
                                                   wedges,
                                                   min_val,
                                                   max_val,
                                                   inc_val,
                                                   strain_dof_key,
                                                   auto_range,
                                                   trim_corners);
    };

    int returncode = insert_configs(enumerator_name,
                                    _primclex,
                                    _config,
                                    constructor,
                                    _filter_expr,
                                    false,
                                    dry_run);

    return returncode;

  }

  ConfigEnumStrain::ConfigEnumStrain(const ConfigEnumInput &_init,
                                     const std::vector<SymRepTools::SubWedge> &_wedges,
                                     Eigen::VectorXd min_val,
                                     Eigen::VectorXd max_val,
                                     Eigen::VectorXd inc_val,
                                     DoFKey const &_strain_key,
                                     bool auto_range,
                                     bool trim_corners) :
    m_strain_key(_strain_key),
    m_trim_corners(trim_corners),
    m_current(_init.config()),
    m_equiv_ind(0),
    m_wedges(_wedges),
    m_shape_factor(Eigen::MatrixXd::Identity(min_val.size(), min_val.size())) {

    //Condition range arrays and build shape_factor matrix
    Index nc = 0;
    if(m_wedges.size() == 0)
      return;
    Index nsub = m_wedges[0].irrep_wedges().size();
    for(Index s = 0; s < nsub; s++) {
      double abs_mag = 0.;
      for(Index i = 0; i < m_wedges[0].irrep_wedges()[s].axes.cols(); i++)
        abs_mag = max<double>(abs_mag, max(abs(max_val(nc)), abs(min_val(nc))));

      for(Index i = 0; i < m_wedges[0].irrep_wedges()[s].axes.cols(); i++, nc++) {
        if(m_wedges[0].irrep_wedges()[s].mult[i] == 1 && auto_range) {
          min_val(nc) = -max_val(nc);
        }

        if(abs_mag > TOL)
          m_shape_factor(nc, nc) /= abs_mag * abs_mag;

        if(inc_val(nc) == 0. || (max_val(nc) - min_val(nc)) / inc_val(nc) > 1e4) {
          max_val(nc) = min_val(nc) + TOL;
          inc_val(nc) = 10 * TOL;
        }
      }
    }

    //Find first valid config
    m_counter = EigenCounter<Eigen::VectorXd>(min_val, max_val, inc_val);

    m_counter.reset();

    //Increment past any invalid values, including those that are outside specified ellipsoid (if trim_corners==true)
    while(m_counter.valid() && (trim_corners && double(m_counter().transpose()*m_wedges[m_equiv_ind].trans_mat().transpose()*m_shape_factor * m_wedges[m_equiv_ind].trans_mat()*m_counter()) > 1.0 + TOL)) {
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
    //Increment past any invalid values, including those that are outside specified ellipsoid (if trim_corners==true)
    while(++m_counter && (m_trim_corners && double(m_counter().transpose()*m_wedges[m_equiv_ind].trans_mat().transpose()*m_shape_factor * m_wedges[m_equiv_ind].trans_mat()*m_counter()) > 1.0 + TOL)) {

    }

    // move to next part of wedge if necessary
    if(!m_counter.valid() && m_equiv_ind + 1 < m_wedges.size()) {
      m_counter.reset();
      ++m_equiv_ind;

      //Increment past any invalid values, including those that are outside specified ellipsoid (if trim_corners==true)
      // this time it's for the new wedge
      while(m_counter &&
            (m_trim_corners && double(m_counter().transpose()*m_wedges[m_equiv_ind].trans_mat().transpose()*m_shape_factor * m_wedges[m_equiv_ind].trans_mat()*m_counter()) > 1.0 + TOL)) {
        ++m_counter;
      }
    }

    if(m_counter.valid()) {

      m_current.configdof().set_global_dof(m_strain_key, m_wedges[m_equiv_ind].trans_mat() * m_counter());

      _increment_step();
    }
    else {
      _invalidate();
    }

    m_current.set_source(this->source(step()));
    return;
  }

}

