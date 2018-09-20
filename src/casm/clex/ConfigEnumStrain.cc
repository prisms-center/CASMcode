#include "casm/clex/ConfigEnumStrain.hh"
#include <algorithm>
#include "casm/crystallography/Structure.hh"
#include "casm/clex/PrimClex.hh"
#include "casm/clex/Supercell.hh"
#include "casm/clex/ConfigEnumStrain.hh"
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

      "  config: string (required) Name of configuration used as strain reference\n"
      "    Must be single configuration. Strains are enumerated as perturbations to\n"
      "    specified configuration. Ex: \"config\" : \"SCEL2_2_1_1_0_0_0\"\n\n"

      "  min: array of doubles (required) Minimum, starting value of grid counter\n"
      "    Dimension must be either 6, or equal to number of rows of \"axes\"\n"
      "    Ex: \"min\" : [-0.1, -0.1, -0.1, -0.01, -0.01, -0.01]\n\n"

      "  max: array of doubles (required) Maximum, final value of grid counter\n"
      "    Dimension must be either 6, or equal to number of rows of \"axes\"\n"
      "    Ex: \"min\" : [0.1, 0.1, 0.1, 0.01, 0.01, 0.01]\n\n"

      "  increment: array of doubles (required) Amount by which to increment counter elements\n"
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

    if(!analysis) {
      if(!(_kwargs.contains("min") && _kwargs.contains("max") && _kwargs.contains("increment")))
        throw std::runtime_error("JSON options for enumeration method 'ConfigEnumStrain' must specify ALL of \"min\", \"max\", and \"increment\"");
      from_json(min_val, _kwargs["min"]);
      from_json(max_val, _kwargs["max"]);
      from_json(inc_val, _kwargs["increment"]);
    }

    if(!_kwargs.contains("config") || !_kwargs["config"].is_string())
      throw std::runtime_error("JSON options for enumeration method 'ConfigEnumStrain' must include exactly one starting configuration, specified by field \"config\"");

    std::vector<std::string> filter_expr = make_enumerator_filter_expr(_kwargs, enum_opt);

    std::string configname = _kwargs["config"].get<std::string>();
    auto it primclex.const_db().find(configname);
    if(it == primclex.const_db().end())
      throw std::runtime_error("In ConfigEnumStrain::run(), specified configuration " + configname + " does not exist in database.\n");

    return run(primclex,
               *it,
               axes,
               min_val,
               max_val,
               inc_val,
               sym_axes,
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
                            bool trim_corners,
                            bool analysis,
                            std::vector<std::string> const &_filter_expr,
                            bool dry_run) {

    std::vector<Eigen::MatrixXd> wedges;
    std::vector<int> dims;
    SymGroup pg = _config.point_group();
    if(!analysis && ! sym_axes)
      wedges.push_back(_axes);
    else {
      auto result = SymRepTools::symrep_wedge_facets(_config.point_group(), _primclex.strainrep_ID());
      wedges = result.second;
      dims = result.first;
    }

    if(analysis) {
      Log &log = primclex.log();
      //PRINT INFO TO LOG:

      return 0;
    }
    auto constructor = [&](const Supercell & scel) {
      return notstd::make_unique<ConfigEnumStrain>(_config,
                                                   wedges,
                                                   min_val,
                                                   max_val,
                                                   inc_val,
                                                   trim_corners,);
    };

    int returncode = insert_unique_canon_configs(enumerator_name,
                                                 _primclex,
                                                 begin,
                                                 end,
                                                 lambda,
                                                 filter_expr,
                                                 dry_run);



  }

  ConfigEnumStrain::ConfigEnumStrain(const Configuration &_init,
                                     const std::vector<Eigen::MatrixXd> &_wedges,
                                     Eigen::Ref<const Eigen::VectorXd> const &min_val,
                                     Eigen::Ref<const Eigen::VectorXd> const &max_val,
                                     Eigen::Ref<const Eigen::VectorXd> const &inc_val,
                                     bool trim_corners) :
    m_current(_init),
    m_equiv_ind(0),
    m_strain_calc(_mode),
    m_perm_begin(_scel.permute_begin()),
    m_perm_end(_scel.permute_end()),
    m_shape_factor(Eigen::MatrixXd::Identity(m_strain_calc.dim(), m_strain_calc.dim())) {

    //Condition range arrays and build shape_factor matrix
    Index nc = 0;
    for(Index s = 0; s < num_sub; s++) {
      double wedgevol = sqrt((wedges[s].transpose() * wedges[s]).determinant());
      Index N = round(pow(linear_partitions[s], wedges[s].cols()));
      //double density = double(N) / pow(magnitudes[s], wedges[s].cols());
      //std::cout << "wedgevol: " << wedgevol << ", N: " << N;// << ", density: " << density;
      N = max(1, (int) ceil(pow(wedgevol * double(N), 1.0 / double(wedges[s].cols())) - TOL));
      //std::cout << ", linearN: " << N << "\n";
      //std::cout << "mult.size() is: " << mult.size() << "  and mult is ";
      //for(auto m : mult)
      //std::cout  << m << "  ";
      //std::cout << "\n";

      for(Index i = 0; i < wedges[s].cols(); i++, nc++) {
        if(mult[nc] == 1 && linear_partitions[s] > 1) {
          init(nc) = -absmags[s];
          //inc(nc)=2*absmags[s]/double(subspace_partitions[s]);
        }
        else {
          init(nc) = 0.0;
          //inc(nc)=absmags[s]/double(subspace_partitions[s]);
        }

        if(absmags[s] > TOL)
          m_shape_factor(nc, nc) /= absmags[s] * absmags[s];

        if(linear_partitions[s] < 2) {
          final(nc) = init(nc);
          inc(nc) = 10 * TOL;
        }
        else {
          final(nc) = absmags[s] + TOL;
          inc(nc) = absmags[s] / double(N);
        }
      }
    }



    //m_shape_factor = m_strain_calc.sop_transf_mat().transpose() * m_shape_factor * m_strain_calc.sop_transf_mat();
    //Find first valid config
    m_counter = EigenCounter<Eigen::VectorXd>(init, final, inc);

    m_counter.reset();
    while(m_counter.valid() && double(m_counter().transpose()*m_trans_mats[m_equiv_ind].transpose()*m_shape_factor * m_trans_mats[m_equiv_ind]*m_counter()) > 1.0 + TOL) {
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

    while(++m_counter && double(m_counter().transpose()*m_trans_mats[m_equiv_ind].transpose()*m_shape_factor * m_trans_mats[m_equiv_ind]*m_counter()) > 1.0 + TOL) {
      //just burning throught the count
    }

    // move to next part of wedge if necessary
    if(!m_counter.valid() && m_equiv_ind + 1 < m_trans_mats.size()) {
      m_counter.reset();
      ++m_equiv_ind;
      std::cout << "INCREMENTING m_equiv_ind to " << m_equiv_ind << "\n";
    }

    while(m_counter && double(m_counter().transpose()*m_trans_mats[m_equiv_ind].transpose()*m_shape_factor * m_trans_mats[m_equiv_ind]*m_counter()) > 1.0 + TOL) {
      //just burning throught the count
      ++m_counter;
    }

    if(m_counter.valid()) {
      m_current.set_deformation(m_strain_calc.unrolled_strain_metric_to_F(m_trans_mats[m_equiv_ind] * m_counter()));
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

