#include "casm/kinetics/DiffusionTransformation.hh"
#include "casm/clex/Configuration.hh"
#include "casm/clex/Supercell.hh"
#include "casm/symmetry/Orbit_impl.hh"
#include "casm/kinetics/DiffTransConfiguration.hh"
#include "casm/clex/Clexulator.hh"
#include "casm/clex/PrimClex.hh"
#include "casm/clex/Supercell.hh"
#include "casm/clex/NeighborList.hh"
#include "casm/database/ConfigDatabase.hh"

namespace CASM {

  namespace Kinetics {


    DiffTransConfiguration::DiffTransConfiguration(const Configuration &_from_config,
                                                   const DiffusionTransformation &_diff_trans) :
      m_diff_trans(_diff_trans), m_from_config(_from_config), Calculable(_from_config.primclex()) {
      ScelPeriodicDiffTransSymCompare symcompare(m_from_config.supercell().prim_grid(),
                                                 m_from_config.supercell().crystallography_tol());
      m_diff_trans = symcompare.prepare(m_diff_trans);
    }

    /// Construct a DiffTransConfiguration from JSON data
    DiffTransConfiguration::DiffTransConfiguration(
      const Supercell &_supercell,
      const jsonParser &_data) :
      m_diff_trans(_supercell.prim()),
      m_from_config(_supercell),
      Calculable(_supercell.primclex()) {

      this->from_json(_data, _supercell);

    }

    /// Construct a DiffTransConfiguration from JSON data
    DiffTransConfiguration::DiffTransConfiguration(
      const PrimClex &_primclex,
      const jsonParser &_data) :
      m_diff_trans(_primclex.prim()),
      m_from_config(Supercell(& _primclex, _primclex.prim().lattice())),
      Calculable(_primclex) {

      this->from_json(_data, _primclex);

    }

    /// \brief sort DiffTransConfiguration in place
    DiffTransConfiguration &DiffTransConfiguration::sort() {
      Configuration to = to_config();
      if(to < m_from_config) {
        m_from_config = to;
        m_diff_trans.reverse();
      }
      return *this;
    }

    /// \brief Returns a sorted version of this DiffTransConfiguration
    DiffTransConfiguration DiffTransConfiguration::sorted() const {
      DiffTransConfiguration tmp {*this};
      return tmp.sort();
    }

    bool DiffTransConfiguration::is_sorted() const {
      Configuration to = to_config();
      return m_from_config < to;
    }

    PermuteIterator DiffTransConfiguration::to_canonical() const {
      // check which supercell factor group operations
      // when applied to m_diff_trans results in the greatest
      // DiffusionTransformation
      std::vector<PermuteIterator> checklist;
      ScelPeriodicDiffTransSymCompare symcompare(m_from_config.supercell().prim_grid(),
                                                 m_from_config.supercell().crystallography_tol());
      DiffusionTransformation greatest = symcompare.prepare(m_diff_trans);
      for(auto it = m_from_config.supercell().permute_begin();
          it != m_from_config.supercell().permute_end(); ++it) {
        DiffusionTransformation tmp = symcompare.prepare(copy_apply(it.sym_op(), m_diff_trans));
        if(tmp == greatest) {
          checklist.push_back(it);
        }
        else if(tmp > greatest) {
          checklist.clear();
          greatest = tmp;
          checklist.push_back(it);
        }
      }

      // of these operations check which one maximizes
      // the result of applying to m_from_config
      auto it = checklist.begin();
      DiffTransConfiguration max_dtc(copy_apply(*it, from_config()), greatest);
      PermuteIterator canon_op_it {*it};
      ++it;
      for(; it != checklist.end(); ++it) {

        Configuration tmp = copy_apply(*it, from_config());

        DiffTransConfiguration dtc_tmp(tmp, greatest);
        if((dtc_tmp >= max_dtc || !max_dtc.has_valid_from_occ()) && dtc_tmp.has_valid_from_occ()) {
          max_dtc = dtc_tmp.sorted();
          canon_op_it = *it;
        }
      }
      // return the operation that transforms this to canonical form
      return canon_op_it;
    }

    DiffTransConfiguration DiffTransConfiguration::canonical_form() const {
      return copy_apply(this->to_canonical(), *this);
    }

    bool DiffTransConfiguration::is_canonical() const {
      return std::all_of(m_from_config.supercell().permute_begin(),
                         m_from_config.supercell().permute_end(),
      [&](const PermuteIterator & p) {
        return (copy_apply(p, *this) <= *this);
      });
    }

    DiffTransConfiguration &DiffTransConfiguration::apply_sym(const PermuteIterator &it) {
      ScelPeriodicDiffTransSymCompare symcompare(m_from_config.supercell().prim_grid(),
                                                 m_from_config.supercell().crystallography_tol());
      m_from_config = apply(it, m_from_config);

      m_diff_trans = apply(it, m_diff_trans);

      m_diff_trans = symcompare.prepare(m_diff_trans);

      return *this;
    }

    std::string DiffTransConfiguration::_generate_name() const {
      return m_orbit_name + "/" + from_config().supercell().name() + "/" + id();
    }

    void DiffTransConfiguration::set_orbit_name(const std::string &orbit_name) {
      m_orbit_name = orbit_name;
    }




    /// Writes the DiffTransConfiguration to JSON
    jsonParser &DiffTransConfiguration::to_json(jsonParser &json) const {
      json.put_obj();
      json["from_configname"] = from_config().name();
      from_config().to_json(json["from_config_data"]);
      CASM::to_json(diff_trans(), json["diff_trans"]);
      json["orbit_name"] = orbit_name();
      json["cache"].put_obj();
      if(cache_updated()) {
        json["cache"] = cache();
      }
      return json;
    }

    /// Reads the DiffTransConfiguration from JSON
    void DiffTransConfiguration::from_json(const jsonParser &json, const Supercell &scel) {
      m_diff_trans = jsonConstructor<Kinetics::DiffusionTransformation>::from_json(json["diff_trans"], scel.prim());
      std::vector<std::string> splt_vec;
      std::string configname = json["from_configname"].get<std::string>();
      boost::split(splt_vec, configname, boost::is_any_of("/"), boost::token_compress_on);
      m_from_config = Configuration(scel, splt_vec[1], json["from_config_data"]);
      set_orbit_name(json["orbit_name"].get<std::string>());
    }

    /// Reads the DiffTransConfiguration from JSON
    void DiffTransConfiguration::from_json(const jsonParser &json, const PrimClex &primclex) {
      m_diff_trans = jsonConstructor<Kinetics::DiffusionTransformation>::from_json(json["diff_trans"], primclex.prim());
      m_from_config = Configuration(primclex, json["from_configname"].get<std::string>(), json["from_config_data"]);
      set_orbit_name(json["orbit_name"].get<std::string>());
    }

    bool DiffTransConfiguration::is_valid_neb() const {
      std::set<Index> unique_indeces;
      for(auto &traj : diff_trans().specie_traj()) {
        Index l = from_config().supercell().linear_index(traj.from.uccoord);
        unique_indeces.insert(l);
      }
      return (diff_trans().specie_traj().size() == unique_indeces.size());
    }

    bool DiffTransConfiguration::has_valid_from_occ() const {
      ScelPeriodicDiffTransSymCompare symcompare(from_config().supercell().prim_grid(),
                                                 from_config().supercell().crystallography_tol());
      if(diff_trans() != symcompare.prepare(diff_trans())) {
        std::cerr << "Diffusion Transformation is not based in this Configuration's supercell!" << std::endl;
      }
      for(auto traj : diff_trans().specie_traj()) {
        Index l = from_config().supercell().linear_index(traj.from.uccoord);
        if(from_config().occ(l) != traj.from.occ) {
          return false;
        }
      }
      return true;
    }


    /// The name of the canonical form of the from config
    std::string DiffTransConfiguration::from_configname() const {
      auto it = primclex().db<Configuration>().scel_range(from_config().supercell().name()).begin();
      for(; it != primclex().db<Configuration>().scel_range(from_config().supercell().name()).end(); ++it) {
        if(it->is_equivalent(from_config())) {
          return it->name();
        }
      }
      return "not_in_project";
    }

    /// The name of the canonical form of the to config
    std::string DiffTransConfiguration::to_configname() const {
      auto it = primclex().db<Configuration>().scel_range(from_config().supercell().name()).begin();
      for(; it != primclex().db<Configuration>().scel_range(from_config().supercell().name()).end(); ++it) {
        if(it->is_equivalent(to_config())) {
          return it->name();
        }
      }
      return "not_in_project";
    }

    /// A permute iterator it such that from_config = copy_apply(it,from_config.canonical_form())
    PermuteIterator DiffTransConfiguration::from_config_from_canonical() const {
      return from_config().from_canonical();
    }

    /// A permute iterator it such that to_config = copy_apply(it,to_config.canonical_form())
    PermuteIterator DiffTransConfiguration::to_config_from_canonical() const {
      return to_config().from_canonical();
    }

    /// \brief prints this DiffTransConfiguration
    std::ostream &operator<<(std::ostream &sout, const DiffTransConfiguration &dtc) {
      sout << dtc.diff_trans();
      sout << dtc.from_config();
      return sout;
    }

    /// \brief returns a copy of bg_config with sites altered such that diff_trans can be placed as is
    Configuration make_attachable(const DiffusionTransformation &diff_trans, const Configuration &bg_config) {
      Configuration result = bg_config;
      ScelPeriodicDiffTransSymCompare symcompare(bg_config.supercell().prim_grid(),
                                                 bg_config.supercell().crystallography_tol());
      if(diff_trans != symcompare.prepare(diff_trans)) {
        std::cerr << "Diffusion Transformation is not based in this Configuration's supercell!" << std::endl;
      }
      for(auto traj : diff_trans.specie_traj()) {
        Index l = bg_config.supercell().linear_index(traj.from.uccoord);
        if(bg_config.occ(l) != traj.from.occ) {
          if(traj.from.occ <= bg_config.supercell().max_allowed_occupation()[l]) {
            result.set_occ(l, traj.from.occ);
          }
          else {
            std::cerr << "Diffusion Transformation does not have valid starting occupant for this position in this Configuration" << std::endl;
          }
        }
      }
      return result;
    }


    /// \brief Returns correlations using 'clexulator'.
    Eigen::VectorXd correlations(const DiffTransConfiguration &dtc, Clexulator &clexulator) {
      Eigen::VectorXd correlations = Eigen::VectorXd::Zero(clexulator.corr_size());

      //do work here

      return correlations;
    }

  }

  Kinetics::DiffTransConfiguration jsonConstructor<Kinetics::DiffTransConfiguration>::from_json(
    const jsonParser &json,
    const PrimClex &primclex) {

    return Kinetics::DiffTransConfiguration(primclex, json);
  }

  Kinetics::DiffTransConfiguration jsonConstructor<Kinetics::DiffTransConfiguration>::from_json(
    const jsonParser &json,
    const Supercell &scel) {

    return Kinetics::DiffTransConfiguration(scel, json);
  }
}
