#include <boost/filesystem.hpp>
#include <boost/filesystem/fstream.hpp>

#include "casm/kinetics/DiffTransConfiguration_impl.hh"

#include "casm/app/DirectoryStructure.hh"
#include "casm/casm_io/VaspIO.hh"
#include "casm/symmetry/Orbit_impl.hh"
#include "casm/clex/Clexulator.hh"


namespace CASM {

  namespace Kinetics {

    // --- DiffTransConfiguration ---

    DiffTransConfiguration::DiffTransConfiguration(const Configuration &_from_config,
                                                   const PrimPeriodicDiffTransOrbit &_dtorbit) :
      m_config_A(_from_config),
      m_config_B(_from_config),
      m_sym_compare(_from_config.supercell()),
      m_diff_trans(m_sym_compare.prepare(_dtorbit.prototype())),
      m_orbit_name(_dtorbit.name()) {

      m_diff_trans.apply_to(m_config_B);
      _sort();
    }

    DiffTransConfiguration::DiffTransConfiguration(const Configuration &_from_config,
                                                   const DiffusionTransformation &_diff_trans) :
      m_config_A(_from_config),
      m_config_B(_from_config),
      m_sym_compare(_from_config.supercell()),
      m_diff_trans(m_sym_compare.prepare(_diff_trans)) {

      m_diff_trans.apply_to(m_config_B);
      _sort();
    }

    /// Construct a DiffTransConfiguration from JSON data
    DiffTransConfiguration::DiffTransConfiguration(
      const Supercell &_supercell,
      const jsonParser &_data) :

      m_config_A(_supercell),
      m_config_B(_supercell),
      m_sym_compare(_supercell),
      m_diff_trans(_supercell.prim()) {

      this->from_json(_data, _supercell);

    }

    /// Construct a DiffTransConfiguration from JSON data
    DiffTransConfiguration::DiffTransConfiguration(
      const PrimClex &_primclex,
      const jsonParser &_data) :

      m_config_A(Supercell(&_primclex, _primclex.prim().lattice())),
      m_config_B(Supercell(&_primclex, _primclex.prim().lattice())),
      m_sym_compare(supercell()),
      m_diff_trans(_primclex.prim()) {

      this->from_json(_data, _primclex);

    }

    const Supercell &DiffTransConfiguration::supercell() const {
      return from_config().supercell();
    }

    /// \brief Returns the initial configuration
    const Configuration &DiffTransConfiguration::from_config() const {
      return *m_from_config;
    }

    /// \brief Returns the final configuration
    const Configuration &DiffTransConfiguration::to_config() const {
      return *m_to_config;
    }

    /// \brief Returns the diffusion transformation that is occurring
    const DiffusionTransformation &DiffTransConfiguration::diff_trans() const {
      return m_diff_trans;
    }

    /// \brief Compare DiffTransConfiguration
    /// Compares m_diff_trans first then
    /// m_from_config if m_diff_trans are equal
    /// - Comparison is made using the sorted forms
    bool DiffTransConfiguration::operator<(const DiffTransConfiguration &B) const {
      return less()(B);
    }

    /// \brief Compare DiffTransConfiguration
    /// Compares m_diff_trans first then
    /// m_from_config if m_diff_trans are equal
    /// - Comparison is made using the sorted forms
    DiffTransConfigCompare DiffTransConfiguration::less() const {
      return DiffTransConfigCompare(*this);
    }

    /// \brief Check if DiffTransConfiguration are equal
    /// Compares m_diff_trans first then
    /// m_from_config if m_diff_trans are equal
    /// - Comparison is made using the sorted forms
    DiffTransConfigIsEqual DiffTransConfiguration::equal_to() const {
      return DiffTransConfigIsEqual(*this);
    }

    /// \brief sort DiffTransConfiguration in place
    DiffTransConfiguration &DiffTransConfiguration::sort() {
      // now always keeping sorted
      return *this;
    }

    /// \brief Returns a sorted version of this DiffTransConfiguration
    DiffTransConfiguration DiffTransConfiguration::sorted() const {
      // now always keeping sorted
      return *this;
    }

    bool DiffTransConfiguration::is_sorted() const {
      // now always keeping sorted
      return true;
    }

    /// \brief Apply symmetry, does not sort
    DiffTransConfiguration &DiffTransConfiguration::apply_sym(const PermuteIterator &it) {
      m_diff_trans = m_sym_compare.prepare(m_diff_trans.apply_sym(it));
      m_config_A.apply_sym(it);
      m_config_B = m_config_A;
      m_diff_trans.apply_to(m_config_B);
      _sort();
      return *this;
    }

    std::string DiffTransConfiguration::generate_name_impl() const {
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
    void DiffTransConfiguration::from_json(const jsonParser &json, const Supercell &_scel) {

      // get cache
      CASM::from_json(cache(), json["cache"]);

      // get config
      std::vector<std::string> splt_vec;
      std::string configname = json["from_configname"].get<std::string>();
      boost::split(splt_vec, configname, boost::is_any_of("/"), boost::token_compress_on);
      m_config_A = Configuration(_scel, splt_vec[1], json["from_config_data"]);
      m_config_B = m_config_A;

      // get diff trans
      m_sym_compare = ScelPeriodicDiffTransSymCompare(_scel);
      m_diff_trans = m_sym_compare.prepare(
                       jsonConstructor<Kinetics::DiffusionTransformation>::from_json(json["diff_trans"], prim()));
      m_diff_trans.apply_to(m_config_B);
      set_orbit_name(json["orbit_name"].get<std::string>());

      _sort();
    }

    /// Reads the DiffTransConfiguration from JSON
    void DiffTransConfiguration::from_json(const jsonParser &json, const PrimClex &_primclex) {

      // get cache
      CASM::from_json(cache(), json["cache"]);

      // get config
      m_config_A = Configuration(
                     _primclex,
                     json["from_configname"].get<std::string>(),
                     json["from_config_data"]);
      m_config_B = m_config_A;

      // get diff trans
      m_sym_compare = ScelPeriodicDiffTransSymCompare(supercell());
      m_diff_trans = m_sym_compare.prepare(
                       jsonConstructor<Kinetics::DiffusionTransformation>::from_json(json["diff_trans"], prim()));
      m_diff_trans.apply_to(m_config_B);
      set_orbit_name(json["orbit_name"].get<std::string>());

      _sort();
    }

    /*
    /// Reads the DiffTransConfiguration from JSON
    void DiffTransConfiguration::from_json(const jsonParser &json, const PrimClex &_primclex) {

      // get cache
      CASM::from_json(cache(), json["cache"]);

      // get from_config
      m_config_A = make_configuration(_primclex, json["from_config"].get<std::string>());
      m_config_B = m_config_A;

      // get diff_trans
      m_sym_compare = ScelPeriodicDiffTransSymCompare(supercell());
      m_diff_trans = m_sym_compare.prepare(
        make_diff_trans(_primclex, json["diff_trans"].get<std::string>()));

      _sort();
    }
    */

    /// Check that DiffTrans does not include a single supercell site twice due
    /// to small supercell size
    bool DiffTransConfiguration::is_valid() const {
      std::set<Index> unique_indices;
      for(auto &traj : diff_trans().specie_traj()) {
        Index l = from_config().supercell().linear_index(traj.from.uccoord);
        unique_indices.insert(l);
      }
      return (diff_trans().specie_traj().size() == unique_indices.size());
    }

    bool DiffTransConfiguration::has_valid_from_occ() const {
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
      auto cache_it = cache().find("from_configname");
      if(cache_it == cache().end()) {
        std::string result = from_config().canonical_form().name();
        cache_insert("from_configname", result);
        return result;
      }
      return cache_it->get<std::string>();
    }

    /// The name of the canonical form of the to config
    std::string DiffTransConfiguration::to_configname() const {
      auto cache_it = cache().find("to_configname");
      if(cache_it == cache().end()) {
        std::string result = to_config().canonical_form().name();
        cache_insert("to_configname", result);
        return result;
      }
      return cache_it->get<std::string>();
    }

    /// A permute iterator it such that from_config = copy_apply(it,from_config.canonical_form())
    PermuteIterator DiffTransConfiguration::from_config_from_canonical() const {
      if(!cache().contains("from_config_from_canonical")) {
        PermuteIterator result = from_config().from_canonical();
        cache_insert("from_config_from_canonical", result);
        return result;
      }
      return cache()["from_config_from_canonical"].get<PermuteIterator>(supercell());
    }

    /// A permute iterator it such that to_config = copy_apply(it,to_config.canonical_form())
    PermuteIterator DiffTransConfiguration::to_config_from_canonical() const {
      if(!cache().contains("to_config_from_canonical")) {
        PermuteIterator result = from_config().from_canonical();
        cache_insert("to_config_from_canonical", result);
        return result;
      }
      return cache()["to_config_from_canonical"].get<PermuteIterator>(supercell());
    }

    void DiffTransConfiguration::write_pos() const {
      const auto &dir = primclex().dir();
      try {
        fs::create_directories(dir.configuration_dir(name()));
      }
      catch(const fs::filesystem_error &ex) {
        std::cerr << "Error in DiffTransConfiguration::write_pos()." << std::endl;
        std::cerr << ex.what() << std::endl;
      }

      fs::ofstream file(dir.POS(name()));
      write_pos(file);
      return;
    }
    std::ostream &DiffTransConfiguration::write_pos(std::ostream &sout) const {
      sout << "Initial POS:" << std::endl;
      VaspIO::PrintPOSCAR from(sorted().from_config());
      from.print(sout);
      sout << std::endl;
      sout << "Final POS:" << std::endl;
      VaspIO::PrintPOSCAR to(sorted().to_config());
      to.print(sout);
      return sout;
    }

    /// \brief prints this DiffTransConfiguration
    std::ostream &operator<<(std::ostream &sout, const DiffTransConfiguration &dtc) {
      sout << dtc.diff_trans();
      sout << dtc.from_config();
      return sout;
    }

    void DiffTransConfiguration::_sort() {

      if(m_config_B < m_config_A) {
        m_from_config = &m_config_B;
        m_to_config = &m_config_A;
      }
      else {
        m_from_config = &m_config_A;
        m_to_config = &m_config_B;
      }
    }

    /// \brief returns a copy of bg_config with sites altered such that diff_trans can be placed as is
    Configuration make_attachable(const DiffusionTransformation &diff_trans, const Configuration &bg_config) {
      Configuration result = bg_config;
      ScelPeriodicDiffTransSymCompare sym_compare(bg_config.supercell());
      if(diff_trans != sym_compare.prepare(diff_trans)) {
        std::cout << "diff_trans: \n" << diff_trans << std::endl;
        std::cout << "sym_compare.prepare(diff_trans): \n" << sym_compare.prepare(diff_trans) << std::endl;

        throw std::runtime_error("Error in make_attachable(const DiffusionTransformation &diff_trans, const Configuration &bg_config)");
      }
      for(auto traj : diff_trans.specie_traj()) {
        Index l = bg_config.supercell().linear_index(traj.from.uccoord);
        if(bg_config.occ(l) != traj.from.occ) {
          if(traj.from.occ <= bg_config.supercell().max_allowed_occupation()[l]) {
            result.set_occ(l, traj.from.occ);
          }
          else {
            throw std::runtime_error("Error in make_attachable(const DiffusionTransformation &diff_trans, const Configuration &bg_config)");
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
