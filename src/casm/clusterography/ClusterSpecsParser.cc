#include "casm/clusterography/ClusterSpecsParser_impl.hh"

namespace CASM {

  OrbitBranchSpecsParser::OrbitBranchSpecsParser(jsonParser &_input, fs::path _path, bool _required) :
    KwargsParser(_input, _path, _required),
    max_branch(max_orbit_branch()) {}

  jsonParser::const_iterator OrbitBranchSpecsParser::branch(int branch_i) const {
    auto it = self.find(branch_to_string(branch_i));
    if(it == self.end()) {
      throw std::runtime_error(std::string("Error in OrbitBranchSpecsParser: ")
                               + "Attempted to access branch '" + branch_to_string(branch_i) + "'");
    }
    return it;
  }

  int OrbitBranchSpecsParser::branch_to_int(jsonParser::const_iterator it) const {
    return boost::lexical_cast<int>(it.name());
  }

  std::string OrbitBranchSpecsParser::branch_to_string(int branch_i) const {
    return boost::lexical_cast<std::string>(branch_i);
  }

  bool OrbitBranchSpecsParser::is_integer(jsonParser::const_iterator it) const {
    try {
      int branch = branch_to_int(it);
      (void) branch;
      return true;
    }
    catch(std::exception &e) {
      //pass
    }
    return false;
  }

  /// warn for non-integer 'branch'
  bool OrbitBranchSpecsParser::warn_non_integer_branch(
    jsonParser::const_iterator it,
    const std::set<std::string> &expected) {
    if(expected.find(it.name()) != expected.end()) {
      return true;
    }
    if(!is_integer(it)) {
      warning.insert(std::string("Warning: Ignoring '")
                     + it.name() + "', is not an integer.");
      return false;
    }
    return true;
  }

  /// add warning if option is unnecessary
  bool OrbitBranchSpecsParser::warn_unnecessary(
    jsonParser::const_iterator it,
    const std::set<std::string> &expected) {
    bool all_necessary = true;
    for(auto opt_it = it->begin(); opt_it != it->end(); ++opt_it) {
      if(expected.find(opt_it.name()) == expected.end()) {
        warning.insert(
          std::string("Warning: Ignoring branch '") + it.name()
          + "' setting '" + opt_it.name() + "' (it is unnecessary).");
        all_necessary = false;
      }
    }
    return all_necessary;
  }

  /// add warning if branch is unnecessary
  bool OrbitBranchSpecsParser::warn_unnecessary_branch(jsonParser::const_iterator it) {
    warning.insert(
      std::string("Warning: Ignoring branch '") + it.name()
      + "' settings (they are unnecessary).");
    return false;
  }

  jsonParser::const_iterator OrbitBranchSpecsParser::previous(jsonParser::const_iterator it) {
    int branch = branch_to_int(it);
    int previous = branch - 1;
    return self.find(branch_to_string(previous));
  }

  int OrbitBranchSpecsParser::max_orbit_branch() const {
    int res = 1; // always assume branch 1
    if(exists()) {
      for(auto it = self.begin(); it != self.end(); ++it) {
        if(is_integer(it)) {
          int curr = branch_to_int(it);
          if(curr > res) {
            res = curr;
          }
        }
      }
    }
    return res;
  }


  // --- PrimPeriodicOrbitBranchSpecsParser ---

  PrimPeriodicOrbitBranchSpecsParser::PrimPeriodicOrbitBranchSpecsParser(
    jsonParser &_input,
    fs::path _path,
    bool _required) :
    OrbitBranchSpecsParser(_input, _path, _required) {

    if(exists()) {
      if(!self.is_obj()) {
        error.insert(std::string("Error: ") + name() + " is not a JSON object");
        return;
      }

      for(auto it = self.begin(); it != self.end(); ++it) {
        if(!warn_non_integer_branch(it)) {
          continue;
        }
        int branch = branch_to_int(it);
        if(branch < 2) {
          warn_unnecessary_branch(it);
          continue;
        }
        if(branch > 1) {
          require_at<double>(fs::path {it.name()} / "max_length");
          warn_unnecessary(it, {"max_length"});
        }
        if(branch > 2) {
          require_nonincreasing<double>(it, "max_length");
        }
      }
    }
  }

  double PrimPeriodicOrbitBranchSpecsParser::max_length(int branch_i) const {
    if(branch_i < 2 || branch_i > max_branch) {
      throw std::invalid_argument(std::string("Error: ")
                                  + "Requested max_length for invalid branch: " + branch_to_string(branch_i));
    }
    auto it = branch(branch_i)->find("max_length");
    if(it == branch(branch_i)->end()) {
      throw std::runtime_error(std::string("Error: ") + "'max_length' not given for branch: " + branch_to_string(branch_i));
    }
    return it->get<double>();
  }


  // --- PrimPeriodicOrbitSpecsParser ---

  PrimPeriodicOrbitSpecsParser::PrimPeriodicOrbitSpecsParser(
    const PrimClex &_primclex,
    const SymGroup &_generating_grp,
    const PrimPeriodicSymCompare<IntegralCluster> &_sym_compare,
    jsonParser &_input,
    fs::path _path,
    bool _required):
    KwargsParser(_input, _path, _required),
    primclex(_primclex),
    custom_generators(_generating_grp, _sym_compare) {

    if(exists()) {

      if(!self.is_array()) {
        error.insert(std::string("Error: ") + name() + " is not a JSON array");
        return;
      }

      // for each custom orbit
      for(Index i = 0; i < self.size(); ++i) {
        const jsonParser &val = self[i];
        KwargsParser tmp {input, path / boost::lexical_cast<std::string>(i), true};
        tmp.warn_unnecessary({"coordinate_mode", "prototype", "sites", "include_subclusters"});
        this->insert(tmp);

        // read orbit generating cluster from specs
        typename OrbitType::Element input_cluster(primclex.prim());
        try {
          from_json(input_cluster, val, primclex.crystallography_tol());
        }
        catch(std::exception &e) {
          std::stringstream msg;
          msg << "Error in '" << name() << "'[" << i << "] reading cluster: "
              << e.what() << std::endl << val << std::endl;
          error.insert(msg.str());
          continue;
        }

        try {
          // check if subclusters should be included (yes by default)
          auto f_it = val.find("include_subclusters");
          if(f_it == val.end() ||
             (f_it != val.end() && f_it->get<bool>())) {
            insert_subcluster_generators(input_cluster, custom_generators, primclex.log());
          }
          else {
            custom_generators.insert(input_cluster);
          }
        }
        catch(std::exception &e) {
          std::stringstream msg;
          msg << "Error: in '" << name() << "'[" << i << "] constructing generators: "
              << e.what() << std::endl << val << std::endl;
          error.insert(msg.str());
          continue;
        }
      }
    }

  }


  // --- PrimPeriodicClustersByMaxLength ---

  std::string PrimPeriodicClustersByMaxLength::cspecs_help() {
    return
      "  cspecs: JSON object \n"
      "    Indicate clusters to enumerate all occupational diffusion transformations. The \n"
      "    JSON item \"cspecs\" should be a cspecs style initialization of cluster number and sizes.\n"
      "    See below.          \n\n";
  }

  PrimPeriodicClustersByMaxLength::PrimPeriodicClustersByMaxLength(
    const PrimClex &_primclex,
    const SymGroup &_generating_grp,
    const PrimPeriodicSymCompare<IntegralCluster> &_sym_compare,
    jsonParser &_input,
    fs::path _path,
    bool _required) :
    InputParser(_input, _path, _required) {

    if(exists()) {
      auto orbit_branch_specs_subparser = std::make_shared<PrimPeriodicOrbitBranchSpecsParser>(
                                            input, relpath("orbit_branch_specs"), false);
      this->insert(orbit_branch_specs_subparser->path, orbit_branch_specs_subparser);

      auto orbit_specs_subparser = std::make_shared<PrimPeriodicOrbitSpecsParser>(
                                     _primclex, _generating_grp, _sym_compare, input, relpath("orbit_specs"), false);
      this->insert(orbit_specs_subparser->path, orbit_specs_subparser);

      // require one of {"orbit_branch_specs", "orbit_specs"}
      if(!orbit_branch_specs_subparser->exists() && !orbit_specs_subparser->exists()) {
        error.insert(std::string("Error: ") + "One of \"orbit_branch_specs\" or \"orbit_specs\" is required.");
      }

      // warn about unnecessary properties
      warn_unnecessary({"orbit_branch_specs", "orbit_specs"});
    }
  }

  int PrimPeriodicClustersByMaxLength::max_branch() const {
    return orbit_branch_specs().max_branch;
  }

  double PrimPeriodicClustersByMaxLength::max_length(int branch_i) const {
    return orbit_branch_specs().max_length(branch_i);
  }

  const OrbitGenerators<PrimPeriodicOrbit<IntegralCluster>> &PrimPeriodicClustersByMaxLength::custom_generators() const {
    return orbit_specs().custom_generators;
  }

  const PrimPeriodicOrbitBranchSpecsParser &PrimPeriodicClustersByMaxLength::orbit_branch_specs() const {
    return static_cast<const PrimPeriodicOrbitBranchSpecsParser &>(*this->find(relpath("orbit_branch_specs"))->second);
  }

  const PrimPeriodicOrbitSpecsParser &PrimPeriodicClustersByMaxLength::orbit_specs() const {
    return static_cast<const PrimPeriodicOrbitSpecsParser &>(*this->find(relpath("orbit_specs"))->second);
  }


  // --- LocalOrbitBranchSpecsParser ---

  LocalOrbitBranchSpecsParser::LocalOrbitBranchSpecsParser(
    jsonParser &_input, fs::path _path, bool _required):
    OrbitBranchSpecsParser(_input, _path, _required) {

    if(exists()) {
      if(!self.is_obj()) {
        error.insert(std::string("Error: ") + name() + " is not a JSON object");
        return;
      }

      for(auto it = self.begin(); it != self.end(); ++it) {
        warn_non_integer_branch(it, {"max_length_including_phenomenal"});
        if(!is_integer(it)) {
          continue;
        }
        int branch = branch_to_int(it);
        if(branch < 1) {
          warn_unnecessary_branch(it);
          continue;
        }
        if(max_length_including_phenomenal()) {
          // branch 1+ requires cutoff_radius & max_length
          if(branch > 0) {
            require_at<double>(fs::path {it.name()} / "max_length");
            require_at<double>(fs::path {it.name()} / "cutoff_radius");
            warn_unnecessary(it, {"max_length", "cutoff_radius"});
          }
          if(branch > 1) {
            require_nonincreasing<double>(it, "max_length");
            require_nonincreasing<double>(it, "cutoff_radius");
          }
        }
        else {
          // branch 1+ requires cutoff_radius, 2+ requires max_length
          if(branch == 1) {
            require_at<double>(fs::path {it.name()} / "cutoff_radius");
            warn_unnecessary(it, {"cutoff_radius"});
          }
          else if(branch == 2) {
            require_at<double>(fs::path {it.name()} / "max_length");
            require_at<double>(fs::path {it.name()} / "cutoff_radius");
            require_nonincreasing<double>(it, "cutoff_radius");
            warn_unnecessary(it, {"max_length", "cutoff_radius"});
          }
          else if(branch > 2) {
            require_at<double>(fs::path {it.name()} / "max_length");
            require_at<double>(fs::path {it.name()} / "cutoff_radius");
            require_nonincreasing<double>(it, "max_length");
            require_nonincreasing<double>(it, "cutoff_radius");
            warn_unnecessary(it, {"max_length", "cutoff_radius"});
          }
        }
      }
    }
  }

  bool LocalOrbitBranchSpecsParser::max_length_including_phenomenal() const {
    bool result;
    self.get_else<bool>(result, "max_length_including_phenomenal", false);
    return result;
  }

  double LocalOrbitBranchSpecsParser::max_length(int branch_i) const {
    if(branch_i < 1 || branch_i > max_branch) {
      throw std::invalid_argument(std::string("Error: ")
                                  + "Requested max_length for invalid branch: " + branch_to_string(branch_i));
    }
    auto it = branch(branch_i)->find("max_length");
    if(it == branch(branch_i)->end()) {
      throw std::runtime_error(std::string("Error: ") + "'max_length' not given for branch: " + branch_to_string(branch_i));
    }
    return it->get<double>();
  }

  double LocalOrbitBranchSpecsParser::cutoff_radius(int branch_i) const {
    if(branch_i < 1 || branch_i > max_branch) {
      throw std::invalid_argument(std::string("Error: ")
                                  + "Requested cutoff_radius for invalid branch: " + branch_to_string(branch_i));
    }
    auto it = branch(branch_i)->find("cutoff_radius");
    if(it == branch(branch_i)->end()) {
      throw std::runtime_error(std::string("Error: ") + "'cutoff_radius' not given for branch: " + branch_to_string(branch_i));
    }
    return it->get<double>();
  }


  // --- LocalOrbitSpecsParser ---

  LocalOrbitSpecsParser::LocalOrbitSpecsParser(
    const PrimClex &_primclex,
    jsonParser &_input,
    fs::path _path,
    bool _required):
    KwargsParser(_input, _path, _required),
    primclex(_primclex) {

    if(exists()) {
      if(!self.is_array()) {
        error.insert(std::string("Error: ") + name() + " is not a JSON array");
        return;
      }

      // for each custom orbit
      std::set<std::string> opt {"coordinate_mode", "prototype", "sites", "include_subclusters"};
      for(Index i = 0; i < self.size(); ++i) {
        const jsonParser &val = self[i];
        KwargsParser tmp {input, path / boost::lexical_cast<std::string>(i), true};
        tmp.warn_unnecessary(opt);
        this->insert(tmp);

        // read orbit generating cluster from specs
        IntegralCluster input_cluster(primclex.prim());
        try {
          from_json(input_cluster, val, primclex.crystallography_tol());
        }
        catch(std::exception &e) {
          std::stringstream msg;
          msg << "Error in '" << name() << "'[" << i << "] reading cluster: "
              << e.what() << std::endl << val << std::endl;
          error.insert(msg.str());
          continue;
        }

        bool include_subclusters;
        try {
          // check if subclusters should be included (yes by default)
          val.get_else(include_subclusters, "include_subclusters", true);
        }
        catch(std::exception &e) {
          std::stringstream msg;
          msg << "Error: in '" << name() << "'[" << i << "] constructing generators: "
              << e.what() << std::endl << val << std::endl;
          error.insert(msg.str());
          continue;
        }

        prototypes.emplace_back(Data{input_cluster, include_subclusters});
      }
    }
  }
}
