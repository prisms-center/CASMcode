#include "casm/monte_carlo/canonical/CanonicalSettings.hh"
#include "casm/monte_carlo/canonical/CanonicalConditions.hh"
#include "casm/monte_carlo/canonical/CanonicalIO.hh"
#include "casm/app/AppIO.hh"

namespace CASM {
  namespace Monte {
    std::string _help() {
      std::string s =
        "For CanonicalConditions, expect a JSON object of form:\n"
        "  {\n"
        "    \"comp\": {                  // option 1: parameteric composition object\n"
        "      \"a\" : 0.3,\n"
        "      ...\n"
        "    },\n"
        "    \"comp\": [0.3, 0.2, ...],   // option 2: parameteric composition array\n"
        "    \"comp_n\": {                // option 3: mol per prim composition object\n"
        "      \"A\" : 1.2,\n"
        "      ...\n"
        "    },\n"
        "    \"comp_n\": [1.2, 0.3, ...], // option 4: mol per prim composition array\n"
        "    \"temperature\" : 350.0,\n"
        "    \"tolerance\" : 0.001\n"
        "  }\n";
      return s;
    }

    /// \brief Construct EquilibriumMonteSettings by reading a settings JSON file
    CanonicalSettings::CanonicalSettings(const PrimClex &_primclex, const fs::path &read_path) :
      EquilibriumMonteSettings(_primclex, read_path) {

      if(!_primclex.has_composition_axes()) {
        throw std::runtime_error("No composition axes selected.");
      }

    }

    // --- CanonicalConditions settings ---------------------

    /// \brief Expects initial_conditions
    CanonicalConditions CanonicalSettings::initial_conditions() const {
      if(drive_mode() == Monte::DRIVE_MODE::INCREMENTAL) {
        return _conditions("initial_conditions");
      }
      else if(drive_mode() == Monte::DRIVE_MODE::CUSTOM) {
        return custom_conditions()[0];
      }
      else {
        throw std::runtime_error("ERROR: Invalid drive mode.");
      }
    }

    /// \brief Expects final_conditions
    CanonicalConditions CanonicalSettings::final_conditions() const {
      return _conditions("final_conditions");
    }

    /// \brief Expects incremental_conditions
    CanonicalConditions CanonicalSettings::incremental_conditions() const {
      return _conditions("incremental_conditions");
    }

    /// \brief Expects incremental_conditions
    std::vector<CanonicalConditions> CanonicalSettings::custom_conditions() const {
      std::string level1 = "driver";
      std::string level2 = "custom_conditions";

      try {
        std::vector<CanonicalConditions> cond;
        const jsonParser &json = (*this)[level1][level2];
        for(auto it = json.begin(); it != json.end(); ++it) {
          cond.push_back(_conditions(*it));
        }
        return cond;
      }
      catch(std::runtime_error &e) {
        Log &err_log = default_err_log();
        err_log.error<Log::standard>("Reading Monte Carlo settings");
        err_log << "Tried to read an array of CanonicalConditions from [\"" << level1 << "\"][\"" << level2 << "\"]" << std::endl;
        err_log << _help() << std::endl;
        throw e;
      }
    }

    // --- Project settings ---------------------

    /// \brief Get formation energy cluster expansion
    ClexDescription CanonicalSettings::formation_energy(const PrimClex &primclex) const {
      const ProjectSettings &set = primclex.settings();
      std::string level1 = "model";
      // deprecated
      if(_is_setting(level1, "clex")) {

        // expect "clex" is "formation_energy"

        std::vector<std::string> var {"clex", "calctype", "ref", "bset", "eci"};
        std::vector<std::string> help {
          "string\n  Names the cluster expansion to be used.\n",
          "string\n  Names the calctype to be used.\n",
          "string\n  Names the reference to be used.\n",
          "string\n  Names the basis set to be used.\n",
          "string\n  Names the ECI to be used.\n"
        };

        return ClexDescription(
                 _get_setting<std::string>(level1, var[0], help[0]),
                 _get_setting<std::string>(level1, var[0], help[0]),
                 _get_setting<std::string>(level1, var[1], help[1]),
                 _get_setting<std::string>(level1, var[2], help[2]),
                 _get_setting<std::string>(level1, var[3], help[3]),
                 _get_setting<std::string>(level1, var[4], help[4]));
      }

      std::string help = "(string, default='formation_energy')\n"
                         "  Names the formation_energy cluster expansion to be used.\n";

      std::string formation_energy = "formation_energy";
      if(_is_setting(level1, "formation_energy")) {
        formation_energy = _get_setting<std::string>(level1, "formation_energy", help);
      }

      if(!set.has_clex(formation_energy)) {
        Log &err_log = default_err_log();
        err_log.error<Log::standard>("Reading Monte Carlo settings");
        err_log << "Error reading [\"model\"][\"formation_energy\"]\n";
        err_log << "[\"model\"][\"formation_energy\"]: (string, optional, default='formation_energy')\n";
        err_log << "  Names the cluster expansion to be used for calculating formation_energy.\n";
        err_log << "No cluster expansion named: '" << formation_energy << "' exists.\n";
      }
      return set.clex(formation_energy);
    }

    // --- Sampler settings ---------------------

    /// \brief Return true if all correlations should be sampled
    bool CanonicalSettings::all_correlations() const {
      std::string level1 = "data";
      std::string level2 = "measurements";
      try {
        const jsonParser &json = (*this)[level1][level2];
        for(auto it = json.cbegin(); it != json.cend(); ++it) {
          if(it->contains("quantity") && (*it)["quantity"].get<std::string>() == "all_correlations") {
            return true;
          }
        }
        return false;
      }
      catch(std::runtime_error &e) {
        Log &err_log = default_err_log();
        err_log.error<Log::standard>("Reading Monte Carlo settings");
        err_log << "Error checking if 'all_correlations' should be sampled.\n";
        err_log << "Error reading settings at [\"" << level1 << "\"][\"" << level2 << "\"]\n";
        err_log << "See 'casm format --monte' for help.\n" << std::endl;
        throw e;
      }

    }

    CanonicalConditions CanonicalSettings::_conditions(std::string name) const {

      std::string level1 = "driver";
      std::string level2 = name;
      try {
        return _conditions((*this)[level1][level2]);
      }
      catch(std::runtime_error &e) {
        Log &err_log = default_err_log();
        err_log.error<Log::standard>("Reading Monte Carlo settings");
        err_log << "Error reading: " << name << std::endl;
        err_log << "Tried to construct CanonicalCondtions from [\"" << level1 << "\"][\"" << level2 << "\"]" << std::endl;
        err_log << _help() << std::endl;
        throw e;
      }
    }

    CanonicalConditions CanonicalSettings::_conditions(const jsonParser &json) const {
      CanonicalConditions result;
      from_json(result, primclex(), json);
      return result;
    }
  }
}

