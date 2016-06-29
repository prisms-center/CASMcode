#include "casm/monte_carlo/grand_canonical/GrandCanonicalSettings.hh"
#include "casm/monte_carlo/grand_canonical/GrandCanonicalConditions.hh"
#include "casm/monte_carlo/grand_canonical/GrandCanonicalIO.hh"
#include "casm/app/AppIO.hh"

namespace CASM {

  std::string _help() {
    std::string s =
      "Expected JSON object of form:\n"
      "  {\n"
      "    \"param_chem_pot\": {\n"
      "    \"a\" : -1.0,\n"
      "    ...\n"
      "  },\n"
      "  \"temperature\" : 350.0,\n"
      "  \"tolerance\" : 0.001\n"
      "}\n";
    return s;
  }

  /// \brief Construct EquilibriumMonteSettings by reading a settings JSON file
  GrandCanonicalSettings::GrandCanonicalSettings(const fs::path &read_path) :
    EquilibriumMonteSettings(read_path) {

    DirectoryStructure dir(root());
    CompositionAxes axes(dir.composition_axes(calctype(), ref()));

    if(axes.has_current_axes) {
      m_comp_converter = axes.curr;
    }
    else {
      throw std::runtime_error("No composition axes selected.");
    }

  }

  // --- GrandCanonicalConditions settings ---------------------

  /// \brief Expects initial_conditions
  GrandCanonicalConditions GrandCanonicalSettings::initial_conditions() const {
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
  GrandCanonicalConditions GrandCanonicalSettings::final_conditions() const {
    return _conditions("final_conditions");
  }

  /// \brief Expects incremental_conditions
  GrandCanonicalConditions GrandCanonicalSettings::incremental_conditions() const {
    return _conditions("incremental_conditions");
  }

  /// \brief Expects incremental_conditions
  std::vector<GrandCanonicalConditions> GrandCanonicalSettings::custom_conditions() const {
    std::string level1 = "driver";
    std::string level2 = "custom_conditions";
    try {
      std::vector<GrandCanonicalConditions> cond;
      const jsonParser &json = (*this)[level1][level2];
      for(auto it = json.begin(); it != json.end(); ++it) {
        cond.push_back(_conditions(*it));
      }
      return cond;
    }
    catch(std::runtime_error &e) {
      std::cerr << "ERROR in GrandCanonicalSettings::custom_conditions" << std::endl;
      std::cerr << "Tried to construct GrandCanonicalCondtions from [\"" << level1 << "\"][\"" << level2 << "\"]" << std::endl;
      std::cerr << _help();
      throw e;
    }
  }

  // --- Project settings ---------------------

  /// \brief Get formation energy cluster expansion
  ClexDescription GrandCanonicalSettings::formation_energy(const PrimClex &primclex) const {
    if(_is_setting("model", "formation_energy") {
    primclex.clex(_get_setting<std::string>(level1, "model", "formation_energy"));
    }
    else {
      // deprecated
      std::vector<std::string> var {"clex", "calctype", "ref", "bset", "eci"};
      std::vector<std::string> help {
        "string\n  Names the cluster expansion to be used.\n",
        "string\n  Names the calctype to be used.\n",
        "string\n  Names the reference to be used.\n",
        "string\n  Names the basis set to be used.\n",
        "string\n  Names the ECI to be used.\n",

        return ClexDescription(_get_setting<std::string>(level1, var[0], help[0],
        _get_setting<std::string>(level1, var[1], help[1],
        _get_setting<std::string>(level1, var[2], help[2],
        _get_setting<std::string>(level1, var[3], help[3],
        _get_setting<std::string>(level1, var[4], help[4],

      }
    }

    /// \brief Given a settings jsonParser figure out what the project clex settings to use are:
    std::string GrandCanonicalSettings::clex() const {
      std::string level1 = "model";
      std::string level2 = "clex";
      std::string help = "string\n"
                         "  Names the cluster expansion to be used.\n";
      return _get_setting<std::string>(level1, level2, help);
    }

    /// \brief Given a settings jsonParser figure out what the project bset settings to use are:
    std::string GrandCanonicalSettings::bset() const {
      std::string level1 = "model";
      std::string level2 = "bset";
      std::string help = "string\n"
                         "  Names the basis set to be used.\n";
      return _get_setting<std::string>(level1, level2, help);
    }

    /// \brief Given a settings jsonParser figure out what the project calctype settings to use are:
    std::string GrandCanonicalSettings::calctype() const {
      std::string level1 = "model";
      std::string level2 = "calctype";
      std::string help = "string\n"
                         "  Names the calctype used.\n";
      return _get_setting<std::string>(level1, level2, help);
    }

    /// \brief Given a settings jsonParser figure out what the project ref settings to use are:
    std::string GrandCanonicalSettings::ref() const {
      std::string level1 = "model";
      std::string level2 = "ref";
      std::string help = "string\n"
                         "  Names the reference used.\n";
      return _get_setting<std::string>(level1, level2, help);
    }

    /// \brief Given a settings jsonParser figure out what the project eci settings to use are:
    std::string GrandCanonicalSettings::eci() const {
      std::string level1 = "model";
      std::string level2 = "eci";
      std::string help = "string\n"
                         "  Names the ECI to be used\n";
      return _get_setting<std::string>(level1, level2, help);
    }

    /// \brief Get current clex key
    ClexKey GrandCanonicalSettings::clex_key() const {
      return ClexKey(clex_name(), calctype(), ref(), bset(), eci());
    }

    // --- Sampler settings ---------------------

    /// \brief Return true if all correlations should be sampled
    bool GrandCanonicalSettings::all_correlations() const {
      if(method() == Monte::METHOD::LTE1) { //hack
        return false;
      }
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
        std::cerr << "ERROR in GrandCanonicalSettings::all_correlations" << std::endl;
        std::cerr << "Expected [\"" << level1 << "\"][\"" << level2 << "\"]" << std::endl;
        throw e;
      }

    }

    GrandCanonicalConditions GrandCanonicalSettings::_conditions(std::string name) const {

      std::string level1 = "driver";
      std::string level2 = name;
      try {
        return _conditions((*this)[level1][level2]);
      }
      catch(std::runtime_error &e) {
        std::cerr << "ERROR in GrandCanonicalSettings::" << name << std::endl;
        std::cerr << "Tried to construct GrandCanonicalCondtions from [\"" << level1 << "\"][\"" << level2 << "\"]" << std::endl;
        std::cerr << _help();
        throw e;
      }
    }

    GrandCanonicalConditions GrandCanonicalSettings::_conditions(const jsonParser & json) const {
      GrandCanonicalConditions result;
      from_json(result, m_comp_converter, json);
      return result;
    }

  }

