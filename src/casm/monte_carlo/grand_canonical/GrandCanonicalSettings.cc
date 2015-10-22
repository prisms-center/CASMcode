#include "casm/monte_carlo/grand_canonical/GrandCanonicalSettings.hh"
#include "casm/monte_carlo/grand_canonical/GrandCanonicalConditions.hh"
#include "casm/monte_carlo/grand_canonical/GrandCanonicalIO.hh"

namespace CASM {
  
  // --- GrandCanonicalConditions settings ---------------------
    
  /// \brief Expects initial_conditions
  GrandCanonicalConditions GrandCanonicalSettings::initial_conditions() const {
    return _conditions("initial_conditions");
  }
  
  /// \brief Expects final_conditions
  GrandCanonicalConditions GrandCanonicalSettings::final_conditions() const {
    return _conditions("final_conditions");
  }
  
  /// \brief Expects incremental_conditions
  GrandCanonicalConditions GrandCanonicalSettings::incremental_conditions() const {
    return _conditions("incremental_conditions");
  }

  // --- Project settings ---------------------
  
  /// \brief Given a settings jsonParser figure out what the project clex settings to use are:
  std::string GrandCanonicalSettings::clex() const {
    std::string level1 = "initialization";
    std::string level2 = "clex";
    try {
      return (*this)[level1][level2].get<std::string>();
    }

    catch(std::runtime_error &e) {
      std::cerr << "ERROR in GrandCanonicalSettings::clex" << std::endl;
      std::cerr << "No clex settings found" << std::endl;
      std::cerr << "Expected [\"" << level1 << "\"][\"" << level2 << "\"]" << std::endl;
      throw e;
    }
  }

  /// \brief Given a settings jsonParser figure out what the project bset settings to use are:
  std::string GrandCanonicalSettings::bset() const {
    std::string level1 = "initialization";
    std::string level2 = "bset";
    try {
      return (*this)[level1][level2].get<std::string>();
    }

    catch(std::runtime_error &e) {
      std::cerr << "ERROR in GrandCanonicalSettings::bset" << std::endl;
      std::cerr << "No bset settings found" << std::endl;
      std::cerr << "Expected [\"" << level1 << "\"][\"" << level2 << "\"]" << std::endl;
      throw e;
    }
  }

  /// \brief Given a settings jsonParser figure out what the project calctype settings to use are:
  std::string GrandCanonicalSettings::calctype() const {
    std::string level1 = "initialization";
    std::string level2 = "calctype";
    try {
      return (*this)[level1][level2].get<std::string>();
    }

    catch(std::runtime_error &e) {
      std::cerr << "ERROR in GrandCanonicalSettings::calctype" << std::endl;
      std::cerr << "No calctype settings found" << std::endl;
      std::cerr << "Expected [\"" << level1 << "\"][\"" << level2 << "\"]" << std::endl;
      throw e;
    }
  }

  /// \brief Given a settings jsonParser figure out what the project ref settings to use are:
  std::string GrandCanonicalSettings::ref() const {
    std::string level1 = "initialization";
    std::string level2 = "ref";
    try {
      return (*this)[level1][level2].get<std::string>();
    }

    catch(std::runtime_error &e) {
      std::cerr << "ERROR in GrandCanonicalSettings::ref" << std::endl;
      std::cerr << "No ref settings found" << std::endl;
      std::cerr << "Expected [\"" << level1 << "\"][\"" << level2 << "\"]" << std::endl;
      throw e;
    }

  }

  /// \brief Given a settings jsonParser figure out what the project eci settings to use are:
  std::string GrandCanonicalSettings::eci() const {
    std::string level1 = "initialization";
    std::string level2 = "eci";
    try {
      return (*this)[level1][level2].get<std::string>();
    }

    catch(std::runtime_error &e) {
      std::cerr << "ERROR in GrandCanonicalSettings::eci" << std::endl;
      std::cerr << "No eci settings found" << std::endl;
      std::cerr << "Expected [\"" << level1 << "\"][\"" << level2 << "\"]" << std::endl;
      throw e;
    }

  }

  GrandCanonicalConditions GrandCanonicalSettings::_conditions(std::string name) const {
    
    std::string level1 = "driver";
    std::string level2 = name;
    try {
      
      DirectoryStructure dir(root());
      CompositionAxes axes(dir.composition_axes(calctype(), ref()));
      CompositionConverter comp_converter;
      
      if(axes.has_current_axes) {
        comp_converter = axes.curr;
      }
      else {
        throw std::runtime_error("No composition axes selected.");
      }
      
      GrandCanonicalConditions result;
      from_json(result, comp_converter, (*this)[level1][level2]);
      return result;
    }
    catch(std::runtime_error &e) {
      std::cerr << "ERROR in GrandCanonicalSettings::" << name << std::endl;
      std::cerr << "Tried construct GrandCanonicalCondtions from [\"" << level1 << "\"][\"" << level2 << "\"]" << std::endl;
      throw e;
    }
  }
  
}

