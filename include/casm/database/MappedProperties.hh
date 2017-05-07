#ifndef CASM_MappedProperties
#define CASM_MappedProperties

#include <string>
#include "casm/casm_io/jsonParser.hh"
#include "casm/casm_io/EnumIO.hh"

namespace CASM {
  namespace DB {

    struct MappedProperties {

      /// The name of the config that was the starting point for relaxation
      /// If no relaxation, from == to
      std::string from;

      /// The name of the config that was the ending point of relaxation
      /// If no relaxation, from == to
      std::string to;

      /// Raw & unmapped calculated properties, as from properties.calc.json
      jsonParser unmapped;

      /// Parsed & mapped properties, as from ConfigMapper::import_structure output
      jsonParser mapped;
    };

    jsonParser &to_json(const MappedProperties &score, jsonParser &json);

    void from_json(MappedProperties &score, const jsonParser &json);


    /// Compare via A.from < B.from
    ///
    /// - Enforces one calculated properties per starting configuration
    inline bool operator<(const MappedProperties &A, const MappedProperties &B) {
      return A.from < B.from;
    }


    /// Resolve mapping conflicts by 'scoring' the MappedProperties structure
    ///
    /// Options:
    /// \code
    /// {
    ///   "method" : "deformation_cost",
    ///   "lattice_weight" : number in range [0,1.0]
    /// }
    /// {
    ///   "method" : "minimum",
    ///   "property" : property name (i.e. "relaxed_energy")
    /// }
    /// {
    ///   "method" : "maximum",
    ///   "property" : property name (i.e. "some property")
    /// }
    /// {
    ///   "method" : "direct_selection",
    ///   "name" : configname to force as 'best' (i.e. "SCEL3_1_1_3_0_0_0/4")
    /// }
    /// \endcode
    ///
    class ScoreMappedProperties {

    public:

      enum class Method {deformation_cost, minimum, maximum, direct_selection};

      /// \brief Default uses minimum relaxed_energy
      ScoreMappedProperties();

      explicit ScoreMappedProperties(const jsonParser &_params);

      double operator()(const MappedProperties &obj) const;

      bool validate(const MappedProperties &obj) const;

      bool operator==(const ScoreMappedProperties &B) const;

      bool operator!=(const ScoreMappedProperties &B) const;

      const jsonParser &params() const;

    private:

      jsonParser m_params;

      Method m_method;

      std::string m_propname;
      double m_lattice_weight;
      std::string m_direct_selection_name;
    };

    jsonParser &to_json(const ScoreMappedProperties &score, jsonParser &json);

    void from_json(ScoreMappedProperties &score, const jsonParser &json);

  }
}

namespace CASM {
  ENUM_IO_DECL(CASM::DB::ScoreMappedProperties::Method)
  ENUM_TRAITS(CASM::DB::ScoreMappedProperties::Method)
}

#endif
