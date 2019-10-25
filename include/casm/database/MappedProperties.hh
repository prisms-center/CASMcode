#ifndef CASM_MappedProperties
#define CASM_MappedProperties

#include <string>
#include <ctime>
#include "casm/CASM_global_Eigen.hh"
#include "casm/casm_io/EnumIO.hh"

namespace CASM {

  /// \brief Data structure stored in PropertiesDatabase
  ///
  /// - Use documentation here as reference for reading/writing 'properties.calc.json'
  ///   files.
  ///
  struct MappedProperties {

    /// The name of the config that was the starting point for relaxation
    /// If no relaxation, from == to
    std::string from;

    /// The name of the config that was the ending point of relaxation
    /// If no relaxation, from == to
    std::string to;

    std::time_t timestamp;

    /// Raw & unmapped calculated properties, as from properties.calc.json:
    ///
    /// For Configuration:
    /// - "relaxed_forces": shape=(Nx3), Final forces
    ///   on atoms.
    /// - "relaxed_lattice": Lattice, shape=(3x3), Final lattice
    ///   vectors, stored column vectors of (3x3) matrix
    /// - "relaxed_basis": Eigen::MatrixXd, shape=(Nx3), Final basis
    ///   coordinates.
    /// - "relaxed_energy": double, shape=(1x1), Final energy.
    /// - "relaxed_magmom": double, shape=(1x1) Final total magnetic moment.
    /// - "relaxed_mag_basis": shape=(Nx1), Final magnetic moment at each
    ///   basis site
    std::map<std::string, Eigen::MatrixXd> global;

    /// Parsed & mapped properties, as from ConfigMapperResult:
    ///
    /// For Configuration:
    /// - From ConfigMapperResult::relaxation_properties:
    ///   - 'lattice_deformation': double, shape=(1x1) , lattice mapping score
    ///   - 'basis_deformation': double, shape=(1x1), basis mapping score
    ///   - 'relaxation_deformation': shape=(3x3), also referred to as F:
    ///     cart_op*L_calc = F*L_ideal, L_calc being the lattice in the imported structure,
    ///     and L_ideal the lattice of the ideal mapped to Configuration, and
    ///     cart_op is a cartesian transformation matrix operation.
    ///   - 'relaxation_displacement': shape=(N,3) Gives the
    ///     the cartesian displacement of each Molecule reference position from the
    ///     ideal basis location. Remember, displacements are applied before
    ///     deformation.
    /// - From ConfigMapperResult::cart_op:
    ///   - "cart_op": shape=(3x3), the cartesian isometry that rotates the
    ///     imported structure onto the coordinate system of the ideal crystal.
    ///
    std::map<std::string, Eigen::MatrixXd> site;

    bool has_scalar(std::string const &_name) const;

    double const &scalar(std::string const &_name) const;

    double &scalar(std::string const &_name);
  };

  //jsonParser &to_json(const MappedProperties &score, jsonParser &json);

  //void from_json(MappedProperties &score, const jsonParser &json);


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

    struct Option {
      Option(Method _method = Method::minimum, std::string _name = "relaxed_energy") :
        Option(_method, _name, 0.) {}


      Option(Method _method, double _lattice_weight = 0.5) :
        Option(_method, "", _lattice_weight) {}

      /// Method for scoring
      Method method;

      /// Property name or configname used for scoring
      std::string name;

      double lattice_weight;

    private:
      Option(Method _method, std::string _name, double _lattice_weight);
    };


    /// \brief Default uses minimum relaxed_energy
    ScoreMappedProperties();

    explicit ScoreMappedProperties(Option _opt = Option(Method::minimum, "relaxed_energy"));

    double operator()(const MappedProperties &obj) const;

    bool validate(const MappedProperties &obj) const;

    bool operator==(const ScoreMappedProperties &B) const;

    bool operator!=(const ScoreMappedProperties &B) const;

    const Option &option() const;

  private:
    Option m_opt;
  };

  //jsonParser &to_json(const ScoreMappedProperties &score, jsonParser &json);

  //void from_json(ScoreMappedProperties &score, const jsonParser &json);
}

namespace CASM {
  ENUM_IO_DECL(CASM::ScoreMappedProperties::Method)
  ENUM_TRAITS(CASM::ScoreMappedProperties::Method)
}

#endif
