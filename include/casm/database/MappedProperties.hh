#ifndef CASM_MappedProperties
#define CASM_MappedProperties

#include <string>
#include "casm/casm_io/jsonParser.hh"
#include "casm/casm_io/EnumIO.hh"

namespace CASM {
  namespace DB {

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

      /// Raw & unmapped calculated properties, as from properties.calc.json:
      ///
      /// For Configuration:
      /// - "atom_type": std::vector<std::string>, List of species names as in
      ///   the PRIM, as ordered in basis.
      /// - "atoms_per_type": Eigen::VectorXi, List of number of atoms of each type,
      ///   corresponds to entries in "atom_type".
      /// - "coord_mode": str, Coordinate mode for basis (always 'direct').
      /// - "relaxed_forces": Eigen::MatrixXd, shape=(N,3), Final forces
      ///   on atoms.
      /// - "relaxed_lattice": Lattice, shape=(3,3), Final lattice
      ///   vectors, stored as list [ [a0, a1, a2], [b0, b1, b2], [c0, c1, c2]].
      /// - "relaxed_basis": Eigen::MatrixXd, shape=(N,3), Final basis
      ///   coordinates.
      /// - "relaxed_energy": double, Final energy.
      /// - "relaxed_magmom": double, Final total magnetic moment.
      /// - "relaxed_mag_basis": Eigen::VectorXd, Final magnetic moment at each
      ///   basis site
      jsonParser unmapped;

      /// Parsed & mapped properties, as from ConfigMapperResult:
      ///
      /// For Configuration:
      /// - From ConfigMapperResult::relaxation_properties:
      ///   - 'lattice_deformation': double, lattice mapping score
      ///   - 'basis_deformation': double, basis mapping score
      ///   - 'volume_relaxation': double, V/V_ideal
      ///   - 'relaxation_deformation': Eigen::Matrix3d, also referred to as F:
      ///     cart_op*L_calc = F*L_ideal, L_calc being the lattice in the imported structure,
      ///     and L_ideal the lattice of the ideal mapped to Configuration, and
      ///     cart_op is a cartesian transformation matrix operation.
      /// - From ConfigMapperResult::best_assignment:
      ///   - "best_assignment": Permutation, Representations the permutation
      ///     that maps sites in the imported structure onto sites of the ideal crystal.
      ///     i.e. config[i] = imported_structure[best_assignment[i]]
      /// - From ConfigMapperResult::cart_op:
      ///   - "cart_op": Eigen::Matrix3d, the cartesian isometry that rotates the
      ///     imported structure onto the coordinate system of the ideal crystal.
      ///
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
