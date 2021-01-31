#ifndef CASM_AppIO
#define CASM_AppIO

#include <boost/filesystem/path.hpp>
#include <map>
#include <memory>
#include <string>

#include "casm/casm_io/Log.hh"
#include "casm/casm_io/json/jsonParser.hh"
#include "casm/clex/CompositionConverter.hh"
#include "casm/crystallography/CoordinateSystems.hh"
#include "casm/global/definitions.hh"
#include "casm/global/enum.hh"

namespace CASM {
namespace xtal {
class COORD_MODE;
class SpeciesAttribute;
class AtomPosition;
class Molecule;
class Site;
class UnitCellCoord;
class BasicStructure;
}  // namespace xtal
using xtal::AtomPosition;
using xtal::BasicStructure;
using xtal::COORD_MODE;
using xtal::Molecule;
using xtal::Site;
using xtal::SpeciesAttribute;
using xtal::UnitCellCoord;

// --- These functions are for casm I/O -----------

class jsonParser;
class HamiltonianModules;
class IntegralCluster;
class ChemicalReference;
class SymGroup;
class Structure;

/** \defgroup ProjectIO
 *
 *  \ingroup Project
 *  \ingroup casmIO
 *
 *  \brief Relates to CASM project input/output
 *
 *  @{
 */

// --------- PrimIO Declarations
// --------------------------------------------------

/// \brief Read SpeciesAttribute from json
jsonParser const &from_json(SpeciesAttribute &_attr, jsonParser const &json);

/// \brief From SpeciesAttribute to json
jsonParser &to_json(SpeciesAttribute const &_attr, jsonParser &json);

/// \brief Print AtomPosition to json after applying affine transformation
/// cart2frac*cart()+trans
jsonParser &to_json(const AtomPosition &apos, jsonParser &json,
                    Eigen::Ref<const Eigen::Matrix3d> const &cart2frac);

/// \brief Read AtomPosition from json and then apply affine transformation
/// cart2frac*cart()
void from_json(AtomPosition &apos, const jsonParser &json,
               Eigen::Ref<const Eigen::Matrix3d> const &frac2cart,
               HamiltonianModules const &_modules);

template <>
struct jsonConstructor<AtomPosition> {
  /// \brief Read from json [b, i, j, k], using 'unit' for AtomPosition::unit()
  static AtomPosition from_json(const jsonParser &json,
                                Eigen::Matrix3d const &f2c_mat,
                                HamiltonianModules const &_modules);
};

jsonParser &to_json(const Molecule &mol, jsonParser &json,
                    Eigen::Ref<const Eigen::Matrix3d> const &c2f_mat);

void from_json(Molecule &mol, const jsonParser &json,
               Eigen::Ref<const Eigen::Matrix3d> const &f2c_mat,
               HamiltonianModules const &_modules);

template <>
struct jsonConstructor<Molecule> {
  static Molecule from_json(const jsonParser &json,
                            Eigen::Ref<const Eigen::Matrix3d> const &f2c_mat,
                            HamiltonianModules const &_modules);
};

template <>
struct jsonConstructor<Site> {
  static Site from_json(const jsonParser &json, Lattice const &_home,
                        COORD_TYPE coordtype,
                        std::map<std::string, Molecule> const &mol_map,
                        HamiltonianModules const &_modules);
};

jsonParser &to_json(const Site &value, jsonParser &json, COORD_TYPE coordtype);

void from_json(Site &value, const jsonParser &json, Lattice const &_home,
               COORD_TYPE coordtype,
               std::map<std::string, Molecule> const &mol_map,
               HamiltonianModules const &_modules);

BasicStructure read_prim(fs::path filename, double xtal_tol,
                         HamiltonianModules const *_modules = nullptr);

BasicStructure read_prim(jsonParser const &json, double xtal_tol,
                         HamiltonianModules const *_modules = nullptr);

/// \brief Write prim.json to file
void write_prim(const BasicStructure &prim, fs::path filename, COORD_TYPE mode,
                bool include_va = false);

/// \brief Write prim.json as JSON
void write_prim(const BasicStructure &prim, jsonParser &json, COORD_TYPE mode,
                bool include_va = false);

/*
/// \brief Write prim.json as JSON
void write_prim(SimpleStructure const &prim_struc,
                std::set<std::string> const &dofs,
                std::set<std::string> const &attributes,
                jsonParser &json,
                COORD_TYPE mode,
                bool include_va=false);
*/

// --------- SymmetryIO Declarations
// --------------------------------------------------

void write_symop(const SymGroup &grp, Index i, jsonParser &j);

void write_symgroup(const SymGroup &grp, jsonParser &json);

// --------- ChemicalReference IO Declarations
// --------------------------------------------------

ChemicalReference read_chemical_reference(fs::path filename,
                                          const BasicStructure &prim,
                                          double tol);

ChemicalReference read_chemical_reference(const jsonParser &json,
                                          const BasicStructure &prim,
                                          double tol);

void write_chemical_reference(const ChemicalReference &chem_ref,
                              fs::path filename);

void write_chemical_reference(const ChemicalReference &chem_ref,
                              jsonParser &json);

void write_chemical_reference(const ChemicalReference &chem_ref,
                              jsonParser &json);

// --------- CompositionAxes Declarations
// --------------------------------------------------

struct CompositionAxes {
  CompositionAxes() {}

  /// \brief Read CompositionAxes from file
  CompositionAxes(fs::path _filename);

  /// \brief Read CompositionAxes from JSON
  CompositionAxes(const jsonParser &json);

  /// \brief Read CompositionAxes from file
  void read(fs::path _filename);

  /// \brief Read CompositionAxes from JSON
  void read(const jsonParser &json);

  /// \brief Iterate over list of CompositionConverter and insert each one as an
  /// enumerated axes set with a unique numerical name
  template <typename IterType>
  void insert_enumerated(IterType begin, IterType end);

  /// \brief Erase all enumerated axes and clear this->enumerated
  void erase_enumerated();

  /// \brief Write CompositionAxes to file
  void write(fs::path _filename) const;

  /// \brief Write CompositionAxes to JSON
  void write(jsonParser &json) const;

  /// \brief Set this->curr using key
  void select(std::string key);

  /// \brief True if curr_key is set
  bool has_current_axes() const { return !curr_key.empty(); }

  std::map<std::string, CompositionConverter> all_axes;
  std::set<std::string> enumerated;
  std::string curr_key;
  CompositionConverter curr;

  int err_code = 0;
  std::string err_message;
};

jsonParser &to_json(const CompositionConverter &f, jsonParser &json);

/// \brief Deserialize CompositionConverter from JSON
void from_json(CompositionConverter &f, const jsonParser &json);

/// \brief Read standard axes from JSON, and output to std::map<std::string,
/// CompositionConverter>
template <typename OutputIterator>
OutputIterator read_composition_axes(OutputIterator result,
                                     const jsonParser &json);

}  // namespace CASM

#endif
