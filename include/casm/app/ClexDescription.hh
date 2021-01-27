#ifndef CASM_ClexDescription
#define CASM_ClexDescription

#include <iosfwd>
#include <string>

namespace CASM {

class jsonParser;
class DirectoryStructure;

/// \brief Specifies a particular cluster expansion
///
/// - Cluster expansions are given an identifying name as a shortcut
/// - Cluster expansions are fully specified via:
///   - property: the property being expanded, for instance 'formation_energy'
///   - calctype: the type of calculations of the property the cluster
///     expansion is fit to
///   - ref: indicates a reference used to calculate the property being
///     expanded
///   - bset: the basis set used
///
struct ClexDescription {
  ClexDescription() {}

  ClexDescription(std::string _name, std::string _property,
                  std::string _calctype, std::string _ref, std::string _bset,
                  std::string _eci)
      : name(_name),
        property(_property),
        calctype(_calctype),
        ref(_ref),
        bset(_bset),
        eci(_eci) {}

  void print(std::ostream &sout, bool is_default, int indent = 0) const;

  std::string name;
  std::string property;
  std::string calctype;
  std::string ref;
  std::string bset;
  std::string eci;
};

/// Create default configuration ClexDescription
///
/// Using:
/// - name,property="formation_energy"
/// - calctype,ref,bset,eci="default"
ClexDescription default_configuration_clex();

/// \brief Compare using name strings: A.name < B.name
bool operator<(const ClexDescription &A, const ClexDescription &B);

jsonParser &to_json(const ClexDescription &desc, jsonParser &json);

void from_json(ClexDescription &desc, const jsonParser &json);

bool clex_exists(const DirectoryStructure &dir, const ClexDescription &desc);

/// Add directories for a particular cluster expansion
///
/// Create bset, calc_settings, ref, and eci directories for the cluster
/// expansion `desc`.
///
/// \returns bool, true if all directories exist afterwards
bool new_dir(const DirectoryStructure &dir, ClexDescription const &desc);

}  // namespace CASM

#endif
