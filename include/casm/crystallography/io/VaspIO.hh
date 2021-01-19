#ifndef CASM_crystallography_io_VaspIO
#define CASM_crystallography_io_VaspIO

#include "casm/crystallography/SimpleStructure.hh"
#include "casm/global/definitions.hh"
#include "casm/global/enum.hh"

namespace CASM {

class Log;

/// \brief Functions and classes related to VASP input/output
namespace VaspIO {

/** \addtogroup VaspIO

    \brief Functions and classes related to VASP input/output

    \ingroup casmIO

*/

/// \brief Print POSCAR with formating options
///
/// Example:
/// \code
/// std::ostream file("POSCAR");
/// Configuration config;
/// PrintPOSCAR printer(config);
/// printer.title("My system");
/// printer.set_cart();
/// printer.sort();
/// printer.print(file);
/// file.close();
/// \endcode
///
/// \ingroup VaspIO
///

class PrintPOSCAR {
 public:
  using AtomName = std::string;

  using SpeciesMode = xtal::SimpleStructure::SpeciesMode;

  /// \brief Construct PrintPOSCAR object
  ///
  /// By default:
  /// - title = ""
  /// - scale = 1.0
  /// - coordinate mode = frac (Direct)
  /// - atom names line is printed
  /// - No selective dynamics
  /// - atom names appended to each coordinate line
  /// - {"Va", "va", "VA"} atoms ignored
  ///
  PrintPOSCAR(xtal::SimpleStructure _struc, std::string _title = "",
              SpeciesMode _mode = SpeciesMode::ATOM);

  /// \brief Set title
  void set_title(std::string title) { m_title = title; }

  /// \brief Set scaling factor
  void set_scale(double s) { m_scale = s; }

  /// \brief Set coordinate mode to Direct (fractional)
  void set_direct() { m_coord_mode = FRAC; }

  /// \brief Set coordinate mode to fractional (Direct)
  void set_frac() { m_coord_mode = FRAC; }

  /// \brief Set coordinate mode to Cartesian
  void set_cart() { m_coord_mode = CART; }

  /// \brief Set coordinate mode
  void set_coord_mode(COORD_TYPE mode) { m_coord_mode = mode; }

  /// \brief Set selective dynamics off
  void set_selective_dynamics_off() { m_sel_dynamics = false; }

  /// \brief Set selective dynamics on
  void set_selective_dynamics_on() { m_sel_dynamics = true; }

  /// \brief Do not print atom names line
  void set_atom_names_off() { m_atom_names = false; }

  /// \brief Print atom names line
  void set_atom_names_on() { m_atom_names = true; }

  /// \brief Do not append atom name to end of each coordinate line
  void set_append_atom_names_off() { m_append_atom_names = false; }

  /// \brief Append atom name to end of each coordinate line
  void set_append_atom_names_on() { m_append_atom_names = false; }

  /// \brief Access set of atom names which should not be printed, such as for
  /// vacancies
  std::set<std::string> &ignore() { return m_ignore; }

  /// \brief const Access set of atom names which should not be printed, such as
  /// for vacancies
  const std::set<std::string> &ignore() const { return m_ignore; }

  /// \brief Default sort is by atom name
  void sort();

  /// \brief Print POSCAR to stream
  void print(std::ostream &sout) const;

  /// \brief Print POSCAR to log (enables indentation)
  void print(Log &sout) const;

 private:
  std::string m_title;
  SpeciesMode m_species_mode;

  xtal::SimpleStructure m_struc;
  std::vector<Index> m_permute;

  double m_scale;
  COORD_TYPE m_coord_mode;
  bool m_atom_names;
  bool m_sel_dynamics;
  bool m_append_atom_names;

  /// \brief List of atom names which should not be printed (primarily for
  /// vacancies)
  std::set<std::string> m_ignore;
};

}  // namespace VaspIO
}  // namespace CASM

#endif
