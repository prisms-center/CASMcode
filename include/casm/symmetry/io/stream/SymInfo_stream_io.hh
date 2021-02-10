#ifndef CASM_symmetry_SymInfo_stream_io
#define CASM_symmetry_SymInfo_stream_io

#include <iomanip>
#include <iostream>
#include <string>

#include "casm/casm_io/enum/stream_io.hh"
#include "casm/global/definitions.hh"
#include "casm/symmetry/SymInfo.hh"

namespace CASM {

class SymGroup;
struct SymInfo;

ENUM_TRAITS(symmetry_type)
ENUM_IO_DECL(symmetry_type)

/// Options for printing SymInfo
struct SymInfoOptions {
  SymInfoOptions(COORD_TYPE _coord_type = FRAC, double _tol = TOL,
                 Index _prec = 7, bool _print_matrix_tau = false)
      : coord_type(_coord_type),
        tol(_tol),
        prec(_prec),
        print_matrix_tau(_print_matrix_tau) {}
  COORD_TYPE coord_type;
  double tol;
  Index prec;
  bool print_matrix_tau;
};

/// \brief Print SymInfo
void print_sym_info(Log &log, const SymInfo &info,
                    SymInfoOptions opt = SymInfoOptions());

/// \brief Print SymInfo to string
std::string to_string(const SymInfo &info,
                      SymInfoOptions opt = SymInfoOptions());

/// \brief Print symmetry symbol to string
std::string to_brief_unicode(const SymInfo &info,
                             SymInfoOptions opt = SymInfoOptions());

/// \brief Print SymInfo to string
std::string description(const SymOp &op, const xtal::Lattice &lat,
                        SymInfoOptions opt = SymInfoOptions());

/// \brief Print SymGroup with matrix / tau
void description(Log &log, const SymGroup &g, const xtal::Lattice &lat,
                 SymInfoOptions opt = SymInfoOptions());

/// \brief Print SymInfo to brief string
std::string brief_description(const SymOp &op, const xtal::Lattice &lat,
                              SymInfoOptions opt = SymInfoOptions());

/// \brief Print SymGroup with brief string
void brief_description(Log &log, const SymGroup &g, const xtal::Lattice &lat,
                       SymInfoOptions opt = SymInfoOptions());

}  // namespace CASM

#endif
