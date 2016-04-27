/*
 *  version.cpp
 *
 */

#include "casm/version/version.hh"

using namespace CASM;

const std::string &CASM::version() {
  static const std::string &ver = "v0.2.X_files";
  return ver;
};
