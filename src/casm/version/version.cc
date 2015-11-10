/*
 *  version.cpp
 *
 */

#include "casm/version/version.hh"

using namespace CASM;

const std::string &CASM::version() {
  static const std::string &ver = "python_api";
  return ver;
};
