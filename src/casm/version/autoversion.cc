/*
 *  version.cpp
 *
 */

#include "casm/version/version.hh"

using namespace CASM;

#ifndef TXT_VERSION
#define TXT_VERSION "unknown"
#endif 

const std::string &CASM::version() {
  static const std::string &ver = TXT_VERSION;
  return ver;
};
