/*
 *  version.cpp
 *
 */

#include "casm/version/version.hh"

using namespace CASM;

const std::string &CASM::version() {
  <<< <<< < HEAD
  static const std::string &ver = "v0.2.0_merge_test";
  == == == =
    static const std::string & ver = "v0.1.1_merge";
  >>> >>> > v0.1.1_merge
  return ver;
};
