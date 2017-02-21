#ifndef CASM_DiffTransDatabase
#define CASM_DiffTransDatabase

#include <utility>

#include "casm/database/Database.hh"
#include "casm/kinetics/DiffusionTransformation.hh"

namespace CASM {

  namespace DB {

    template<>
    class Database<PrimPeriodicDiffTransOrbit> :
      public ValDatabase<PrimPeriodicDiffTransOrbit, std::string> {

    };

  }
}

#endif
