#ifndef CASM_DiffTransDatabase
#define CASM_DiffTransDatabase

#include <utility>

#include "casm/database/Database.hh"
#include "casm/kinetics/DiffusionTransformation.hh"

namespace CASM {

  namespace DB {

    template<>
    class Database<Kinetics::PrimPeriodicDiffTransOrbit> :
      public ValDatabase<Kinetics::PrimPeriodicDiffTransOrbit> {

    };

  }
}

#endif
