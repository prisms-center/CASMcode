#ifndef CASM_DiffTransDatabase
#define CASM_DiffTransDatabase

#include <utility>

#include "casm/database/Database.hh"
#include "casm/kinetics/DiffusionTransformation.hh"

namespace {

  namespace DB {

    class DiffTransDatabase : public Database<PrimPeriodicDiffTransOrbit> {

    };

  }
}

#endif
