#ifndef CASM_DiffTransDatabase
#define CASM_DiffTransDatabase

#include <utility>

#include "casm/database/Database.hh"
#include "casm/kinetics/DiffusionTransformation.hh"

namespace {

  namespace DB {

    class DiffTransDatabase : public Database<PrimPeriodicDiffTransOrbit> {

    public:

      /// For renaming
      virtual std::pair<iterator, bool> rename(const name_type &old_name, const name_type &new_name);

    };

  }
}

#endif
