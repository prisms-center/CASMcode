#include "casm/database/DatabaseHandler.hh"
#include "casm/clex/PrimClex.hh"

#include "casm/database/json/jsonDatabase.hh"


namespace CASM {

  namespace DB {

    /// Constructor
    ///
    /// - Uses PrimClex for constructing Supercell, etc.
    /// - Uses PrimClex.settings().db_type() to determine default database type
    DatabaseHandler::DatabaseHandler(const PrimClex &_primclex) :
      m_primclex(&_primclex),
      m_default_db_name(m_primclex->settings().db_name()) {

      jsonDB::insert(m_db);
    }

  }

}
