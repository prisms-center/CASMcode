#include "casm/app/DirectoryStructure.hh"
#include "casm/casm_io/Log.hh"

namespace CASM {

  /// \brief Remove files recursively
  void recurs_rm_files(fs::path p, bool dry_run, Log &log) {
    if(!fs::exists(p)) {
      return;
    }

    auto it = fs::directory_iterator(p);
    auto end = fs::directory_iterator();

    for(; it != end; ++it) {
      if(fs::is_regular_file(*it)) {
        log << "rm " << *it << std::endl;
        if(!dry_run) {
          fs::remove(*it);
        }
      }
      else {
        recurs_rm_files(*it, dry_run, log);
        log << "rm " << *it << std::endl;
      }
    }

    log << "rm " << p << std::endl;
    if(!dry_run) {
      fs::remove(p);
    }
  }

}

