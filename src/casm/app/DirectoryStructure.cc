#include "casm/app/DirectoryStructure.hh"
#include "casm/clex/Configuration.hh"
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

  void recurs_cp_files(const fs::path &from_dir, const fs::path &to_dir, bool dry_run, Log &log) {
    auto it = fs::directory_iterator(from_dir);
    auto end = fs::directory_iterator();
    for(; it != end; ++it) {
      if(fs::is_regular_file(*it)) {
        log << "cp " << *it << " " << to_dir << std::endl;
        if(!dry_run) {
          fs::copy_file(*it, to_dir / it->path().filename());
        }
      }
      else {
        fs::path new_to_dir = to_dir / it->path().filename();
        if(!dry_run) {
          fs::create_directories(new_to_dir);
        }
        recurs_cp_files(*it, new_to_dir, dry_run, log);
      }
    }
  }

}

