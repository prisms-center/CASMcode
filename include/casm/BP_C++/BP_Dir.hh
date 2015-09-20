#ifndef BP_Dir_HH
#define BP_Dir_HH

#include <iostream>
#include <istream>
#include <fstream>
#include <string>
#include <cstdio>
#include <sstream>
#include "casm/BP_C++/BP_Parse.hh"
#include "casm/BP_C++/BP_basic.hh"
#include <errno.h>
#include <dirent.h>
#include <sys/stat.h>
#include <unistd.h>

namespace BP {

  /// \brief BP_Dir, a class for managing directories and files
  /// - list, create, copy, rename (move), delete directories
  /// - list, copy, rename (move), delete files
  /// - use absolute or relative filenames
  ///
  class BP_Dir {
  private:
    std::string _name;
    std::string _abs_name;
    std::string _joint;

    BP_Vec<std::string> files_list;
    BP_Vec<std::string> dirs_list;

  public:

    // // Initialize with current directory:
    // BP_Dir dir;
    //
    // // Get the current directory (just directory name):
    // cout << "The current directory: " << dir.name() << endl;
    //
    // // Get the current directory (include full path):
    // cout << "The current directory: " << dir.abs_name() << endl;
    //
    //
    // // get a list of files or directories in current directory:
    // BP_Vec<std::string> files = dir.files();
    // BP_Vec<std::string> dirs = dir.dirs();
    //
    // // or just use directly:
    // cout << "There are " << dir.files().size() << " files in this directory" << endl;
    //
    // // Change the current directory (of the BP_Dir object):
    // dir.cd( dir.dirs()[0]);
    //
    // // or directly:
    // dir.cd( "the/path/other_directory");
    //
    // // to be more general you can do this:
    // dir.cd( "the" + dir.joint() + "path" + dir.joint() + "other_directory");
    //
    // // Write a file in the new current directory using BP_Write:
    // BP_Write outfile(dir.abs_name("new_file.txt"));
    // outfile.newfile();
    // outfile << "This is a new file" << endl;



    BP_Dir();

    std::string					name() const {
      return _name;
    };
    std::string					abs_name() const {
      return _abs_name;
    };

    void					refresh();
    bool					cd(std::string s);
    bool					mkdir(std::string s);
    bool					rmdir(std::string s);
    bool					rmdir(BP_Vec<std::string> s);
    bool					recurs_rmdir(std::string s);
    bool					recurs_rmdir(BP_Vec<std::string> s);
    bool					rm(std::string s);
    bool					rm(BP_Vec<std::string> s);
    bool					rename(std::string s1, std::string s2);
    bool					copy(std::string s1, std::string s2);
    bool					copydir(std::string s1, std::string s2);
    bool					cat(std::string s1, std::string s2, std::string s3);

    bool					print(std::string s, std::ostream &sout) const;
    const BP_Vec<std::string>	&files() const {
      return files_list;
    };
    const BP_Vec<std::string>	&dirs() const {
      return dirs_list;
    };
    bool					is_file(std::string s) const;
    bool					is_dir(std::string s) const;
    bool					exists(std::string s) const;

    void					set_joint(std::string s) {
      _joint = s;
    };
    std::string					joint() const {
      return _joint;
    };
    std::string					without_path(std::string s) const;

    bool					is_absolute_path(std::string s) const;
    void					make_absolute_path(std::string &s) const;
    std::string					abs_name(std::string s) const;


  private:
    bool generate_files_dirs();
    bool remove_file(std::string s);
    bool remove_dir(std::string s);
    bool recurs_remove_dir(std::string s);
    bool change_directory(std::string s);

  };


}

#endif // BP_Dir_HH

