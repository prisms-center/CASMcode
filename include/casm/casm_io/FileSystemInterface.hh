#ifndef FileSystemInterface_HH
#define FileSystemInterface_HH

#include "casm/CASM_global_definitions.hh"

namespace CASM {

  // Note on the following functions: Some of the functions are adapted/copied from the boost filesystem tutorial on the boost website

  //dirExist Returns if the given path is a directory:
  //      0-> doesn't exist
  //     -1-> Its a file
  //      1-> DIRECTORY
  //     -2-> unknown type, but exists
  // Quits if the read operation fails
  int dirExist(boost::filesystem::path fPath);

  // Writes an error file into path/IO_ERROR. Error filename can be changed here
  // if we decide to call it something else.
  void writeErrorFile(boost::filesystem::path dirPath, std::string errorDir);

  // Returns the path with the last level removed. Better(?) way to do this:
  // append .. to path, and then get canonical.
  boost::filesystem::path popEnd(boost::filesystem::path dirPath);

  // Makes a directory and catches for any write errors
  // If directory exists, the IO_ERROR file in the directory
  // is updated. If a file of the same name exists, the IO_ERROR
  // in the parent directory is updated.
  // Quits out if there is an error more serious than that.
  bool makeDirectory(boost::filesystem::path dirPath, std::ostream &stream = std::cerr, bool errorWarn = true);

  /*
   * Recursively copy a directory using boost::filesystem. This routine is meant
   * to be safe by checking that nothing gets overwritten.
   */

  bool copyDirectory(fs::path const &source, fs::path const &destination);

  //Some example uses of the file
  // To compile on UNIX: g++ FileSystemInterface.cc -o fs_debug -lboost_system -lboost_filesystem
  // int main(int argc, char *argv[]){
  //     boost::filesystem::path p;
  //     p=boost::filesystem::current_path();
  //     std::cout<<"End: "<<(*(--p.end()));
  //     for (boost::filesystem::path::iterator it = p.begin(); it != p.end(); ++it)
  //         cout << " " << *it << '\n';
  // #ifdef BOOST_POSIX_API
  //     std::cout<<p.string();
  // #else // BOOST_WINDOWS_API
  //     std::wcout<<p.wstring();
  // #endif
  //     std::cout<<"Trying to create a directory in the above folder...\n";
  //     p/="boost_dir";
  //     makeDirectory(p);
  //     return 0;
  // }

}

#endif

