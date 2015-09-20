#include "casm/casm_io/FileSystemInterface.hh"

#include <iostream>

#define BOOST_NO_SCOPED_ENUMS
#define BOOST_NO_CXX11_SCOPED_ENUMS
#include <boost/filesystem.hpp>
#include <iterator>
#include <algorithm>
#include <fstream>

#include "casm/CASM_global_definitions.hh"

namespace CASM {

  // Note on the following functions: Some of the functions are adapted/copied from the boost filesystem tutorial on the boost website

  //dirExist Returns if the given path is a directory:
  //      0-> doesn't exist
  //     -1-> Its a file
  //      1-> DIRECTORY
  //     -2-> unknown type, but exists
  // Quits if the read operation fails
  int dirExist(boost::filesystem::path fPath) {
    try {
      if(boost::filesystem::exists(fPath)) {   // does p actually exist?
        if(boost::filesystem::is_regular_file(fPath)) {        // is p a regular file?
          return TYPEFILE;
        }

        else if(boost::filesystem::is_directory(fPath)) {     // is p a directory?
          return TYPEDIR;
        }
        else {
          return TYPEOTHER;
        }
      }
      else {
        return IOERR;
      }

    }

    catch(const boost::filesystem::filesystem_error &ex) {
      std::cerr << ex.what() << '\n';
      exit(666);
    }

  }

  // Writes an error file into path/IO_ERROR. Error filename can be changed here
  // if we decide to call it something else.
  void writeErrorFile(boost::filesystem::path dirPath, std::string errorDir) {
    std::ofstream errorFile;
    dirPath /= "IO_ERROR";
    errorFile.open((dirPath.string()).c_str(), std::ios::out | std::ios::app);
    errorFile << "There was an error trying to write: " << errorDir << std::endl;
    errorFile.close();
  }

  // Returns the path with the last level removed. Better(?) way to do this:
  // append .. to path, and then get canonical.
  boost::filesystem::path popEnd(boost::filesystem::path dirPath) {
    boost::filesystem::path::iterator almostEnd;
    almostEnd = --dirPath.end();
    boost::filesystem::path poppedPath;
    for(boost::filesystem::path::iterator it = dirPath.begin(); it != almostEnd; ++it)
      poppedPath /= *it;
    return poppedPath;
  }

  // Makes a directory and catches for any write errors
  // If directory exists, the IO_ERROR file in the directory
  // is updated. If a file of the same name exists, the IO_ERROR
  // in the parent directory is updated.
  // Quits out if there is an error more serious than that.
  bool makeDirectory(boost::filesystem::path dirPath, std::ostream &stream, bool errorWarn) {
    if(boost::filesystem::create_directory(dirPath)) {
      return true;
    }
    else {
      if(errorWarn) {
        if(dirExist(dirPath) == TYPEDIR) {
          stream << "Directory already exists.\n";
          stream << "Culprit was " << dirPath.string() << std::endl;
          writeErrorFile(dirPath, (--dirPath.end())->filename().string());
        }
        else if(dirExist(dirPath) == TYPEFILE) {
          stream << "There is a file by the same name. Check error file in parent directory \n";
          stream << "Culprit was " << dirPath.string() << std::endl;
          std::string fileName = (--dirPath.end())->filename().string();
          dirPath = popEnd(dirPath);
          writeErrorFile(dirPath, fileName);
        }
        else if(dirExist(dirPath) == TYPEOTHER) {
          stream << "Something weird is happening, it looks like something of the same name already exists but is neither a file nor a directory\nQUITTING!";
          stream << "Culprit was " << dirPath.string() << std::endl;
          exit(666);
        }
        else if(dirExist(dirPath) == IOERR) {
          stream << "Something weird is happening, it looks like there is an IO error\nQUITTING!";
          stream << "Culprit was " << dirPath.string() << std::endl;
          exit(666);
        }
      }
      return false;
    }
  }

  /*
   * Recursively copy a directory using boost::filesystem. This routine is meant
   * to be safe by checking that nothing gets overwritten.
   */

  bool copyDirectory(fs::path const &source, fs::path const &destination) {
    try {
      //Check if function call is valid
      if(!fs::exists(source) || !fs::is_directory(source)) {
        std::cerr << "WARNING in copyDirectory" << std::endl;
        std::cerr << source.string() << " does not exist or is not a directory." << std::endl;
        return false;
      }

      if(fs::exists(destination)) {
        std::cerr << "WARNING in copyDirectory" << std::endl;
        std::cerr << "Destination directory " << destination.string() << " already exists. Will not copy." << std::endl;
        return false;
      }
      //Make directory
      if(!fs::create_directory(destination)) {
        std::cerr << "WARNING in copyDirectory" << std::endl;
        std::cerr << "Unable to create destination directory " << destination.string() << std::endl;
        return false;
      }
    }
    catch(fs::filesystem_error const &e) {
      std::cerr << e.what() << std::endl;
      return false;
    }

    //If you got here all is good and you want to recursively copy every file
    for(fs::directory_iterator file(source); file != fs::directory_iterator(); ++file) {
      try {
        fs::path current(file->path());
        //Hit a directory. Use recursion
        if(fs::is_directory(current)) {
          if(!copyDirectory(current, destination / current.filename())) {
            return false;
          }
        }
        //Hit a file. Copy
        else {
          fs::copy_file(current, destination / current.filename());
        }
      }
      catch(fs::filesystem_error const &e) {
        std::cerr << e.what() << std::endl;
      }
    }
    return true;
  }

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
