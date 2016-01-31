#include "casm_functions.hh"

#include "casm/clex/PrimClex.hh"
#include "casm/clex/ConfigSelection.hh"
#include "casm/external/gzstream/gzstream.h"

namespace CASM {

  /// \brief If !_primclex, construct new PrimClex stored in uniq_primclex, then 
  ///        return reference to existing or constructed PrimClex
  ///
  /// \param _primclex Pointer to possibly existing PrimClex
  /// \param uniq_primclex Reference to null std::unique_ptr<PrimClex> to manage PrimClex, if it is constructed
  /// \param root fs::path with path to CASM project passed to PrimClex constructor, if it is constructed
  /// \param sout std::ostream to be passed to PrimClex constructor, if it is constructed
  ///
  /// \returns reference to PrimClex (either newly constructed managed by 
  ///          uniq_primclex, or existing pointed at by _primclex)
  ///
  PrimClex& make_primclex_if_not(PrimClex* _primclex, std::unique_ptr<PrimClex>& uniq_primclex, fs::path root, std::ostream& sout) {
    if(!_primclex) {
      sout << "Initialize primclex: " << root << std::endl << std::endl;
      uniq_primclex.reset(new PrimClex(root, sout));
      sout << "  DONE." << std::endl << std::endl;
      return *uniq_primclex;
    }
    return *_primclex;
  }
  
  /// \brief Return a reference to proper std::ostream
  ///
  /// \param output Output mode: False: use 'sout', True: check 'out_path' and 'gzip' to decide
  /// \param sout stream to use if not writing to file
  /// \param fout will be given an open file if writing to file
  /// \param out_path, where to write if 'output': if "STDOUT", use 'sout'; otherwise filename
  /// \param gzip: if true, write to gzip file
  ///
  /// \return reference to stream to use
  ///
  std::ostream& make_ostream_if(bool output, std::ostream& sout, std::unique_ptr<std::ostream>& fout, fs::path out_path, bool gzip) {
  
    if(output) {
      if(out_path.string() == "STDOUT") {
        return sout;
      }
      
      out_path = fs::absolute(out_path);
      
      if(gzip) {
        fout.reset(new gz::ogzstream(out_path.string().c_str()));
        return *fout;
      } 
      
      fout.reset(new fs::ofstream(out_path));
      return *fout;
    }
    else {
      return sout;
    }
    
  }
  
  ConstConfigSelection make_config_selection(std::string selection_str, const PrimClex& primclex) {
    if(selection_str == "MASTER") {
      return ConstConfigSelection(primclex);
    }
    else if(selection_str == "ALL") {
      ConstConfigSelection selection(primclex);
      for(auto it=primclex.config_cbegin(); it!=primclex.config_cend(); ++it) {
        selection.set_selected(it->name(), true);
      }
      return selection;
    }
    else if(selection_str == "CALCULATED") {
      ConstConfigSelection selection(primclex);
      for(auto it=primclex.config_cbegin(); it!=primclex.config_cend(); ++it) {
        selection.set_selected(it->name(), is_calculated(*it));
      }
      return selection;
    }
    else {
      fs::path selection_path = fs::absolute(fs::path(selection_str));
      if(!fs::exists(selection_path)) {
        std::stringstream ss;
        ss << "ERROR in parsing configuation selection name. \n"
           << "  Expected <filename>, 'ALL', 'CALCULATED', or 'MASTER' <--default \n"
           << "  Received: '" << selection_str << "'\n"
           << "  No file named '" << selection_path << "'.";
        throw std::runtime_error(ss.str());
      }
      return ConstConfigSelection(primclex, selection_path);
    }
  }
  


}
