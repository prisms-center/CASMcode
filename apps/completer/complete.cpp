#include <iostream>
#include <iterator>
#include "casm/completer/Complete.hh"
#include "casm/completer/Handlers.hh"

using namespace CASM;

namespace CASM {
  namespace Completer {
    Engine build_casm_engine() {
      Completer::Engine casm_engine;

      QueryOption dumbquery;
      casm_engine.push_back(Option(dumbquery.tag(), dumbquery.desc()));

      MonteOption dumbmonte;
      casm_engine.push_back(Option(dumbmonte.tag(), dumbmonte.desc()));

      RunOption dumbrun;
      casm_engine.push_back(Option(dumbrun.tag(), dumbrun.desc()));

      BsetOption dumbbset;
      casm_engine.push_back(Option(dumbbset.tag(), dumbbset.desc()));

      CompositionOption dumbcomposition;
      casm_engine.push_back(Option(dumbcomposition.tag(), dumbcomposition.desc()));

      RefOption dumbref;
      casm_engine.push_back(Option(dumbref.tag(), dumbref.desc()));

      FilesOption dumbfiles;
      casm_engine.push_back(Option(dumbfiles.tag(), dumbfiles.desc()));

      FormatOption dumbformat;
      casm_engine.push_back(Option(dumbformat.tag(), dumbformat.desc()));

      ImportOption dumbimport;
      casm_engine.push_back(Option(dumbimport.tag(), dumbimport.desc()));

      InitOption dumbinit;
      casm_engine.push_back(Option(dumbinit.tag(), dumbinit.desc()));

      PerturbOption dumbperturb;
      casm_engine.push_back(Option(dumbperturb.tag(), dumbperturb.desc()));

      SelectOption dumbselect;
      casm_engine.push_back(Option(dumbselect.tag(), dumbselect.desc()));

      SettingsOption dumbsettings;
      casm_engine.push_back(Option(dumbsettings.tag(), dumbsettings.desc()));

      StatusOption dumbstatus;
      casm_engine.push_back(Option(dumbstatus.tag(), dumbstatus.desc()));

      SuperOption dumbsuper;
      casm_engine.push_back(Option(dumbsuper.tag(), dumbsuper.desc()));

      SymOption dumbsym;
      casm_engine.push_back(Option(dumbsym.tag(), dumbsym.desc()));

      UpdateOption dumbupdate;
      casm_engine.push_back(Option(dumbupdate.tag(), dumbupdate.desc()));

      ViewOption dumbview;
      casm_engine.push_back(Option(dumbview.tag(), dumbview.desc()));

      EnumOption dumbenum;
      casm_engine.push_back(Option(dumbenum.tag(), dumbenum.desc()));

      return casm_engine;
    }
  }
}

int main(int argc, char *argv[]) {


  Completer::Engine casm_engine = Completer::build_casm_engine();


  if(argc == 1) {
    std::cout << casm_engine.probe_options();
  }

  if(argc == 2) {
    std::cout << casm_engine.probe_suboptions(argv[1]);
  }

  if(argc == 3) {
    std::cout << casm_engine.probe_arguments(argv[1], argv[2]);
  }

  return 0;
}
