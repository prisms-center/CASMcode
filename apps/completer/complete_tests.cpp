#include <iostream>
#include <iterator>
#include "casm/completer/Complete.hh"
#include "completer_functions.hh"
#include "casm/casm_io/DataFormatter.hh"
#include "casm/clex/Configuration.hh"

namespace CASM {
  using namespace Completer;
  typedef ArgHandler::ARG_TYPE ARG_TYPE;

  std::string recast(ARG_TYPE atype) {
    switch(atype) {
    case ARG_TYPE::VOID:
      return "VOID";
    case ARG_TYPE::PATH:
      return "PATH";
    case ARG_TYPE::COMMAND:
      return "COMMAND";
    case ARG_TYPE::SCELNAME:
      return "SCELNAME";
    case ARG_TYPE::QUERY:
      return "QUERY";
    case ARG_TYPE::OPERATOR:
      return "OPERATOR";
    default:
      return "FUCKUP";
    }
  }

  /*
     template <typename T>
     std::ostream &operator<< (std::ostream &out, const std::vector<T> &v) {
     if(!v.empty()) {
     out << '[';
     std::copy(v.begin(), v.end(), std::ostream_iterator<T>(out, ", "));
     out << "\b\b]";
     }
     return out;
     }
     */

  void subopt_test() {
    Suboption sub0("--void", "-v", ARG_TYPE::VOID);

    std::cout << sub0.long_tag() << std::endl;
    std::cout << sub0.short_tag() << std::endl;
    std::cout << sub0.matches("--void") << sub0.matches("-v") << sub0.matches("void") << std::endl;
    std::cout << recast(sub0.argument_type()) << std::endl;
  }

  void opt_test() {
    Suboption sub0("--void", "-v", ARG_TYPE::VOID);
    Suboption sub1("--path", "-p", ARG_TYPE::PATH);
    Suboption sub2("--command", "-c", ARG_TYPE::COMMAND);
    Suboption sub3("--scel", "-s", ARG_TYPE::SCELNAME);
    Suboption sub4("--query", "-q", ARG_TYPE::QUERY);
    Suboption sub5("--operator", "-o", ARG_TYPE::OPERATOR);

    std::vector<Suboption> subopts;
    subopts.push_back(sub0);
    subopts.push_back(sub1);
    subopts.push_back(sub2);
    subopts.push_back(sub3);
    subopts.push_back(sub4);
    subopts.push_back(sub5);

    Option opt("fake", subopts);

    std::cout << opt.tag() << std::endl;
    std::cout << opt.probe_suboptions() << std::endl;
    std::cout << recast(opt.probe_argument_type("--query")) << std::endl;
    std::cout << recast(opt.probe_argument_type("-o")) << std::endl;
    std::cout << opt.matches("fake") << opt.matches("real") << std::endl;

    return;
  }

  void opt_test2() {
    int min_vol, max_vol;
    std::vector<std::string> scellname_list, filter_expr;

    po::options_description desc("'casm enum' usage");
    desc.add_options()
    ("help,h", "Write help documentation")
    ("min", po::value<int>(&min_vol), "Min volume")
    ("max", po::value<int>(&max_vol), "Max volume")
    ("filter,f", po::value<std::vector<std::string> >(&filter_expr)->multitoken(), "Filter configuration enumeration so that")
    ("scellname,n", po::value<std::vector<std::string> >(&scellname_list)->multitoken(), "Enumerate configs for given supercells")
    ("all,a", "Enumerate configurations for all supercells")
    ("supercells,s", "Enumerate supercells")
    ("configs,c", "Enumerate configurations");

    Option opt("testopt", desc);


    std::cout << opt.tag() << std::endl;
    std::cout << opt.probe_suboptions() << std::endl;
    std::cout << recast(opt.probe_argument_type("--help")) << std::endl;
    std::cout << recast(opt.probe_argument_type("-f")) << std::endl;
    std::cout << opt.matches("testopt") << opt.matches("real") << std::endl;

    return;
  }

  void po_test() {
    int min_vol, max_vol;
    std::vector<std::string> scellname_list, filter_expr;

    po::options_description desc("'casm enum' usage");
    desc.add_options()
    ("help,h", "Write help documentation")
    ("min", po::value<int>(&min_vol), "Min volume")
    ("max", po::value<int>(&max_vol), "Max volume")
    ("filter,f", po::value<std::vector<std::string> >(&filter_expr)->multitoken(), "Filter configuration enumeration so that")
    ("scellname,n", po::value<std::vector<std::string> >(&scellname_list)->multitoken(), "Enumerate configs for given supercells")
    ("all,a", "Enumerate configurations for all supercells")
    ("supercells,s", "Enumerate supercells")
    ("configs,c", "Enumerate configurations");

    po::option_description testop("filter,x", po::value<int>(&min_vol)->default_value(5)->value_name("lol"), "Filter configuration enumeration so that");

    std::cout << "FORMAT:" << std::endl;
    std::cout << testop.format_parameter() << std::endl;

    std::cout << testop.canonical_display_name(po::command_line_style::allow_dash_for_short) << std::endl;
    std::cout << testop.canonical_display_name(po::command_line_style::allow_long) << std::endl;

    for(auto it = desc.options().begin(); it != desc.options().end(); ++it) {
      std::cout << (*it)->canonical_display_name(po::command_line_style::allow_dash_for_short) << "    ";
      std::cout << (*it)->canonical_display_name(po::command_line_style::allow_long) << std::endl;
    }

    return;
  }

  Option generate_option(std::string postfix, char beginshort) {
    int phony;

    po::options_description desc("phony target " + postfix);

    desc.add_options()
    ((("aaaa" + postfix + ",") + std::string(1, char(beginshort + 1))).c_str(),
     po::value<int>(&phony)->value_name(ArgHandler::supercell()),
     ("aaaa" + postfix + " info").c_str())

    ((("bbbb" + postfix + ",") + std::string(1, char(beginshort + 2))).c_str(),
     po::value<int>(&phony)->default_value(9)->value_name(ArgHandler::query()),
     ("bbbb" + postfix + " info").c_str())

    ((("cccc" + postfix + ",") + std::string(1, char(beginshort + 3))).c_str(),
     po::value<int>(&phony)->value_name(ArgHandler::operation()),
     ("cccc" + postfix + " info").c_str())

    ((("dddd" + postfix + ",") + std::string(1, char(beginshort + 4))).c_str(),
     po::value<int>(&phony),
     ("dddd" + postfix + " info").c_str())

    ((("eeee" + postfix + ",") + std::string(1, char(beginshort + 5))).c_str(),
     po::value<int>(&phony),
     ("eeee" + postfix + " info").c_str());

    return Option("fake" + postfix, desc);
  }

  void engine_test() {
    Option gen0 = generate_option("zz", 'a');
    Option gen1 = generate_option("yy", 'h');

    Engine testengine;
    testengine.push_back(gen0);
    testengine.push_back(gen1);
    testengine.push_back(generate_option("xx", 'l'));

    std::cout << testengine.probe_options() << std::endl;
    std::cout << testengine.probe_suboptions("fakezz") << std::endl;
    std::cout << testengine.probe_suboptions("fakeyy") << std::endl;
    std::cout << testengine.probe_suboptions("fakexx") << std::endl;
    std::cout << testengine.probe_suboptions("yy") << std::endl;
    //std::cout << recast(testengine.probe_argument_type("fakexx", "--aaaaxx")) << std::endl;
    //std::cout << recast(testengine.probe_argument_type("fakexx", "-n")) << std::endl;
    //std::cout << recast(testengine.probe_argument_type("fakexx", "--ccccxx")) << std::endl;
    //std::cout << recast(testengine.probe_argument_type("fakexx", "--ddddxx")) << std::endl;
    return;
  }

  void argtype_test() {
    int min_vol;
    std::string test = "asdf";
    po::option_description testop
    ("filter,x",
     po::value<int>(&min_vol)->default_value(5)->value_name(ArgHandler::path()),
     "Filter configuration enumeration so that");

    ARG_TYPE determined_type = ArgHandler::determine_type(testop);
    std::cout << "DETERMINED: " << recast(determined_type) << std::endl;
    return;
  }

  void property_test() {
    //I think this is s a DataFormatterDictionary<Configuration>
    DataFormatterDictionary<Configuration> dict = make_dictionary<Configuration>();

    for(auto it = dict.begin(); it != dict.cend(); ++it) {
      if(it->type() ==  BaseDatumFormatter<Configuration>::Property)
        std::cout << it->name() << "    ";
    }

    std::cout << std::endl;

    return;
  }

}
using namespace CASM;

int main(int argc, char *argv[]) {

  //subopt_test();

  //std::cout << std::endl;

  //opt_test();

  //std::cout << std::endl;

  //po_test();

  //std::cout << std::endl;

  //opt_test2();

  //std::cout << std::endl;

  //argtype_test();

  //std::cout << std::endl;

  //engine_test();

  //std::cout <<std::endl;

  property_test();


  return 0;
}
