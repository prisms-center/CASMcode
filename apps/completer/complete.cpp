#include <iostream>
#include <iterator>
#include "complete.hh"

using namespace Completer;

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

template <typename T>
std::ostream &operator<< (std::ostream &out, const std::vector<T> &v) {
  if(!v.empty()) {
    out << '[';
    std::copy(v.begin(), v.end(), std::ostream_iterator<T>(out, ", "));
    out << "\b\b]";
  }
  return out;
}

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
  ("filter,f", po::value<std::vector<Completer::query_str> >(&filter_expr)->multitoken(), "Filter configuration enumeration so that")
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
  ("filter,f", po::value<std::vector<Completer::query_str> >(&filter_expr)->multitoken(), "Filter configuration enumeration so that")
  ("scellname,n", po::value<std::vector<std::string> >(&scellname_list)->multitoken(), "Enumerate configs for given supercells")
  ("all,a", "Enumerate configurations for all supercells")
  ("supercells,s", "Enumerate supercells")
  ("configs,c", "Enumerate configurations");

  po::option_description testop("filter,x", po::value<int>(&min_vol), "Filter configuration enumeration so that");

  std::cout << testop.canonical_display_name(po::command_line_style::allow_dash_for_short) << std::endl;
  std::cout << testop.canonical_display_name(po::command_line_style::allow_long) << std::endl;

  for(auto it = desc.options().begin(); it != desc.options().end(); ++it) {
    std::cout << (*it)->canonical_display_name(po::command_line_style::allow_dash_for_short) << "    ";
    std::cout << (*it)->canonical_display_name(po::command_line_style::allow_long) << std::endl;
  }

  return;
}


int main() {

  subopt_test();

  std::cout << std::endl;

  opt_test();

  std::cout << std::endl;

  po_test();

  std::cout << std::endl;

  opt_test2();

  return 0;
}
