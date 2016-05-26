#include <iostream>
#include "complete.hh"

using namespace Completer;

void subopt_test() {
  Suboption sub0("--void", "-v", ARG_TYPE::VOID);

  std::cout << sub0.long_tag() << std::endl;
  std::cout << sub0.short_tag() << std::endl;
  std::cout << sub0.matches("--void") << sub0.matches("-v") << sub0.matches("void") << std::endl;
  std::cout << static_cast<std::underlying_type<ARG_TYPE>::type>(sub0.argument_type()) << std::endl;
}


int main() {
  Suboption sub0("--void", "-v", ARG_TYPE::VOID);
  Suboption sub1("--path", "-p", ARG_TYPE::VOID);
  Suboption sub2("--command", "-c", ARG_TYPE::VOID);
  Suboption sub3("--scel", "-s", ARG_TYPE::VOID);
  Suboption sub4("--query", "-q", ARG_TYPE::VOID);
  Suboption sub5("--operator", "-o", ARG_TYPE::VOID);

  std::vector<Suboption> subopts(6);
  subopts.push_back(sub0);
  subopts.push_back(sub1);
  subopts.push_back(sub2);
  subopts.push_back(sub3);
  subopts.push_back(sub4);
  subopts.push_back(sub5);

  subopt_test();

  return 0;
}
