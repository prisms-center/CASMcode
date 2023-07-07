#ifndef CASM_InfoInterface_impl
#define CASM_InfoInterface_impl

#include <boost/algorithm/string/classification.hpp>
#include <boost/algorithm/string/split.hpp>
#include <string>

#include "casm/casm_io/Log.hh"
#include "casm/casm_io/dataformatter/DataFormatterTools_impl.hh"
#include "casm/casm_io/dataformatter/DataFormatter_impl.hh"

namespace CASM {

template <typename DataObject>
void print_info_desc(DataFormatterDictionary<DataObject> const &dict,
                     std::ostream &sout) {
  auto it_begin = dict.cbegin();
  auto it_end = dict.cend();
  std::string::size_type len(0);
  for (auto it = it_begin; it != it_end; ++it) {
    if (it->type() == DatumFormatterClass::Operator) continue;
    len = max(len, it->name().size());
  }

  int verbosity = Log::standard;
  bool show_clock = false;
  int indent_space = 4;
  Log log(sout, verbosity, show_clock, indent_space);
  log.set_width(60 + indent_space);
  log.increase_indent();
  for (auto it = it_begin; it != it_end; ++it) {
    if (it->type() == DatumFormatterClass::Operator) continue;
    log.indent() << it->name() << std::endl;
    log.increase_indent();

    std::vector<std::string> paragraphs;
    boost::split(paragraphs, it->description(), boost::is_any_of("\n"),
                 boost::token_compress_on);
    for (auto const &paragraph : paragraphs) {
      log.paragraph(paragraph);
    }

    log.decrease_indent();
    log << std::endl;
  }
}

}  // namespace CASM

#endif
