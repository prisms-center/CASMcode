#include "casm/casm_io/Log.hh"
#include <boost/algorithm/string/split.hpp>
#include <boost/algorithm/string/classification.hpp>
#include "casm/external/MersenneTwister/MersenneTwister.h"

namespace CASM {

  void LogParagraph::print(Log &log) const {
    log.paragraph(text);
  }


  void LogVerbatim::print(Log &log) const {
    log.verbatim(text, indent_first_line);
  }


  Log::Log(std::ostream &_ostream, int _verbosity, bool _show_clock, int _indent_space) :
    m_verbosity(_verbosity),
    m_show_clock(_show_clock),
    m_indent_space(_indent_space),
    m_indent_level(0),
    m_indent_spaces(0),
    m_paragraph_width(100),
    m_justification(JustificationType::Left),
    m_ostream(&_ostream) {

    begin_section();
  }

  void Log::restart_clock() {
    m_start_time = boost::chrono::steady_clock::now();
  }

  void Log::show_clock() {
    m_show_clock = true;
  }

  void Log::hide_clock() {
    m_show_clock = false;
  }

  double Log::time_s() const {
    using namespace boost::chrono;
    auto curr_time = steady_clock::now();
    return duration_cast<duration<double> >(curr_time - m_start_time).count();
  }


  void Log::begin_lap() {
    m_lap_start_time = boost::chrono::steady_clock::now();
  }

  double Log::lap_time() const {
    using namespace boost::chrono;
    auto curr_time = steady_clock::now();
    return duration_cast<duration<double> >(curr_time - m_lap_start_time).count();
  }

  int Log::verbosity() const {
    return m_verbosity;
  }

  void Log::set_verbosity(int _verbosity) {
    m_verbosity = _verbosity;
  }


  void Log::reset(std::ostream &_ostream) {
    m_ostream = &_ostream;
  }

  /// \brief Choose c random unique numbers in [0,n)
  std::vector<int> rand_unique(int n, int c, MTRand &mtrand) {
    std::vector<int> index;
    for(int i = 0; i < n; ++i) {
      index.push_back(i);
    }
    using std::swap;
    int choice;
    int s = index.size();
    std::vector<int> res;
    for(int i = 0; i < c; ++i) {
      choice = mtrand.randInt(s - 1);
      res.push_back(index[choice]);
      swap(index[s - 1], index[choice]);
      s--;
    }

    return res;
  }

  void Log::_print_justified_line(std::vector<std::string> &line, int curr_width) {
    // treat too-long case as left justified
    if(justification() == JustificationType::Left || curr_width + line.size() - 1 >= width()) {
      _print_left_justified_line(line, curr_width);
    }
    else if(justification() == JustificationType::Right) {
      _print_right_justified_line(line, curr_width);
    }
    else if(justification() == JustificationType::Center) {
      _print_center_justified_line(line, curr_width);
    }
    else if(justification() == JustificationType::Full) {
      _print_full_justified_line(line, curr_width);
    }
    else {
      throw std::runtime_error("Log print justification error");
    }
  }

  void Log::_print_left_justified_line(std::vector<std::string> &line, int curr_width) {
    indent();
    for(int i = 0; i < line.size(); ++i) {
      if(i != 0) {
        *this << " ";
      }
      *this << line[i];
    }
    *this << std::endl;
  }

  void Log::_print_right_justified_line(std::vector<std::string> &line, int curr_width) {
    indent();
    std::stringstream ss;
    for(int i = 0; i < line.size(); ++i) {
      if(i != 0) {
        ss << " ";
      }
      ss << line[i];
    }
    if(indent_str().size() + ss.str().size() >= width()) {
      *this << std::string(width() - indent_str().size() - ss.str().size(), ' ') << ss.str() << std::endl;
    }
  }

  void Log::_print_center_justified_line(std::vector<std::string> &line, int curr_width) {
    indent();
    std::stringstream ss;
    for(int i = 0; i < line.size(); ++i) {
      if(i != 0) {
        ss << " ";
      }
      ss << line[i];
    }
    std::string str = ss.str();
    std::string before = std::string((width() - str.size()) / 2, ' ');
    std::string after = std::string(width() - before.size() - str.size(), ' ');
    *this << before  << str << after << std::endl;
  }

  void Log::_print_full_justified_line(std::vector<std::string> &line, int curr_width) {
    indent();
    // add ' ' evenly as much as possible
    while(width() - curr_width >= line.size() - 1) {
      for(int i = 0; i < line.size() - 1; ++i) {
        line[i] += ' ';
        curr_width++;
      }
    }
    // add extra uneven ' ' using random number generator to choose locations
    // but seed based on curr_width to give consistent results
    MTRand mtrand(curr_width);
    std::vector<int> index = rand_unique(line.size() - 1, width() - curr_width, mtrand);
    for(int i = 0; i < index.size(); ++i) {
      line[i] += ' ';
    }
    // print words (which now include spaces)
    for(auto &word : line) {
      *this << word;
    }
    *this << std::endl;
  }

  /// \brief Print indented paragraph with wrapping at Log::width()
  Log &Log::paragraph(std::string text) {
    std::vector<std::string> words;
    boost::split(words, text, boost::is_any_of(" "), boost::token_compress_on);

    // 'curr_width' includes indent and words, but not spaces between them
    int curr_width = indent_str().size();
    std::vector<std::string> line;
    for(int i = 0; i < words.size(); ++i) {
      if(line.size() == 0 || curr_width + line.size() + words[i].size() <= width()) {
        line.push_back(words[i]);
        curr_width += words[i].size();
      }
      else {
        // print not-last line
        _print_justified_line(line, curr_width);

        // begin next line
        line.clear();
        line.push_back(words[i]);
        curr_width = indent_str().size() + words[i].size();
      }
    }
    // print last line
    if(justification() == JustificationType::Full) {
      _print_left_justified_line(line, curr_width);
    }
    else {
      _print_justified_line(line, curr_width);
    }

    return *this;
  }

  /// Print verbatim, but with indentation (optional on first line)
  Log &Log::verbatim(std::string text, bool indent_first_line) {
    std::istringstream input;
    input.str(text);
    std::string first_line;
    if(std::getline(input, first_line)) {
      if(indent_first_line) {
        *this << indent_str();
      }
      *this << first_line << std::endl;
      for(std::string line; std::getline(input, line);) {
        indent() << line << std::endl;
      }
    }
    return *this;
  }

  Log::operator std::ostream &() {
    return ostream();
  }

  std::string Log::invalid_verbosity_msg(std::string s) {
    return std::string("Error: Received '") + s +
           "', expected one of 'none', 'quiet', 'standard', 'verbose', 'debug', "
           "or an int in range [0, 100]";
  }

  /// \brief Read verbosity level from a string
  ///
  /// \returns result, a pair of bool,int
  ///          result.first == true if successfully read,
  ///          and result.second is the verbosity level
  ///
  std::pair<bool, int> Log::verbosity_level(std::string s) {

    auto is_int = [](std::string s) {
      int val;
      if(s.empty() || !isdigit(s[0])) {
        return std::make_pair(false, val);
      }
      char *p;
      val = strtol(s.c_str(), &p, 10);
      return std::make_pair(*p == 0 && val >= 0 && val <= 100, val);
    };

    auto res = is_int(s);
    if(res.first) {
      return res;
    }
    else if(s == "none") {
      return std::make_pair(true, 0);
    }
    else if(s == "quiet") {
      return std::make_pair(true, 5);
    }
    else if(s == "standard") {
      return std::make_pair(true, 10);
    }
    else if(s == "verbose") {
      return std::make_pair(true, 20);
    }
    else if(s == "debug") {
      return std::make_pair(true, 100);
    }
    else {
      return std::make_pair(false, 0);
    }

  };

  bool Log::print() const {
    return _print();
  }


  void Log::_add_time() {
    if(m_show_clock) {
      ostream() << "Time: " << time_s() << " (s)";
    }
  }

  bool Log::_print() const {
    return m_print;
  }


  Log &operator<<(Log &log, std::ostream & (*fptr)(std::ostream &)) {
    if(log._print()) {
      fptr(static_cast<std::ostream &>(log));
    }
    return log;
  }

}
