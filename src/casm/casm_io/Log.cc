#include "casm/casm_io/Log.hh"

namespace CASM {

  Log::Log(std::ostream &_ostream, int _verbosity, bool _show_clock) {
    reset(_ostream, _verbosity, _show_clock);
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


  void Log::reset(std::ostream &_ostream, int _verbosity, bool _show_clock) {
    m_verbosity = _verbosity;
    m_print = true;
    m_show_clock = _show_clock;
    m_start_time = boost::chrono::steady_clock::now();
    m_stream = &_ostream;
  }


  Log::operator std::ostream &() {
    return *m_stream;
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


  void Log::_add_time() {
    if(m_show_clock) {
      std::cout << "Time: " << time_s() << " (s)";
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