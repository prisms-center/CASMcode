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
  
  
  Log::operator std::ostream&() {
    return *m_stream;
  }
    
  void Log::_add_time() {
    if(m_show_clock) {
      using namespace boost::chrono;
      auto curr_time = steady_clock::now();
      double s = duration_cast<duration<double> >(curr_time - m_start_time).count();
      std::cout << "Time: " << s << " (s)";
    }
  }
  
  bool Log::_print() const {
    return m_print;
  }
  
    
  Log& operator<<(Log& log, std::ostream& (*fptr)(std::ostream&)) {
    if(log._print()) {
      fptr(static_cast<std::ostream&>(log));
    }
    return log;
  }

}