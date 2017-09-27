#ifndef CASM_Log
#define CASM_Log

#include <iostream>
#include "casm/external/boost.hh"

namespace CASM {

  class Log {

  public:

    static const int none = 0;
    static const int quiet = 5;
    static const int standard = 10;
    static const int verbose = 20;
    static const int debug = 100;

    /// \brief Construct a Log
    ///
    /// \param ostream The stream to print to
    /// \param verbosity The amount to be printed
    ///
    /// For verbosity:
    /// - 0: print nothing
    /// - 10: print all standard output
    /// - 100: print all possible output
    Log(std::ostream &_ostream = std::cout, int _verbosity = standard, bool _show_clock = false);

    template<int _required_verbosity = standard>
    void calculate(const std::string &what) {
      _add<_required_verbosity>("Calculate", what);
    }

    template<int _required_verbosity = standard>
    void construct(const std::string &what) {
      _add<_required_verbosity>("Construct", what);
    }

    template<int _required_verbosity = standard>
    void generate(const std::string &what) {
      _add<_required_verbosity>("Generate", what);
    }

    template<int _required_verbosity = standard>
    void set(const std::string &what) {
      _add<_required_verbosity>("Set", what);
    }

    template<int _required_verbosity = standard>
    void check(const std::string &what) {
      _add<_required_verbosity>("Check", what);
    }

    template<int _required_verbosity = standard>
    void results(const std::string &what) {
      _add<_required_verbosity>("Results", what);
    }

    template<int _required_verbosity = standard>
    void read(const std::string &what) {
      _add<_required_verbosity>("Read", what);
    }

    template<int _required_verbosity = standard>
    void write(const std::string &what) {
      _add<_required_verbosity>("Write", what);
    }

    template<int _required_verbosity = standard>
    void begin(const std::string &what) {
      _add<_required_verbosity>("Begin", what);
    }

    template<int _required_verbosity = standard>
    void end(const std::string &what) {
      _add<_required_verbosity>("End", what);
    }

    template<int _required_verbosity = standard>
    void warning(const std::string &what) {
      _add<_required_verbosity>("Warning", what);
    }

    template<int _required_verbosity = standard>
    void error(const std::string &what) {
      _add<_required_verbosity>("Error", what);
    }

    template<int _required_verbosity = standard>
    void compiling(const std::string &what) {
      _add<_required_verbosity>("Compiling", what);
    }

    template<int _required_verbosity = standard>
    void custom(const std::string &what) {
      static_assert(_required_verbosity >= none && _required_verbosity <= debug, "CASM::Log _required_verbosity must be <= 100");
      m_print = (m_verbosity >= _required_verbosity);
      if(_print()) {
        *m_stream << "-- " << what << " -- ";
        _add_time();
        *m_stream << std::endl;
      }
    }

    template<int _required_verbosity = standard>
    void custom(const std::string &type, const std::string &what) {
      _add<_required_verbosity>(type, what);
    }


    void restart_clock();

    void show_clock();

    void hide_clock();

    double time_s() const;


    void begin_lap();

    double lap_time() const;


    int verbosity() const;

    void set_verbosity(int _verbosity);

    template<int _required_verbosity>
    Log &require() {
      static_assert(_required_verbosity >= none && _required_verbosity <= debug, "CASM::Log _required_verbosity must be <= 100");
      m_print = (m_verbosity >= _required_verbosity);
      return *this;
    }


    void reset(std::ostream &_ostream = std::cout, int _verbosity = standard, bool _show_clock = false);


    template<typename T>
    friend Log &operator<<(Log &log, const T &msg_details);

    friend Log &operator<<(Log &log, std::ostream & (*fptr)(std::ostream &));

    operator std::ostream &();

    explicit operator bool () {
      return m_print;
    }

    /// \brief Read verbosity level from a string
    static std::pair<bool, int> verbosity_level(std::string s);


  private:

    template<int _required_verbosity = standard>
    void _add(const std::string &type, const std::string &what) {
      static_assert(_required_verbosity >= none && _required_verbosity <= debug, "CASM::Log _required_verbosity must be <= 100");
      m_print = (m_verbosity >= _required_verbosity);
      if(_print()) {
        *m_stream << "-- " << type << ": " << what << " -- ";
        _add_time();
        *m_stream << std::endl;
      }
    }

    void _add_time();

    bool _print() const;


    /// If m_verbosity >= required verbosity, then print
    int m_verbosity;

    /// Whether to print
    bool m_print;

    bool m_show_clock;

    boost::chrono::steady_clock::time_point m_start_time;

    boost::chrono::steady_clock::time_point m_lap_start_time;

    std::ostream *m_stream;

  };

  template<typename T>
  Log &operator<<(Log &log, const T &msg_details) {
    if(log._print()) {
      static_cast<std::ostream &>(log) << msg_details;
    }
    return log;
  }

  Log &operator<<(Log &log, std::ostream & (*fptr)(std::ostream &));


  inline Log &default_log() {
    static Log log;
    return log;
  }

  inline Log &default_err_log() {
    static Log log(std::cerr);
    return log;
  }

  inline Log &null_log() {
    static std::ostream nullout(nullptr);
    static Log log(nullout);
    return log;
  }

  class OStringStreamLog : public Log {

  public:

    /// \brief Construct a StringStreamLog
    ///
    /// \param verbosity The amount to be printed
    ///
    /// For verbosity:
    /// - 0: print nothing
    /// - 10: print all standard output
    /// - 100: print all possible output
    OStringStreamLog(int _verbosity = standard, bool _show_clock = false) :
      Log(m_ss, _verbosity, _show_clock) {}

    std::ostringstream &ss() {
      return m_ss;
    };

    const std::ostringstream &ss() const {
      return m_ss;
    };

  private:

    std::ostringstream m_ss;
  };


  class Logging {

  public:

    Logging(Log &log = default_log(), Log &debug_log = default_log(), Log &err_log = default_err_log()) :
      m_log(&log),
      m_debug_log(&debug_log),
      m_err_log(&err_log) {}

    Log &log() const {
      return *m_log;
    }

    Log &debug_log() const {
      return *m_debug_log;
    }

    Log &err_log() const {
      return *m_err_log;
    }

    static Logging null() {
      return Logging(null_log(), null_log(), null_log());
    }

  private:

    Log *m_log;
    Log *m_debug_log;
    Log *m_err_log;

  };

}

#endif
