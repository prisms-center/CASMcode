#ifndef CASM_Log
#define CASM_Log

#include <iostream>
#include <vector>
#include <boost/chrono.hpp>
#include "casm/misc/cloneable_ptr.hh"

namespace CASM {

  class Log;
  struct LogText;
  struct LogParagraph;
  struct LogVerbatim;

  enum class JustificationType {
    Left, Right, Center, Full
  };

  struct LogText : public notstd::Cloneable {

    ABSTRACT_CLONEABLE(LogText)

    explicit LogText(std::string _text) : text(_text) {}

    std::string text;

    virtual void print(Log &log) const = 0;
  };

  struct LogParagraph : public LogText {

    CLONEABLE(LogParagraph)

    explicit LogParagraph(std::string _text) : LogText(_text) {}

    virtual void print(Log &log) const override;
  };

  struct LogVerbatim : public LogText {

    CLONEABLE(LogVerbatim)

    explicit LogVerbatim(std::string _text, bool _indent_first_line = true) :
      LogText(_text), indent_first_line(_indent_first_line) {}

    bool indent_first_line;

    virtual void print(Log &log) const override;
  };

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
    Log(std::ostream &_ostream = std::cout, int _verbosity = standard, bool _show_clock = false, int _indent_space = 2);

    // --- List of section types ---
    // Each adds a header, ends the previous section, and begins a new section
    // Each section is associated with a verbosity level which must be met or
    //   exceded in order for the section to be printed to the stream

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
        ostream() << indent_str() << "-- " << what << " -- ";
        _add_time();
        ostream() << std::endl;
      }
    }

    template<int _required_verbosity = standard>
    void custom(const std::string &type, const std::string &what) {
      _add<_required_verbosity>(type, what);
    }

    /// \brief Begin a section, without header
    template<int _required_verbosity = standard>
    void begin_section() {
      m_required_verbosity.push_back(_required_verbosity);
      m_print = (m_verbosity >= m_required_verbosity.back());
    }

    /// \brief Create a subsection
    ///
    /// - Note creates a subsection, but remembers the previous section, so that
    ///   when 'end_section' is called the previous section's verbosity level
    ///   becomes active again.
    ///
    /// Example:
    /// \code
    /// log.begin("Section A");
    /// log << stuff << std::endl;
    /// log << std::endl;
    ///
    /// log.increase_indent();
    /// log.subsection().begin("Section A.1");
    /// log << indent_str() << subsection_stuff << std::endl;
    /// log << std::endl;
    /// log.end_section();
    /// log.decrease_indent();
    ///
    /// log << "continue in section A" << std::endl;
    /// log << std::endl;
    /// \endcode
    ///
    Log &subsection() {
      begin_section<none>();
      return *this;
    }

    /// \brief End a section
    void end_section() {
      if(m_required_verbosity.size()) {
        m_required_verbosity.pop_back();
      }
    }

    // --- Timing ---
    // Enables timing by section

    void restart_clock();

    void show_clock();

    void hide_clock();

    double time_s() const;


    void begin_lap();

    double lap_time() const;


    // --- Verbosity ---

    int verbosity() const;

    void set_verbosity(int _verbosity);

    template<int _required_verbosity>
    Log &require() {
      static_assert(_required_verbosity >= none && _required_verbosity <= debug, "CASM::Log _required_verbosity must be <= 100");
      m_print = (m_verbosity >= _required_verbosity);
      return *this;
    }


    void reset(std::ostream &_ostream = std::cout);


    // --- Paragraph printing ---

    void set_width(int width) {
      m_paragraph_width = width;
    }

    int width() const {
      return m_paragraph_width;
    }

    void set_justification(JustificationType justification) {
      m_justification = justification;
    }

    JustificationType justification() {
      return m_justification;
    }


    Log &paragraph(std::string text);

    Log &verbatim(std::string text, bool indent_first_line = true);


    // --- List printing

    /// \brief Print a list
    template<typename OutputIterator>
    Log &verbatim_list(OutputIterator begin, OutputIterator end, std::string sep = "- ");


    // --- Stream operations ---

    template<typename T>
    friend Log &operator<<(Log &log, const T &msg_details);

    friend Log &operator<<(Log &log, std::ostream & (*fptr)(std::ostream &));

    operator std::ostream &();

    std::ostream &ostream() {
      return *m_ostream;
    }

    bool print() const;

    explicit operator bool () {
      return m_print;
    }


    // --- Indentation ---
    // Indentation is not coupled to sectioning

    int indent_space() const {
      return m_indent_space;
    }

    std::string indent_str() const {
      return std::string(m_indent_space * m_indent_level + m_indent_spaces, ' ');
    }

    void increase_indent() {
      m_indent_level++;
    }

    void decrease_indent() {
      if(m_indent_level > 0) {
        m_indent_level--;
      }
    }

    void increase_indent_spaces(int n) {
      m_indent_spaces += n;
    }

    void decrease_indent_spaces(int n) {
      m_indent_spaces -= n;
    }

    Log &indent() {
      ostream() << indent_str();
      return *this;
    }

    /// Same as verbatim, but uses stringstream to convert to string first
    template<typename T>
    Log &indent(const T &t) {
      std::stringstream ss;
      ss << t;
      return verbatim(ss.str());
    }



    static std::string invalid_verbosity_msg(std::string s);

    /// \brief Read verbosity level from a string
    static std::pair<bool, int> verbosity_level(std::string s);


  private:

    template<int _required_verbosity = standard>
    void _add(const std::string &type, const std::string &what) {
      static_assert(_required_verbosity >= none && _required_verbosity <= debug, "CASM::Log _required_verbosity must be <= 100");
      end_section();
      begin_section<_required_verbosity>();
      if(_print()) {
        ostream() << indent_str() << "-- " << type << ": " << what << " -- ";
        _add_time();
        ostream() << std::endl;
      }
    }

    void _print_justified_line(std::vector<std::string> &line, int curr_width);
    void _print_left_justified_line(std::vector<std::string> &line, int curr_width);
    void _print_right_justified_line(std::vector<std::string> &line, int curr_width);
    void _print_center_justified_line(std::vector<std::string> &line, int curr_width);
    void _print_full_justified_line(std::vector<std::string> &line, int curr_width);


    void _add_time();

    bool _print() const;

    std::vector<int> m_required_verbosity;

    /// If m_verbosity >= required verbosity, then print
    int m_verbosity;

    /// Whether to print
    bool m_print;

    bool m_show_clock;

    /// indent_str = m_indent_space*m_indent_level + m_indent_spaces
    int m_indent_space;
    int m_indent_level;
    int m_indent_spaces;

    // for paragraph writing
    int m_paragraph_width;
    JustificationType m_justification;

    boost::chrono::steady_clock::time_point m_start_time;

    boost::chrono::steady_clock::time_point m_lap_start_time;

    std::ostream *m_ostream;

  };

  /// \brief Print a list
  ///
  /// - Prints each element in vector to a stringstream, then uses Log::verbatim
  ///   to print into the list.
  /// - Indentation is set to the length of the "sep" string
  ///
  /// Example, with initial indent of 2 spaces, and sep="-- ":
  /// \code
  /// A list:
  ///   -- first value
  ///   -- some value
  ///      that prints
  ///      on multiple lines
  ///   -- last value
  /// \endcode
  template<typename OutputIterator>
  Log &Log::verbatim_list(OutputIterator begin, OutputIterator end, std::string sep) {

    int n_indent_spaces = sep.size();

    bool indent_first_line = false;
    for(auto it = begin; it != end; ++it) {
      indent() << sep;
      std::stringstream ss;
      ss << *it;
      increase_indent_spaces(n_indent_spaces);
      verbatim(ss.str(), indent_first_line);
      decrease_indent_spaces(n_indent_spaces);
    }
    return *this;
  }

  template<typename T>
  Log &operator<<(Log &log, const T &msg_details) {
    if(log._print()) {
      static_cast<std::ostream &>(log) << msg_details;
    }
    return log;
  }

  Log &operator<<(Log &log, std::ostream & (*fptr)(std::ostream &));

  /// A Log whose underlying ostream* cannot be reset
  class FixedLog : public Log {
  public:
    explicit FixedLog(std::ostream &_ostream);
    FixedLog(FixedLog const &) = delete;
    FixedLog &operator=(FixedLog const &RHS) = delete;

  private:
    using Log::reset;

  };

  inline Log &default_log() {
    static Log log {std::cout};
    return log;
  }

  inline Log &default_err_log() {
    static Log log {std::cerr};
    return log;
  }

  inline Log &log() {
    return CASM::default_log();
  }

  inline Log &err_log() {
    return CASM::default_err_log();
  }

  inline Log &cout_log() {
    static FixedLog log {std::cout};
    return log;
  }

  inline Log &cerr_log() {
    static FixedLog log {std::cerr};
    return log;
  }

  inline Log &null_log() {
    static std::ostream nullout {nullptr};
    static FixedLog log {nullout};
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
    OStringStreamLog(int _verbosity = standard, bool _show_clock = false):
      Log(std::cout, _verbosity, _show_clock) {
      reset(m_ss);
    }

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
      m_err_log(&err_log) {
    }

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
