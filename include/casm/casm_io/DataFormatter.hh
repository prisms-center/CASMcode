#ifndef DATAFORMATTER_HH
#define DATAFORMATTER_HH

#include <iostream>
#include <string>
#include <sstream>
#include <vector>
#include <functional>
#include <boost/tokenizer.hpp>
#include "casm/CASM_global_definitions.hh"
#include "casm/misc/CASM_math.hh"
#include "casm/casm_io/jsonParser.hh"
#include "casm/casm_io/DataStream.hh"
#include "casm/casm_io/FormatFlag.hh"


namespace CASM {

  class DataStream;

  template<typename DataObject>
  class BaseDatumFormatter;

  // Given expression string
  //   "subexpr1 subexpr2(subsub1) subexpr3(subsub2(subsubsub1))"
  // splits expresion so that
  //    tag_names = {"subexpr1", "subexpr2", "subexpr3"}
  //    sub_exprs = {"","subsub1", "subsub2(subsubsub1)"}
  // whitespace and commas are ignored
  void split_formatter_expression(const std::string &input_expr,
                                  std::vector<std::string> &tag_names,
                                  std::vector<std::string> &sub_exprs);

  /* A DataFormatter performs extraction of disparate types of data from a 'DataObject' class that contains
   * various types of unasociated 'chunks' of data.
   * The DataFormatter is composed of one or more 'DatumFormatters', with each DatumFormatter knowing how to
   * access and format a particular type of data from the 'DataObject'.
   * BaseDatumFormatter<DataObject> is a virtual class from which all DatumFormatters that access the particular
   * DataObject derive.
   *
   * As an example consider a Configuration object, which has numerous attributes (name, energy, composition, etc).
   * These attributes can be accessed and formatted using a DataFormatter<Configuration>, which contains a number of
   *
   * DatumFormatter<Configuration> objects.  For example
   *      NameConfigDatumFormatter        (derived from DatumFormatter<Configuration>)
   *      EnergyConfigDatumFormatter      (derived from DatumFormatter<Configuration>)
   *      CompositionConfigDatumFormatter (derived from DatumFormatter<Configuration>)
   *
   *  Example usage case:
   *        DataFormatter<DataObject> my_data_formatter;
   *        initialize_data_formatter_from_user_input(my_data_formatter, user_input);
   *        std::cout << my_data_formatter(primclex.config_begin(), primclex.config_end());
   *        my_json_parser = my_data_formatter(primclex.config_begin(), primclex.config_end());
   *
   */
  template<typename DataObject>
  class DataFormatter {
    // These are private classes that sweeten the I/O syntax
    // they have to be forward declared ahead of DataFormatter::operator()
    template<typename IteratorType>
    class FormattedIteratorPair;

    class FormattedObject;
  public:
    DataFormatter(int _sep = 4, int _precision = 12, std::string _comment = "#") :
      m_initialized(false), m_prec(_precision), m_sep(_sep), m_indent(0), m_comment(_comment) {}

    template<typename ...Args>
    DataFormatter(const Args &... formatters)
      : DataFormatter() {
      push_back(formatters...);
    }
    DataFormatter(const DataFormatter<DataObject> &RHS) :
      m_initialized(false), m_col_sep(RHS.m_col_sep), m_col_width(RHS.m_col_width),
      m_prec(RHS.m_prec), m_sep(RHS.m_sep), m_indent(0), m_comment(RHS.m_comment) {

      auto it(RHS.m_data_formatters.cbegin()), it_end(RHS.m_data_formatters.cend());
      for(; it != it_end; ++it)
        m_data_formatters.push_back((*it)->clone());
    }

    ~DataFormatter() {
      clear();
    }

    bool empty() const {
      return m_data_formatters.size() == 0;
    }

    DataFormatter &operator=(const DataFormatter<DataObject> &RHS) {
      if(&RHS == this)
        return *this;

      clear();
      m_indent = RHS.m_indent;
      m_col_sep = RHS.m_col_sep;
      m_col_width = RHS.m_col_width;
      m_prec = RHS.m_prec;
      m_sep = RHS.m_sep;
      m_comment = RHS.m_comment;
      auto it(RHS.m_data_formatters.cbegin()), it_end(RHS.m_data_formatters.cend());
      for(; it != it_end; ++it)
        m_data_formatters.push_back((*it)->clone());
      return *this;
    }

    void clear() {
      for(Index i = 0; i < m_data_formatters.size(); i++)
        delete m_data_formatters[i];
      m_data_formatters.clear();
    }

    template<typename IteratorType>
    FormattedIteratorPair<IteratorType> operator()(IteratorType begin, IteratorType end) const {
      return FormattedIteratorPair<IteratorType>(this, begin, end);
    }

    void set_indent(int _indent) {
      if(m_col_width.size() == 0)
        return;
      m_col_width[0] += _indent;
      m_col_sep[0] += _indent;
      //m_indent=_indent;
    }

    FormattedObject operator()(const DataObject &data_obj) const {
      return FormattedObject(this, data_obj);
    }

    /// Verify that _obj has valid data for all portions of query
    bool validate(const DataObject &_obj) const;

    ///Output selected data from DataObject to DataStream
    void inject(const DataObject &_obj, DataStream &_stream) const;

    ///Output data as specified by *this of the given DataObject
    void print(const DataObject &_obj, std::ostream &_stream) const;

    ///Output data as specified by *this of the given DataObject
    jsonParser &to_json(const DataObject &_obj, jsonParser &json) const;

    ///print the header, using _tmplt_obj to inspect array sizes, etc.
    void print_header(const DataObject &_tmplt_obj, std::ostream &_stream) const;

    /// Add a particular BaseDatumFormatter to *this
    /// If the previous Formatter matches the new formatter, try to just parse the new args into it
    void push_back(const BaseDatumFormatter<DataObject> &new_formatter, const std::string &args) {
      //If the last formatter matches new_formatter, try to parse the new arguments into it
      if(m_data_formatters.size() > 0 && m_data_formatters.back()->name() == new_formatter.name()) {
        if(m_data_formatters.back()->parse_args(args))
          return;
      }

      m_data_formatters.push_back(new_formatter.clone());
      m_col_sep.push_back(0);
      m_col_width.push_back(0);
      m_data_formatters.back()->parse_args(args);

    }

    void push_back(const BaseDatumFormatter<DataObject> &new_formatter) {

      m_data_formatters.push_back(new_formatter.clone());
      m_col_sep.push_back(0);
      m_col_width.push_back(0);
    }

    template<typename ...Args>
    void push_back(const BaseDatumFormatter<DataObject> &new_formatter, const Args &... formatters) {
      push_back(new_formatter);
      push_back(formatters...);
    }

    DataFormatter<DataObject> &operator<<(const BaseDatumFormatter<DataObject> &new_formatter) {
      push_back(new_formatter);
      return *this;
    }

    void set_header_prefix(const std::string &_prefix) {
      m_comment += _prefix;
    }
  private:
    mutable bool m_initialized;
    //List of all the ConfigFormatter objects you want outputted
    std::vector<BaseDatumFormatter<DataObject> *> m_data_formatters;
    mutable std::vector<Index> m_col_sep;
    mutable std::vector<Index> m_col_width;
    //Decimal precision
    int m_prec;
    //number of spaces between columns
    int m_sep;

    int m_indent;
    //comment prefix -- default to "#"
    std::string m_comment;

    void _initialize(const DataObject &_tmplt) const;
  };

  /*
   * BaseDatumFormatter<DataObject> This is an abstract base class from which all other
   * DatumFormatter<DataObject> classes will inherit. The job of a DatumFormatter is to
   * access and format a particular type of data that is stored in a <DataObject> class,
   * which is a template paramter.
   *
   */

  template<typename DataObject>
  class BaseDatumFormatter {
  public:

    enum FormatterType {Property, Operator};
    typedef long difference_type;


    BaseDatumFormatter(const std::string &_init_name, const std::string &_desc) :
      m_name(_init_name), m_description(_desc) {}

    /// Allow polymorphic deletion
    virtual ~BaseDatumFormatter() {};

    /// \brief Returns a name for the formatter, which becomes the tag used for parsing
    const std::string &name() const {
      return m_name;
    }

    /// \brief Returns a short description of the formatter and its allowed arguments (if any).
    /// This description is used to automatically generate help screens
    const std::string &description() const {
      return m_description;
    }

    virtual FormatterType type() const {
      return Property;
    }

    /// \brief Make an exact copy of the formatter (including any initialized members)
    ///
    virtual BaseDatumFormatter *clone() const = 0;
    /**{ return new DerivedDatumFormatter(*this);}**/


    virtual void init(const DataObject &_template_obj) const {

    };

    ///\brief Returns true if _data_obj has valid values for requested data
    virtual bool validate(const DataObject &_data_obj) const {
      return true;
    }

    ///\brief Returns a long expression for each scalar produced by the formatter
    /// parsing the long_header should reproduce the exact query described by the formatter
    /// Ex: "clex(formation_energy)" or "comp(a)    comp(c)"
    virtual std::string long_header(const DataObject &_template_obj) const {
      return short_header(_template_obj);
    };

    ///\brief Returns a short expression for the formatter
    /// parsing the short_header should allow the formatter to be recreated
    /// (but the short header does not specify a subset of the elements)
    /// Ex: "clex(formation_energy)" or "comp"
    virtual std::string short_header(const DataObject &_template_obj) const {
      return name();
    };

    /// If data must be printed on multiple rows, returns number of rows needed to output all data from _data_obj
    /// DataFormatter class will subsequently pass over _data_obj multiple times to complete printing (if necessary)
    virtual Index num_passes(const DataObject &_data_obj) const {
      return 1;
    }

    /// Print formatted data from _data_obj to _stream, while specifying which output pass is requested
    /// If implementation does not depend on pass_index, it may safely be ignored
    virtual void print(const DataObject &_data_obj, std::ostream &_stream, Index pass_index = 0) const = 0;

    /// Stream selected data from _data_obj to _stream, while specifying which output pass is requested
    /// If implementation does not depend on pass_index, it may safely be ignored
    virtual void inject(const DataObject &_data_obj, DataStream &_stream, Index pass_index = 0) const = 0;

    /// Assumes that 'json' object is simply assigned, and it is the job of DataFormatter (or  some other managing entity)
    /// to pass the correct 'json' object.
    ///     Ex:  DerivedDatumFormatter my_formatter;
    ///          initialize(my_formatter); // does some set of initialization steps
    ///          jsonParser my_big_data_object;
    ///          my_formatter.to_json(my_data_object, my_big_data_object["place_to_write"]["my_formatter_data"]);
    virtual jsonParser &to_json(const DataObject &_data_obj, jsonParser &json)const = 0;

    /// If DatumFormatter accepts arguments, parse them here.  Arguments are assumed to be passed from the command line
    /// via:         formattername(argument1,argument2,...)
    ///
    /// from which DerivedDatumFormatter::parse_args() receives the string "argument1,argument2,..."
    /// Returns true if parse is successful, false if not (e.g., takes no arguments, already initialized, malformed input, etc).
    virtual bool parse_args(const std::string &args) {
      return args.size() == 0;
    }
  protected:
    typedef multivector<Index>::X<2> IndexContainer;

    /// Derived DatumFormatters have some optional functionality for parsing index expressions in order to make it easy to handle
    ///       formatter_name(3,4:8)
    /// in which case, DerivedDatumFormatter::parse_args() is called with the string "3,4:8"
    /// by dispatching that string to BaseDatumFormatter::_parse_index_expression(), m_index_rules will be populated with
    /// {{3,4},{3,5},{3,6},{3,7},{3,8}}
    void _parse_index_expression(const std::string &_expr);
    void _add_rule(const std::vector<Index> &new_rule) const {
      m_index_rules.push_back(new_rule);
    }
    const IndexContainer &_index_rules() const {
      return m_index_rules;
    }
  private:
    std::string m_name;
    std::string m_description;
    mutable IndexContainer  m_index_rules;
  };

  template<typename T>
  bool always_true(const T &) {
    return true;
  };

  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  // Couldn't figure out how to get compiler to correctly substitute names
  // virtual "FormattedPrintable" class is a workaround
  class FormattedPrintable {
  public:
    virtual ~FormattedPrintable() {}
    virtual void inject(DataStream &stream)const = 0;
    virtual void print(std::ostream &stream)const = 0;
    virtual jsonParser &to_json(jsonParser &json)const = 0;
  };


  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  template<typename DataObject> template<typename IteratorType>
  class DataFormatter<DataObject>::FormattedIteratorPair : public FormattedPrintable {
    DataFormatter<DataObject> const *m_formatter_ptr;
    IteratorType m_begin_it, m_end_it;
  public:
    FormattedIteratorPair(const DataFormatter<DataObject> *_formatter_ptr, IteratorType _begin, IteratorType _end):
      m_formatter_ptr(_formatter_ptr), m_begin_it(_begin), m_end_it(_end) {}

    void inject(DataStream &_stream) const {
      if(m_begin_it == m_end_it)
        return;
      // hack: always print header to initialize things, like Clexulator, but in this case throw it away
      std::stringstream _ss;
      m_formatter_ptr->print_header(*m_begin_it, _ss);
      for(IteratorType it(m_begin_it); it != m_end_it; ++it)
        m_formatter_ptr->inject(*it, _stream);
    }

    void print(std::ostream &_stream) const {
      if(m_begin_it == m_end_it)
        return;

      FormatFlag format(_stream);
      if(format.print_header()) {
        m_formatter_ptr->print_header(*m_begin_it, _stream);
      }
      else { // hack: always print header to initialize things, like Clexulator, but in this case throw it away
        std::stringstream _ss;
        m_formatter_ptr->print_header(*m_begin_it, _ss);
      }
      format.print_header(false);
      _stream << format;
      for(IteratorType it(m_begin_it); it != m_end_it; ++it)
        m_formatter_ptr->print(*it, _stream);
    }

    jsonParser &to_json(jsonParser &json) const {
      json.put_array();
      for(IteratorType it(m_begin_it); it != m_end_it; ++it)
        json.push_back((*m_formatter_ptr)(*it));
      return json;
    }

  };

  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  template<typename DataObject>
  class DataFormatter<DataObject>::FormattedObject : public FormattedPrintable {
    DataFormatter<DataObject> const *m_formatter_ptr;
    DataObject const *m_obj_ptr;
  public:
    FormattedObject(const DataFormatter<DataObject> *_formatter_ptr, const DataObject &_obj):
      m_formatter_ptr(_formatter_ptr), m_obj_ptr(&_obj) {}

    void inject(DataStream &_stream) const {
      m_formatter_ptr->inject(*m_obj_ptr, _stream);
    }

    void print(std::ostream &_stream) const {
      FormatFlag format(_stream);
      if(format.print_header())
        m_formatter_ptr->print_header(*m_obj_ptr, _stream);
      format.print_header(false);
      _stream << format;

      m_formatter_ptr->print(*m_obj_ptr, _stream);
    }

    jsonParser &to_json(jsonParser &json) const {
      m_formatter_ptr->to_json(*m_obj_ptr, json);
      return json;
    }
  };

  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  /// Parsing dictionary for constructing a DataFormatter<DataObject> object.
  template<typename DataObject>
  class DataFormatterDictionary {
  public:
    DataFormatterDictionary() {};
    DataFormatterDictionary(const DataFormatterDictionary &_dict) {
      typename container::const_iterator it(_dict.m_formatter_map.cbegin()), it_end(_dict.m_formatter_map.cbegin());
      for(; it != it_end; ++it) {
        add_formatter(*(it->second));
      }
    }

    DataFormatterDictionary &operator=(const DataFormatterDictionary &_dict) {
      m_formatter_map.clear();
      typename container::const_iterator it(_dict.m_formatter_map.cbegin()), it_end(_dict.m_formatter_map.cbegin());
      for(; it != it_end; ++it) {
        add_formatter(*(it->second));
      }
    }

    //DataFormatterDictionary(const std::function<void(DataFormatterDictionary<DataObject>&) > &initializer) {
    //initializer(*this);
    //}

    void init(const std::function<void(DataFormatterDictionary<DataObject>&) > &initializer) {
      initializer(*this);
    }

    const BaseDatumFormatter<DataObject>  &lookup(const std::string &_name) const {
      BaseDatumFormatter<DataObject> const *bdf_ptr;
      if(contains(_name, bdf_ptr)) {
        return *bdf_ptr;
      }
      else {
        throw std::runtime_error("CRITICAL ERROR: Invalid format flag \"" + _name + "\" specified.\n"
                                 + "                Did you mean \"" + bdf_ptr->name() + "\"?\n");

      }

    }


    DataFormatterDictionary &add_formatter(const BaseDatumFormatter<DataObject> &new_formatter) {
      if(m_formatter_map.find(new_formatter.name()) != m_formatter_map.end())
        throw std::runtime_error("DataFormatter " + new_formatter.name() + " already exists in parsing dictionary.\nDuplicates are not allowed.\n");
      else
        m_formatter_map[new_formatter.name()].reset(new_formatter.clone());
      return *this;
    }

    bool contains(std::string key, BaseDatumFormatter<DataObject> const *&result_ptr) const {

      typename container::const_iterator it, it_end(m_formatter_map.end());

      it = m_formatter_map.find(key);
      if(it != it_end) {
        result_ptr = (it->second).get();
        return true;
      }
      std::string lkey(key);
      //convert 'key' to lower case and check that
      std::transform(key.begin(), key.end(), lkey.begin(), tolower);
      it = m_formatter_map.find(key);
      if(it != it_end) {
        result_ptr = (it->second).get();
        return true;
      }

      // If no match, try to use demerescau-levenshtein distance to make a helpful suggestion
      it = m_formatter_map.begin();
      int min_dist(-1);
      for(; it != it_end; ++it) {
        int dist = dl_string_dist(key, it->first);
        if(min_dist < 0 || dist < min_dist) {
          min_dist = dist;
          std::cout << "New best: \"" << it->first << "\" aka \"" << it->second->name() << "\"\n";
          result_ptr = (it->second).get();
        }
      }
      return false;
    }

    BaseDatumFormatter<DataObject> const *find(std::string key) const {
      typename container::const_iterator it(m_formatter_map.find(key));
      if(it == m_formatter_map.cend())
        return NULL;
      return (it->second).get();
    }

    void print_help(std::ostream &_stream,
                    typename BaseDatumFormatter<DataObject>::FormatterType ftype,
                    int width, int separation) const;

    DataFormatter<DataObject> parse(const std::string &input)const;
    DataFormatter<DataObject> parse(const std::vector<std::string> &input)const;
    Index size() const {
      return m_formatter_map.size();
    }
  private:
    typedef std::map<std::string, std::unique_ptr<BaseDatumFormatter<DataObject> > > container;
    container m_formatter_map;
    static void _parse(const std::string &input, std::vector<std::string> &format_tags, std::vector<std::string> &format_args);
  };

  //******************************************************************************
  template<typename DataObject>
  struct DataFormatterParser {
  public:
    static int init(const std::function<void(DataFormatterDictionary<DataObject>&) > &initializer) {
      dictionary().init(initializer);
      return dictionary().size();
    }

    static const BaseDatumFormatter<DataObject> &lookup(const std::string &_name) {
      return dictionary().lookup(_name);
    }

    static DataFormatter<DataObject> parse(const std::string &input) {
      return dictionary().parse(input);
    }

    static DataFormatter<DataObject> parse(const std::vector<std::string> &input) {
      return dictionary().parse(input);
    }

    static void add_custom_formatter(const BaseDatumFormatter<DataObject> &new_formatter) {
      dictionary().add_formatter(new_formatter);
    }

    static void load_aliases(const fs::path &alias_path);

    static void print_help(std::ostream &_stream,
                           typename BaseDatumFormatter<DataObject>::FormatterType ftype = BaseDatumFormatter<DataObject>::Property,
                           int width = 60, int separation = 8) {
      dictionary().print_help(_stream, ftype, width, separation);
    }
  private:
    static DataFormatterDictionary<DataObject> &dictionary();

  };

  //******************************************************************************
  template<typename DataObject>
  DataFormatterDictionary<DataObject> &DataFormatterParser<DataObject>::dictionary() {
    static DataFormatterDictionary<DataObject>
    m_dict;
    return m_dict;
  }

  //******************************************************************************

  inline jsonParser &to_json(const FormattedPrintable &_obj, jsonParser &json) {
    return _obj.to_json(json);
  }

  //******************************************************************************
  inline std::ostream &operator<<(std::ostream &_stream, const FormattedPrintable &_formatted) {
    _formatted.print(_stream);
    return _stream;
  }

  //******************************************************************************
  inline DataStream &operator<<(DataStream &_stream, const FormattedPrintable &_formatted) {
    _formatted.inject(_stream);
    return _stream;
  }


}

#include "casm/casm_io/DataFormatter_impl.hh"

#endif
