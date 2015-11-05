#ifndef DATAFORMATTERTOOLS_HH
#define DATAFORMATTERTOOLS_HH
#include <numeric>
#include <iterator>
#include <regex>
#include <boost/algorithm/string.hpp>
#include "casm/casm_io/DataFormatter.hh"
#include "casm/misc/CASM_TMP.hh"

namespace CASM {

  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


  template<typename ValueType, typename ArgType, typename DataObject>
  class DataFormatterOperator : public BaseDatumFormatter<DataObject> {
  public:

    using BaseDatumFormatter<DataObject>::name;
    using Evaluator = std::function<ValueType(const std::vector<ArgType> &)>;
    using Parser = std::function<const BaseDatumFormatter<DataObject>&(const std::string &)>;
    //validator probably not necessary?
    //using Validator = std::function<bool(const DataObject &)>;

    DataFormatterOperator(const std::string &_init_name, const std::string &_desc, Evaluator evaluator, Parser parser = DataFormatterParser<DataObject>::lookup) :
      BaseDatumFormatter<DataObject>(_init_name, _desc), m_evaluator(evaluator), m_parser(parser) {}

    virtual BaseDatumFormatter<DataObject> *clone() const override {
      return new DataFormatterOperator(*this);
    }

    typename BaseDatumFormatter<DataObject>::FormatterType type() const override {
      return BaseDatumFormatter<DataObject>::Operator;
    }

    std::string short_header(const DataObject &_template_obj) const override;

    bool validate(const DataObject &_data_obj) const override {
      return m_arg_formatter.validate(_data_obj);
    }

    bool parse_args(const std::string &_args) override;

    void inject(const DataObject &_data_obj, DataStream &_stream, Index pass_index) const override {
      _stream << _evaluate(_data_obj);
      if(!validate(_data_obj)) {
        _stream << DataStream::failbit;
      }
    }

    void print(const DataObject &_data_obj, std::ostream &_stream, Index pass_index) const override {
      _stream.flags(std::ios::showpoint | std::ios::fixed | std::ios::right);
      _stream.precision(8);
      if(validate(_data_obj))
        _stream << _evaluate(_data_obj);
      else
        _stream << "unknown";
    }

    jsonParser &to_json(const DataObject &_data_obj, jsonParser &json)const override {
      if(validate(_data_obj))
        json = _evaluate(_data_obj);
      return json;
    }

  protected:
    ValueType _evaluate(const DataObject &_data_obj) const {
      VectorDataStream<ArgType> vec_stream;
      vec_stream << m_arg_formatter(_data_obj);
      return m_evaluator(vec_stream.vector());
    }
  private:
    Evaluator m_evaluator;
    Parser m_parser;
    DataFormatter<DataObject> m_arg_formatter;
  };

  template<typename DataObject>
  DataFormatterOperator<double, double, DataObject> format_operator_add() {
    return DataFormatterOperator<double, double, DataObject>("add", "Add two or more numbers",
    [](const std::vector<double> &vec)->double{
      return std::accumulate(vec.cbegin(),
      vec.cend(),
      0.0);
    });

  }

  template<typename DataObject>
  DataFormatterOperator<double, double, DataObject> format_operator_sub() {
    return DataFormatterOperator<double, double, DataObject>("sub", "Subtract two numbers",
    [](const std::vector<double> &vec)->double{
      if(vec.size() != 2)
        throw std::runtime_error("Subtraction operator must receive exactly two values!");
      return vec[0] - vec[1];
    });

  }

  template<typename DataObject>
  DataFormatterOperator<double, double, DataObject> format_operator_mult() {
    return DataFormatterOperator<double, double, DataObject>("mult", "Add two or more numbers",
    [](const std::vector<double> &vec)->double{
      return std::accumulate(vec.cbegin(),
      vec.cend(),
      1.0,
      [](double a, double b)->double{
        return a *b;
      });
    });

  }

  template<typename DataObject>
  DataFormatterOperator<double, double, DataObject> format_operator_div() {
    return DataFormatterOperator<double, double, DataObject>("div", "Divide two numbers",
    [](const std::vector<double> &vec)->double{
      if(vec.size() != 2)
        throw std::runtime_error("Division operator must receive exactly two values!");
      return vec[0] - vec[1];
    });

  }

  template<typename DataObject>
  DataFormatterOperator<double, double, DataObject> format_operator_max() {
    return DataFormatterOperator<double, double, DataObject>("max", "Max value of two or more numbers",
    [](const std::vector<double> &vec)->double{
      return (*std::max_element(vec.cbegin(),
      vec.cend()));
    });
  }

  template<typename DataObject>
  DataFormatterOperator<double, double, DataObject> format_operator_min() {
    return DataFormatterOperator<double, double, DataObject>("min", "Min value of two or more numbers",
    [](const std::vector<double> &vec)->double{
      return (*std::min_element(vec.cbegin(),
      vec.cend()));
    });
  }

  template<typename DataObject>
  DataFormatterOperator<long, double, DataObject> format_operator_imax() {
    return DataFormatterOperator<long, double, DataObject>("imax", "Index (from 0) of max value in array of two or more numbers",
    [](const std::vector<double> &vec)->long{
      auto it = std::max_element(vec.cbegin(), vec.cend());
      return std::distance(vec.cbegin(), it);
    });
  }

  template<typename DataObject>
  DataFormatterOperator<long, double, DataObject> format_operator_imin() {
    return DataFormatterOperator<long, double, DataObject>("imin", "Index (from 0) of min value in array of two or more numbers",
    [](const std::vector<double> &vec)->long{
      auto it = std::min_element(vec.cbegin(), vec.cend());
      return std::distance(vec.cbegin(), it);
    });
  }

  template<typename DataObject>
  DataFormatterOperator<double, double, DataObject> format_operator_exp() {
    return DataFormatterOperator<double, double, DataObject>("exp", "Exponential function",
    [](const std::vector<double> &vec)->double{
      if(vec.size() != 1)
        throw std::runtime_error("Exponent operator must receive exactly one value!");
      return exp(vec[0]);
    });

  }


  template<typename DataObject>
  DataFormatterOperator<bool, std::string, DataObject> format_operator_regex_match() {
    return DataFormatterOperator<bool, std::string, DataObject>("re", "Check if string matches regular expression. Ex: re('input_string','regex_pattern')",
    [](const std::vector<std::string> &vec)->bool{
      if(vec.size() != 2)
        throw std::runtime_error("Operator re('input_string','regex_pattern') must receive exactly 2 values!");
      std::regex e(boost::trim_copy_if(vec[1], boost::is_any_of(" '")));
      return std::regex_match(vec[0], e);
    });
  }

  template<typename DataObject>
  DataFormatterOperator<bool, std::string, DataObject> format_operator_regex_search() {
    return DataFormatterOperator<bool, std::string, DataObject>("rs", "Check if string contains regular expression. Ex: rs('input_string','regex_pattern')",
    [](const std::vector<std::string> &vec)->bool{
      if(vec.size() != 2)
        throw std::runtime_error("Operator re('input_string','regex_pattern') must receive exactly 2 values!");
      std::regex e(boost::trim_copy_if(vec[1], boost::is_any_of(" '")));
      return std::regex_search(vec[0], e);
    });
  }

  template<typename DataObject>
  DataFormatterOperator<double, double, DataObject> format_operator_sq() {
    return DataFormatterOperator<double, double, DataObject>("sq", "Square of a number",
    [](const std::vector<double> &vec)->double{
      if(vec.size() != 1)
        throw std::runtime_error("Square operator must receive exactly one value!");
      return vec[0] * vec[0];
    });

  }

  template<typename DataObject>
  DataFormatterOperator<double, double, DataObject> format_operator_sqrt() {
    return DataFormatterOperator<double, double, DataObject>("sqrt", "Square root of a number",
    [](const std::vector<double> &vec)->double{
      if(vec.size() != 1)
        throw std::runtime_error("Square-root operator must receive exactly one value!");
      return sqrt(vec[0]);
    });

  }

  template<typename DataObject>
  DataFormatterOperator<double, double, DataObject> format_operator_neg() {
    return DataFormatterOperator<double, double, DataObject>("neg", "Negative of a number",
    [](const std::vector<double> &vec)->double{
      if(vec.size() != 1)
        throw std::runtime_error("Negation operator must receive exactly one value!");
      return -vec[0];
    });

  }

  template<typename DataObject>
  DataFormatterOperator<bool, bool, DataObject> format_operator_and() {
    return DataFormatterOperator<bool, bool, DataObject>("and", "Boolean AND for sequence of boolean values",
    [](const std::vector<bool> &vec)->bool{
      return std::accumulate(vec.cbegin(),
      vec.cend(),
      true,
      [](bool a, bool b)->bool{
        return a && b;
      });
    });

  }

  template<typename DataObject>
  DataFormatterOperator<bool, bool, DataObject> format_operator_or() {
    return DataFormatterOperator<bool, bool, DataObject>("or", "Boolean OR for sequence of boolean values",
    [](const std::vector<bool> &vec)->bool{
      return std::accumulate(vec.cbegin(),
      vec.cend(),
      false,
      [](bool a, bool b)->bool{
        return a || b;
      });
    });

  }

  template<typename DataObject>
  DataFormatterOperator<bool, bool, DataObject> format_operator_xor() {
    return DataFormatterOperator<bool, bool, DataObject>("xor", "Boolean XOR for for two boolean values",
    [](const std::vector<bool> &vec)->bool{
      if(vec.size() != 2)
        throw std::runtime_error("Boolean XOR operator expects exactly two values!");
      return (vec[0] && !vec[1]) || (!vec[0] && vec[1]);});
  }

  template<typename DataObject>
  DataFormatterOperator<bool, bool, DataObject> format_operator_not() {
    return DataFormatterOperator<bool, bool, DataObject>("not", "Boolean NOT for a single boolean value",
    [](const std::vector<bool> &vec)->bool{
      if(vec.size() != 1)
        throw std::runtime_error("Boolean NOT operator must receive exactly one value!");
      return !vec[0];
    });

  }

  template<typename DataObject>
  DataFormatterOperator<bool, double, DataObject> format_operator_eq() {
    return DataFormatterOperator<bool, double, DataObject>("eq", "Equality comparison for two values",
    [](const std::vector<double> &vec)->bool{
      if(vec.size() != 2)
        throw std::runtime_error("Greater-than operator must receive exactly two values!");
      return almost_equal(vec[0], vec[1]);
    });

  }

  template<typename DataObject>
  DataFormatterOperator<bool, double, DataObject> format_operator_lt() {
    return DataFormatterOperator<bool, double, DataObject>("lt", "Less-than comparison for two values",
    [](const std::vector<double> &vec)->bool{
      if(vec.size() != 2)
        throw std::runtime_error("Less-than operator must receive exactly two values!");
      return vec[0] < vec[1];
    });

  }

  template<typename DataObject>
  DataFormatterOperator<bool, double, DataObject> format_operator_le() {
    return DataFormatterOperator<bool, double, DataObject>("le", "Less-than-or-equal comparison for two values",
    [](const std::vector<double> &vec)->bool{
      if(vec.size() != 2)
        throw std::runtime_error("Less-than-or-equal operator must receive exactly two values!");
      return vec[0] <= vec[1];
    });

  }

  template<typename DataObject>
  DataFormatterOperator<bool, double, DataObject> format_operator_gt() {
    return DataFormatterOperator<bool, double, DataObject>("gt", "Greater-than comparison for two values",
    [](const std::vector<double> &vec)->bool{
      if(vec.size() != 2)
        throw std::runtime_error("Greater-than operator must receive exactly two values!");
      return vec[0] > vec[1];
    });

  }

  template<typename DataObject>
  DataFormatterOperator<bool, double, DataObject> format_operator_ge() {
    return DataFormatterOperator<bool, double, DataObject>("ge", "Greater-than-or-equal comparison for two values",
    [](const std::vector<double> &vec)->bool{
      if(vec.size() != 2)
        throw std::runtime_error("Greater-than-or-equal operator must receive exactly two values!");
      return vec[0] >= vec[1];
    });

  }


  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  template<typename DataObject>
  class DatumFormatterAlias: public BaseDatumFormatter<DataObject> {
  public:
    using BaseDatumFormatter<DataObject>::name;

    DatumFormatterAlias(const std::string &_name, const std::string &_command) :
      BaseDatumFormatter<DataObject> (_name, "User-specified alias for '" + _command + "'") {

      split_formatter_expression(_command, m_format_tags, m_subexprs);
      if(m_format_tags.size() != 1)
        throw std::runtime_error("Expression '" + _command + "' is either empty or consists of multiple expressions.\n");

    }

    DatumFormatterAlias(const std::string &_name, const BaseDatumFormatter<DataObject> &_rhs) :
      BaseDatumFormatter<DataObject> (_name, ""), m_formatter(_rhs.clone()) {    }

    DatumFormatterAlias(const DatumFormatterAlias &_rhs) :
      BaseDatumFormatter<DataObject> (_rhs), m_format_tags(_rhs.m_format_tags), m_subexprs(_rhs.m_subexprs) {
      if(_rhs.m_formatter)
        m_formatter.reset(_rhs.m_formatter->clone());
    }

    DatumFormatterAlias &operator=(const DatumFormatterAlias &_rhs) {

      BaseDatumFormatter<DataObject> ::operator=(_rhs);
      if(_rhs.m_formatter)
        m_formatter.reset(_rhs.m_formatter->clone());
      return *this;

    }

    typename BaseDatumFormatter<DataObject>::FormatterType type() const {
      return BaseDatumFormatter<DataObject>::Property;
    }

    /// \brief Make an exact copy of the formatter (including any initialized members)
    ///
    BaseDatumFormatter<DataObject> *clone() const override {
      return new DatumFormatterAlias(*this);
    }


    void init(const DataObject &_template_obj) const  override {
      m_formatter->init(_template_obj);
    }

    ///\brief Returns true if _data_obj has valid values for requested data
    bool validate(const DataObject &_data_obj) const  override {
      return m_formatter->validate(_data_obj);
    }

    ///\brief Returns a long expression for each scalar produced by the formatter
    /// parsing the long_header should reproduce the exact query described by the formatter
    /// Ex: "clex(formation_energy)" or "comp(a)    comp(c)"
    std::string long_header(const DataObject &_template_obj) const  override {
      CountDataStream tcount;
      m_formatter->inject(_template_obj, tcount);
      if(tcount.count() == 1)
        return name();

      std::stringstream t_ss;
      for(Index i = 0; i < tcount.count(); i++)
        t_ss << "       " << name() << '(' << i << ')';

      return t_ss.str();
    }

    ///\brief Returns a short expression for the formatter
    /// parsing the short_header should allow the formatter to be recreated
    /// (but the short header does not specify a subset of the elements)
    /// Ex: "clex(formation_energy)" or "comp"
    std::string short_header(const DataObject &_template_obj) const  override {
      return name();
    }

    /// If data must be printed on multiple rows, returns number of rows needed to output all data from _data_obj
    /// DataFormatter class will subsequently pass over _data_obj multiple times to complete printing (if necessary)
    Index num_passes(const DataObject &_data_obj) const  override {
      return m_formatter->num_passes(_data_obj);
    }

    /// Print formatted data from _data_obj to _stream, while specifying which output pass is requested
    /// If implementation does not depend on pass_index, it may safely be ignored
    void print(const DataObject &_data_obj, std::ostream &_stream, Index pass_index = 0) const override {
      m_formatter->print(_data_obj, _stream, pass_index);
    }

    /// Stream selected data from _data_obj to _stream, while specifying which output pass is requested
    /// If implementation does not depend on pass_index, it may safely be ignored
    void inject(const DataObject &_data_obj, DataStream &_stream, Index pass_index = 0) const override {
      m_formatter->inject(_data_obj, _stream, pass_index);
    }


    /// Assumes that 'json' object is simply assigned, and it is the job of DataFormatter (or  some other managing entity)
    /// to pass the correct 'json' object.
    ///     Ex:  DerivedDatumFormatter my_formatter;
    ///          initialize(my_formatter); // does some set of initialization steps
    ///          jsonParser my_big_data_object;
    ///          my_formatter.to_json(my_data_object, my_big_data_object["place_to_write"]["my_formatter_data"]);
    jsonParser &to_json(const DataObject &_data_obj, jsonParser &json)const override {
      return m_formatter->to_json(_data_obj, json);
    }
    /// If DatumFormatter accepts arguments, parse them here.  Arguments are assumed to be passed from the command line
    /// via:         formattername(argument1,argument2,...)
    ///
    /// from which DerivedDatumFormatter::parse_args() receives the string "argument1,argument2,..."
    /// Returns true if parse is successful, false if not (e.g., takes no arguments, already initialized, malformed input, etc).
    bool parse_args(const std::string &args)  override {
      m_formatter.reset(DataFormatterParser<DataObject>::lookup(m_format_tags[0]).clone());
      //Parse the arguments of of the command now.  Later, we may want to do expression substitution
      // (e.g.,  "comp_plus = add(comp(a),$1)" where $1 specifies an argument)
      if(!m_formatter->parse_args(m_subexprs[0])) {
        throw std::runtime_error("Invalid arguments passed to '" + m_format_tags[0] + "'. Cannot accept expression '" + m_subexprs[0] + "'\n");
      }

      return args.size() == 0;
      //return m_formatter->parse_args(args);
    }
  private:
    std::vector<std::string> m_format_tags, m_subexprs;
    std::unique_ptr<BaseDatumFormatter<DataObject> > m_formatter;

  };

  template<typename DataObject>
  DatumFormatterAlias<DataObject> datum_formatter_alias(const std::string &_name, const std::string &_command) {
    return DatumFormatterAlias<DataObject>(_name, _command);
  }

  template<typename DataObject>
  DatumFormatterAlias<DataObject> datum_formatter_alias(const std::string &_name, const BaseDatumFormatter<DataObject> &_inside) {
    return DatumFormatterAlias<DataObject>(_name, _inside);
  }

  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  /* ConstantValueFormatter prints a string value specified at construction.  A header string can also be passed.
   */
  template<typename ValueType, typename DataObject>
  class ConstantValueFormatter: public BaseDatumFormatter<DataObject> {
  public:
    using BaseDatumFormatter<DataObject>::name;

    ConstantValueFormatter(const std::string &_header, const ValueType &_value, bool _print_json = false) :
      BaseDatumFormatter<DataObject>(_header, "Constant value that will be printed on each line"), m_value(_value), m_print_json(_print_json) {}

    BaseDatumFormatter<DataObject> *clone() const override {
      return new ConstantValueFormatter(*this);
    }

    void inject(const DataObject &_data_obj, DataStream &_stream, Index pass_index) const {
      _stream << m_value;
    }

    void print(const DataObject &_data_obj, std::ostream &_stream, Index pass_index) const override {
      _stream << m_value;
    }

    jsonParser &to_json(const DataObject &_data_obj, jsonParser &json)const override {
      if(m_print_json)
        json = m_value;
      return json;
    }

  private:
    ValueType m_value;
    bool m_print_json;
  };

  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  template<typename ValueType, typename DataObject>
  class BaseValueFormatter: public BaseDatumFormatter<DataObject> {
  public:
    using BaseDatumFormatter<DataObject>::name;

    BaseValueFormatter(const std::string &_header, const std::string &_comment) :
      BaseDatumFormatter<DataObject>(_header, _comment) {}


    void inject(const DataObject &_data_obj, DataStream &_stream, Index pass_index) const {
      _stream << _evaluate(_data_obj);
      if(!_validate(_data_obj))
        _stream << DataStream::failbit;
    }

    void print(const DataObject &_data_obj, std::ostream &_stream, Index pass_index) const {
      _stream.flags(std::ios::showpoint | std::ios::fixed | std::ios::right);
      _stream.precision(8);
      if(_validate(_data_obj))
        _stream << _evaluate(_data_obj);
      else
        _stream << "unknown";
    }

    jsonParser &to_json(const DataObject &_data_obj, jsonParser &json)const {
      if(_validate(_data_obj))
        json = _evaluate(_data_obj);
      return json;
    }
  protected:

    virtual bool _validate(const DataObject &_data_obj) const = 0;
    virtual ValueType _evaluate(const DataObject &_data_obj) const = 0;
  };

  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


  template<typename ValueType, typename DataObject>
  class GenericDatumFormatter : public BaseDatumFormatter<DataObject> {
  public:
    using BaseDatumFormatter<DataObject>::name;
    using Evaluator = std::function<ValueType(const DataObject &)>;
    using Validator = std::function<bool(const DataObject &)>;

    GenericDatumFormatter(const std::string &_init_name, const std::string &_desc, Evaluator evaluator, Validator validator = always_true<DataObject>) :
      BaseDatumFormatter<DataObject>(_init_name, _desc), m_evaluate(evaluator), m_validate(validator) {}

    virtual BaseDatumFormatter<DataObject> *clone() const override {
      return new GenericDatumFormatter(*this);
    }

    void inject(const DataObject &_data_obj, DataStream &_stream, Index pass_index) const override {
      _stream << m_evaluate(_data_obj);
      if(!validate(_data_obj))
        _stream << DataStream::failbit;
    }

    bool validate(const DataObject &_data_obj)const override {
      return m_validate(_data_obj);
    }

    void print(const DataObject &_data_obj, std::ostream &_stream, Index pass_index) const override {
      _stream.flags(std::ios::showpoint | std::ios::fixed | std::ios::right);
      _stream.precision(8);
      if(validate(_data_obj))
        _stream << m_evaluate(_data_obj);
      else
        _stream << "unknown";
    }

    jsonParser &to_json(const DataObject &_data_obj, jsonParser &json)const override {
      if(validate(_data_obj))
        json = m_evaluate(_data_obj);
      return json;
    }

  private:
    Evaluator m_evaluate;
    Validator m_validate;

  };

  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  template<typename _Container,
           typename DataObject,
           typename _value_type = typename _Container::value_type,
           typename _size_type = typename _Container::value_type,
           typename Access = CASM_TMP::BracketAccess<_Container, _value_type, _size_type> >
  class Generic1DDatumFormatter : public BaseDatumFormatter<DataObject> {
  public:
    using BaseDatumFormatter<DataObject>::name;
    using BaseDatumFormatter<DataObject>::_add_rule;
    using BaseDatumFormatter<DataObject>::_index_rules;
    using BaseDatumFormatter<DataObject>::_parse_index_expression;
    using Container = _Container;
    using Evaluator = std::function<Container(const DataObject &)>;
    using Sizer = std::function<Index(const Container &)>;
    using Validator = std::function<bool(const DataObject &)>;

    Generic1DDatumFormatter(const std::string &_init_name,
                            const std::string &_desc,
                            Evaluator evaluator,
                            Validator validator = always_true<DataObject>,
                            Sizer sizer =
    [](const Container &cont)->Index{
      return cont.size();
    }) :
      BaseDatumFormatter<DataObject>(_init_name, _desc), m_evaluate(evaluator), m_validate(validator), m_size(sizer) {}

    BaseDatumFormatter<DataObject> *clone() const override {
      return new Generic1DDatumFormatter(*this);
    }

    void init(const DataObject &_template_obj) const override;

    std::string long_header(const DataObject &_template_obj) const override;

    bool validate(const DataObject &_data_obj)const override {
      return m_validate(_data_obj);
    }

    bool parse_args(const std::string &args) override {
      _parse_index_expression(args);
      return true;
    }

    void inject(const DataObject &_data_obj, DataStream &_stream, Index pass_index) const override {
      Container val = m_evaluate(_data_obj);
      auto it(_index_rules().cbegin()), end_it(_index_rules().cend());
      if(!validate(_data_obj))
        _stream << DataStream::failbit;
      for(; it != end_it; ++it) {
        _stream << Access::at(val, (*it)[0]);
      }
    }

    void print(const DataObject &_data_obj, std::ostream &_stream, Index pass_index) const override {
      _stream.flags(std::ios::showpoint | std::ios::fixed | std::ios::right);
      _stream.precision(8);
      bool known = validate(_data_obj);
      Container val;
      if(known)
        val = m_evaluate(_data_obj);
      auto it(_index_rules().cbegin()), end_it(_index_rules().cend());
      for(; it != end_it; ++it) {
        if(known)
          _stream << "    " << Access::at(val, (*it)[0]);
        else
          _stream << "    unknown";

      }
    }

    jsonParser &to_json(const DataObject &_data_obj, jsonParser &json)const override {
      if(validate(_data_obj))
        json = m_evaluate(_data_obj);
      return json;
    }

  private:
    Evaluator m_evaluate;
    Validator m_validate;
    Sizer m_size;
  };

  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  template<typename _Container,
           typename DataObject,
           typename _value_type = typename _Container::value_type,
           typename _size_type = typename _Container::value_type,
           typename Access = CASM_TMP::BracketAccess<_Container, _value_type, _size_type> >
  class Generic2DDatumFormatter : public BaseDatumFormatter<DataObject> {
  public:
    using BaseDatumFormatter<DataObject>::name;
    using BaseDatumFormatter<DataObject>::_add_rule;
    using BaseDatumFormatter<DataObject>::_index_rules;
    using BaseDatumFormatter<DataObject>::_parse_index_expression;
    using Container = _Container;
    using Evaluator = std::function<Container(const DataObject &)>;
    using Sizer = std::function<std::vector<Index>(const Container &)>;
    using Validator = std::function<bool(const DataObject &)>;

    Generic2DDatumFormatter(const std::string &_init_name,
                            const std::string &_desc,
                            Evaluator evaluator,
                            Validator validator = always_true<DataObject>,
                            Sizer sizer =
    [](const Container &cont)->std::vector<Index> {
      std::vector<Index> tsize(2, 0);
      tsize[0] = cont.size();
      if(tsize[0] > 0)
        tsize[1] = cont[0].size();
      return tsize;
    }) :
      BaseDatumFormatter<DataObject>(_init_name, _desc), m_evaluate(evaluator), m_validate(validator), m_size(sizer) {}

    BaseDatumFormatter<DataObject> *clone() const override {
      return new Generic2DDatumFormatter(*this);
    }

    void init(const DataObject &_template_obj) const override;

    std::string long_header(const DataObject &_template_obj) const override;
    std::string short_header(const DataObject &_template_obj) const override;

    bool validate(const DataObject &_data_obj) const override {
      return m_validate(_data_obj);
    }

    bool parse_args(const std::string &args) override {
      _parse_index_expression(args);
      return true;
    }

    void inject(const DataObject &_data_obj, DataStream &_stream, Index pass_index) const override {
      Container val = m_evaluate(_data_obj);
      auto it(_index_rules().cbegin()), end_it(_index_rules().cend());
      if(!validate(_data_obj))
        _stream << DataStream::failbit;
      for(; it != end_it; ++it) {
        _stream << Access::at(val, (*it)[0], (*it)[1]);
      }
    }

    void print(const DataObject &_data_obj, std::ostream &_stream, Index pass_index) const override {
      _stream.flags(std::ios::showpoint | std::ios::fixed | std::ios::right);
      _stream.precision(8);
      bool known = validate(_data_obj);
      Container val;
      if(known)
        val = m_evaluate(_data_obj);
      auto it(_index_rules().cbegin()), end_it(_index_rules().cend());
      for(; it != end_it; ++it) {
        if(known)
          _stream << "    " << Access::at(val, (*it)[0], (*it)[1]);
        else
          _stream << "    unknown";
      }
    }

    jsonParser &to_json(const DataObject &_data_obj, jsonParser &json)const override {
      if(validate(_data_obj))
        json = m_evaluate(_data_obj);
      return json;
    }

  private:
    Evaluator m_evaluate;
    Validator m_validate;
    Sizer m_size;
  };


}

#include "casm/casm_io/DataFormatterTools_impl.hh"
#endif
