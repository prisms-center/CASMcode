#ifndef CASM_DataFormatterTools
#define CASM_DataFormatterTools
#include <numeric>
#include <iterator>
#include <boost/regex.hpp>
#include <boost/algorithm/string.hpp>
#include "casm/casm_io/dataformatter/DataFormatter.hh"
#include "casm/casm_io/container/json_io.hh"
#include "casm/container/ContainerTraits.hh"
#include "casm/misc/CASM_TMP.hh"

namespace CASM {

  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  /** \defgroup DataFormatterOperator

      \brief Operators on other DatumFormatters

      \ingroup DataFormatter
  */


  /// \brief DataFormatters that operate on the results of other DataFormatters
  ///
  /// \ingroup DataFormatterOperator
  ///
  template<typename ValueType, typename ArgType, typename DataObject>
  class DataFormatterOperator : public BaseDatumFormatter<DataObject> {
  public:

    using BaseDatumFormatter<DataObject>::name;
    using Evaluator = std::function<ValueType(const std::vector<ArgType> &)>;
    //validator probably not necessary?
    //using Validator = std::function<bool(const DataObject &)>;

    DataFormatterOperator(
      const std::string &_init_name,
      const std::string &_desc,
      Evaluator evaluator) :
      BaseDatumFormatter<DataObject>(_init_name, _desc), m_evaluate(evaluator) {}

    std::unique_ptr<DataFormatterOperator> clone() const {
      return std::unique_ptr<DataFormatterOperator>(this->_clone());
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
      if(!validate(_data_obj)) {
        _stream << DataStream::failbit << ValueType();
      }
      else {
        _stream << _evaluate(_data_obj);
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
      return m_evaluate(vec_stream.vector());
    }
  private:

    DataFormatterOperator *_clone() const override {
      return new DataFormatterOperator(*this);
    }

    Evaluator m_evaluate;
    DataFormatter<DataObject> m_arg_formatter;
  };



  /// \brief Makes a DataFormatterOperator that adds two or more numbers
  ///
  /// \ingroup DataFormatterOperator
  ///
  template<typename DataObject>
  DataFormatterOperator<double, double, DataObject> format_operator_add() {
    return DataFormatterOperator<double, double, DataObject>("add", "Add two or more numbers",
    [](const std::vector<double> &vec)->double {
      return std::accumulate(vec.cbegin(),
                             vec.cend(),
                             0.0);
    });

  }

  /// \brief Makes a DataFormatterOperator that subtracts two numbers
  ///
  /// \ingroup DataFormatterOperator
  ///
  template<typename DataObject>
  DataFormatterOperator<double, double, DataObject> format_operator_sub() {
    return DataFormatterOperator<double, double, DataObject>("sub", "Subtract two numbers",
    [](const std::vector<double> &vec)->double {
      if(vec.size() != 2)
        throw std::runtime_error("Subtraction operator must receive exactly two values!");
      return vec[0] - vec[1];
    });

  }

  /// \brief Makes a DataFormatterOperator that returns the product of two or more numbers
  ///
  /// \ingroup DataFormatterOperator
  ///
  template<typename DataObject>
  DataFormatterOperator<double, double, DataObject> format_operator_mult() {
    return DataFormatterOperator<double, double, DataObject>("mult", "Multiply two or more numbers",
    [](const std::vector<double> &vec)->double {
      return std::accumulate(vec.cbegin(),
                             vec.cend(),
                             1.0,
      [](double a, double b)->double{
        return a *b;
      });
    });

  }

  /// \brief Makes a DataFormatterOperator that divides two numbers
  ///
  /// \ingroup DataFormatterOperator
  ///
  template<typename DataObject>
  DataFormatterOperator<double, double, DataObject> format_operator_div() {
    return DataFormatterOperator<double, double, DataObject>("div", "Divide two numbers",
    [](const std::vector<double> &vec)->double {
      if(vec.size() != 2)
        throw std::runtime_error("Division operator must receive exactly two values!");
      return vec[0] / vec[1];
    });

  }

  /// \brief Makes a DataFormatterOperator that returns the root-mean-square value of 0 or more elements
  ///
  /// \ingroup DataFormatterOperator
  ///
  template<typename DataObject>
  DataFormatterOperator<double, double, DataObject> format_operator_rms() {
    return DataFormatterOperator<double, double, DataObject>("rms", "Root mean square of 0 or more numerical values",
    [](const std::vector<double> &vec)->double {
      return sqrt(std::accumulate(vec.cbegin(),
                                  vec.cend(),
                                  0.0,
      [](double a, double b)->double{
        return a + b *b;
      }));
    });
  }

  /// \brief Makes a DataFormatterOperator that returns the root-mean-square value of 0 or more elements
  ///
  /// \ingroup DataFormatterOperator
  ///
  template<typename DataObject>
  DataFormatterOperator<double, double, DataObject> format_operator_pnorm() {
    return DataFormatterOperator<double, double, DataObject>("pnorm", "Vector p-norm of zero or more elements Ex: 'pnorm(p, elem1, elem2)' evaluates (elem1^p+elem2^p)^(1/p)",
    [](const std::vector<double> &vec)->double {
      if(vec.empty())
        throw std::runtime_error("pnorm query operator must receive at least one value!");
      double p = vec[0];
      return pow(std::accumulate(++vec.cbegin(),
                                 vec.cend(),
                                 0.0,
      [p](double a, double b)->double{
        return a + pow(b, p);
      }), (1 / p));
    });
  }

  /// \brief Makes a DataFormatterOperator that returns the maximum of two or more numbers
  ///
  /// \ingroup DataFormatterOperator
  ///
  template<typename DataObject>
  DataFormatterOperator<double, double, DataObject> format_operator_max() {
    return DataFormatterOperator<double, double, DataObject>("max", "Max value of two or more numbers",
    [](const std::vector<double> &vec)->double {
      return (*std::max_element(vec.cbegin(),
                                vec.cend()));
    });
  }

  /// \brief Makes a DataFormatterOperator that returns the minimum of two or more numbers
  ///
  /// \ingroup DataFormatterOperator
  ///
  template<typename DataObject>
  DataFormatterOperator<double, double, DataObject> format_operator_min() {
    return DataFormatterOperator<double, double, DataObject>("min", "Min value of two or more numbers",
    [](const std::vector<double> &vec)->double {
      return (*std::min_element(vec.cbegin(),
                                vec.cend()));
    });
  }

  /// \brief Makes a DataFormatterOperator that returns the index (from 0) of the maximum of two or more numbers
  ///
  /// \ingroup DataFormatterOperator
  ///
  template<typename DataObject>
  DataFormatterOperator<long, double, DataObject> format_operator_imax() {
    return DataFormatterOperator<long, double, DataObject>("imax", "Index (from 0) of max value in array of two or more numbers",
    [](const std::vector<double> &vec)->long {
      auto it = std::max_element(vec.cbegin(), vec.cend());
      return std::distance(vec.cbegin(), it);
    });
  }

  /// \brief Makes a DataFormatterOperator that returns the index (from 0) of the minimum of two or more numbers
  ///
  /// \ingroup DataFormatterOperator
  ///
  template<typename DataObject>
  DataFormatterOperator<long, double, DataObject> format_operator_imin() {
    return DataFormatterOperator<long, double, DataObject>("imin", "Index (from 0) of min value in array of two or more numbers",
    [](const std::vector<double> &vec)->long {
      auto it = std::min_element(vec.cbegin(), vec.cend());
      return std::distance(vec.cbegin(), it);
    });
  }

  /// \brief Makes a DataFormatterOperator that returns the exponential of a number
  ///
  /// \ingroup DataFormatterOperator
  ///
  template<typename DataObject>
  DataFormatterOperator<double, double, DataObject> format_operator_exp() {
    return DataFormatterOperator<double, double, DataObject>("exp", "Exponential function",
    [](const std::vector<double> &vec)->double {
      if(vec.size() != 1)
        throw std::runtime_error("Exponent operator must receive exactly one value!");
      return exp(vec[0]);
    });

  }


  /// \brief Makes a DataFormatterOperator that checks if a string matches a regular expression.
  ///
  /// Ex: re('input_string','regex_pattern')"
  ///
  /// \ingroup DataFormatterOperator
  ///
  template<typename DataObject>
  DataFormatterOperator<bool, std::string, DataObject> format_operator_re() {
    return DataFormatterOperator<bool, std::string, DataObject>("re", "Check if string matches regular expression. Ex: re('input_string','regex_pattern')",
    [](const std::vector<std::string> &vec)->bool {
      if(vec.size() != 2)
        throw std::runtime_error("Operator re('input_string','regex_pattern') must receive exactly 2 values!");
      boost::regex e(boost::trim_copy_if(vec[1], boost::is_any_of(" '")));
      return boost::regex_match(vec[0], e);
    });
  }

  /// \brief Makes a DataFormatterOperator that checks if a string contains a regular expression.
  ///
  /// Ex: rs('input_string','regex_pattern')
  ///
  /// \ingroup DataFormatterOperator
  ///
  template<typename DataObject>
  DataFormatterOperator<bool, std::string, DataObject> format_operator_rs() {
    return DataFormatterOperator<bool, std::string, DataObject>("rs", "Check if string contains regular expression. Ex: rs('input_string','regex_pattern')",
    [](const std::vector<std::string> &vec)->bool {
      if(vec.size() != 2)
        throw std::runtime_error("Operator re('input_string','regex_pattern') must receive exactly 2 values!");
      boost::regex e(boost::trim_copy_if(vec[1], boost::is_any_of(" '")));
      return boost::regex_search(vec[0], e);
    });
  }

  /// \brief Makes a DataFormatterOperator that returns the square of a number
  ///
  /// \ingroup DataFormatterOperator
  ///
  template<typename DataObject>
  DataFormatterOperator<double, double, DataObject> format_operator_sq() {
    return DataFormatterOperator<double, double, DataObject>("sq", "Square of a number",
    [](const std::vector<double> &vec)->double {
      if(vec.size() != 1)
        throw std::runtime_error("Square operator must receive exactly one value!");
      return vec[0] * vec[0];
    });

  }

  /// \brief Makes a DataFormatterOperator that returns the square root of a number
  ///
  /// \ingroup DataFormatterOperator
  ///
  template<typename DataObject>
  DataFormatterOperator<double, double, DataObject> format_operator_sqrt() {
    return DataFormatterOperator<double, double, DataObject>("sqrt", "Square root of a number",
    [](const std::vector<double> &vec)->double {
      if(vec.size() != 1)
        throw std::runtime_error("Square-root operator must receive exactly one value!");
      return sqrt(vec[0]);
    });

  }

  /// \brief Makes a DataFormatterOperator that returns the negative of a number
  ///
  /// \ingroup DataFormatterOperator
  ///
  template<typename DataObject>
  DataFormatterOperator<double, double, DataObject> format_operator_neg() {
    return DataFormatterOperator<double, double, DataObject>("neg", "Negative of a number",
    [](const std::vector<double> &vec)->double {
      if(vec.size() != 1)
        throw std::runtime_error("Negation operator must receive exactly one value!");
      return -vec[0];
    });

  }

  /// \brief Makes a DataFormatterOperator that returns the boolean AND for a sequence of boolean values
  ///
  /// \ingroup DataFormatterOperator
  ///
  template<typename DataObject>
  DataFormatterOperator<bool, bool, DataObject> format_operator_and() {
    return DataFormatterOperator<bool, bool, DataObject>("and", "Boolean AND for sequence of boolean values",
    [](const std::vector<bool> &vec)->bool {
      return std::accumulate(vec.cbegin(),
                             vec.cend(),
                             true,
      [](bool a, bool b)->bool{
        return a && b;
      });
    });

  }

  /// \brief Makes a DataFormatterOperator that returns the boolean OR for a sequence of boolean values
  ///
  /// \ingroup DataFormatterOperator
  ///
  template<typename DataObject>
  DataFormatterOperator<bool, bool, DataObject> format_operator_or() {
    return DataFormatterOperator<bool, bool, DataObject>("or", "Boolean OR for sequence of boolean values",
    [](const std::vector<bool> &vec)->bool {
      return std::accumulate(vec.cbegin(),
                             vec.cend(),
                             false,
      [](bool a, bool b)->bool{
        return a || b;
      });
    });

  }

  /// \brief Makes a DataFormatterOperator that returns the boolean NOT for a single boolean value
  ///
  /// \ingroup DataFormatterOperator
  ///
  template<typename DataObject>
  DataFormatterOperator<bool, bool, DataObject> format_operator_xor() {
    return DataFormatterOperator<bool, bool, DataObject>("xor", "Boolean XOR for for two boolean values",
    [](const std::vector<bool> &vec)->bool {
      if(vec.size() != 2)
        throw std::runtime_error("Boolean XOR operator expects exactly two values!");
      return (vec[0] && !vec[1]) || (!vec[0] && vec[1]);
    });
  }

  template<typename DataObject>
  DataFormatterOperator<bool, bool, DataObject> format_operator_not() {
    return DataFormatterOperator<bool, bool, DataObject>("not", "Boolean NOT for a single boolean value",
    [](const std::vector<bool> &vec)->bool {
      if(vec.size() != 1)
        throw std::runtime_error("Boolean NOT operator must receive exactly one value!");
      return !vec[0];
    });

  }

  /// \brief Makes a DataFormatterOperator for equality comparison of two numbers
  ///
  /// \ingroup DataFormatterOperator
  ///
  template<typename DataObject>
  DataFormatterOperator<bool, double, DataObject> format_operator_eq() {
    return DataFormatterOperator<bool, double, DataObject>("eq", "Equality comparison for two values",
    [](const std::vector<double> &vec)->bool {
      if(vec.size() != 2)
        throw std::runtime_error("Equality operator must receive exactly two values!");
      return almost_equal(vec[0], vec[1]);
    });

  }

  /// \brief Makes a DataFormatterOperator for less-than comparison of two numbers
  ///
  /// \ingroup DataFormatterOperator
  ///
  template<typename DataObject>
  DataFormatterOperator<bool, double, DataObject> format_operator_lt() {
    return DataFormatterOperator<bool, double, DataObject>("lt", "Less-than comparison for two values",
    [](const std::vector<double> &vec)->bool {
      if(vec.size() != 2)
        throw std::runtime_error("Less-than operator must receive exactly two values!");
      return vec[0] < vec[1];
    });

  }

  /// \brief Makes a DataFormatterOperator for less-than-or-equal comparison of two numbers
  ///
  /// \ingroup DataFormatterOperator
  ///
  template<typename DataObject>
  DataFormatterOperator<bool, double, DataObject> format_operator_le() {
    return DataFormatterOperator<bool, double, DataObject>("le", "Less-than-or-equal comparison for two values",
    [](const std::vector<double> &vec)->bool {
      if(vec.size() != 2)
        throw std::runtime_error("Less-than-or-equal operator must receive exactly two values!");
      return vec[0] <= vec[1];
    });

  }

  /// \brief Makes a DataFormatterOperator for greater-than comparison of two numbers
  ///
  /// \ingroup DataFormatterOperator
  ///
  template<typename DataObject>
  DataFormatterOperator<bool, double, DataObject> format_operator_gt() {
    return DataFormatterOperator<bool, double, DataObject>("gt", "Greater-than comparison for two values",
    [](const std::vector<double> &vec)->bool {
      if(vec.size() != 2)
        throw std::runtime_error("Greater-than operator must receive exactly two values!");
      return vec[0] > vec[1];
    });

  }

  /// \brief Makes a DataFormatterOperator for greater-than-or-equal comparison of two numbers
  ///
  /// \ingroup DataFormatterOperator
  ///
  template<typename DataObject>
  DataFormatterOperator<bool, double, DataObject> format_operator_ge() {
    return DataFormatterOperator<bool, double, DataObject>("ge", "Greater-than-or-equal comparison for two values",
    [](const std::vector<double> &vec)->bool {
      if(vec.size() != 2)
        throw std::runtime_error("Greater-than-or-equal operator must receive exactly two values!");
      return vec[0] >= vec[1];
    });

  }

  /// \brief Dictionary of all DatumFormatterOperator
  template<typename DataObject>
  DataFormatterDictionary<DataObject> make_operator_dictionary() {
    DataFormatterDictionary<DataObject> dict;

    dict.insert(
      format_operator_add<DataObject>(),
      format_operator_sub<DataObject>(),
      format_operator_mult<DataObject>(),
      format_operator_div<DataObject>(),
      format_operator_exp<DataObject>(),
      format_operator_sq<DataObject>(),
      format_operator_sqrt<DataObject>(),
      format_operator_neg<DataObject>(),
      format_operator_and<DataObject>(),
      format_operator_or<DataObject>(),
      format_operator_not<DataObject>(),
      format_operator_min<DataObject>(),
      format_operator_max<DataObject>(),
      format_operator_imin<DataObject>(),
      format_operator_imax<DataObject>(),
      format_operator_eq<DataObject>(),
      format_operator_lt<DataObject>(),
      format_operator_le<DataObject>(),
      format_operator_gt<DataObject>(),
      format_operator_ge<DataObject>(),
      format_operator_re<DataObject>(),
      format_operator_rs<DataObject>()
    );

    return dict;
  }


  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  /// \brief Implements a DatumFormatter that is an alias for a combination of others
  ///
  /// \ingroup DataFormatterTypes
  ///
  template<typename DataObject>
  class DatumFormatterAlias: public BaseDatumFormatter<DataObject> {
  public:
    using BaseDatumFormatter<DataObject>::name;

    DatumFormatterAlias(const std::string &_name,
                        const std::string &_command,
                        const DataFormatterDictionary<DataObject> &_dict,
                        const std::string &_help = "") :
      BaseDatumFormatter<DataObject> (_name, _help.empty() ? "User-specified alias for '" + _command + "'" : _help) {

      split_formatter_expression(_command, m_format_tags, m_subexprs);
      if(m_format_tags.size() != 1)
        throw std::runtime_error("Expression '" + _command + "' is either empty or consists of multiple expressions.\n");
      m_formatter = _dict.lookup(m_format_tags[0])->clone();
    }

    DatumFormatterAlias(const std::string &_name, const BaseDatumFormatter<DataObject> &_rhs, const std::string &_help = "") :
      BaseDatumFormatter<DataObject> (_name,  _help.empty() ? "User-specified alias for '" + _rhs.name() + "'" : _help),
      m_formatter(_rhs.clone()) {}

    typename BaseDatumFormatter<DataObject>::FormatterType type() const override {
      return BaseDatumFormatter<DataObject>::Property;
    }

    std::unique_ptr<DatumFormatterAlias> clone() const {
      return std::unique_ptr<DatumFormatterAlias>(this->_clone());
    }

    bool init(const DataObject &_template_obj) const  override {
      return m_formatter->init(_template_obj);
    }

    ///\brief Returns true if _data_obj has valid values for requested data
    bool validate(const DataObject &_data_obj) const  override {
      return m_formatter->validate(_data_obj);
    }

    ///\brief Returns a std::vector<std::string> with each column header
    ///
    /// - Default returns: `$NAME` if only 1 column, `$NAME(i)` if >1 column
    std::vector<std::string> col_header(const DataObject &_template_obj) const override {
      std::vector<std::string> _col;
      CountDataStream tcount;
      m_formatter->inject(_template_obj, tcount);
      if(tcount.count() == 1) {
        _col.push_back(name());
        return _col;
      }

      for(Index i = 0; i < tcount.count(); i++) {
        std::stringstream t_ss;
        t_ss << name() << '(' << i << ')';
        _col.push_back(t_ss.str());
      }

      return _col;
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

      if(!m_formatter) {
        throw std::runtime_error("ERROR: DataFormatterAlias has no formatter");
      }
      //Parse the arguments of of the command now.  Later, we may want to do expression substitution
      // (e.g.,  "comp_plus = add(comp(a),$1)" where $1 specifies an argument)
      if(m_subexprs.size() && !m_formatter->parse_args(m_subexprs[0])) {
        throw std::runtime_error("Invalid arguments passed to '" + m_format_tags[0] + "'. Cannot accept expression '" + m_subexprs[0] + "'\n");
      }

      return args.size() == 0;
      //return m_formatter->parse_args(args);
    }
  private:

    /// \brief Make an exact copy of the formatter (including any initialized members)
    ///
    DatumFormatterAlias *_clone() const override {
      return new DatumFormatterAlias(*this);
    }

    std::vector<std::string> m_format_tags, m_subexprs;
    notstd::cloneable_ptr<BaseDatumFormatter<DataObject> > m_formatter;

  };

  /// \brief Make a DatumFormatterAlias
  ///
  /// \ingroup DataFormatter
  ///
  template<typename DataObject>
  DatumFormatterAlias<DataObject> datum_formatter_alias(
    const std::string &_name,
    const std::string &_command,
    const DataFormatterDictionary<DataObject> &_dict,
    const std::string &_help = "") {
    return DatumFormatterAlias<DataObject>(_name, _command, _dict, _help);
  }

  /// \brief Make a DatumFormatterAlias
  ///
  /// \ingroup DataFormatter
  ///
  template<typename DataObject>
  DatumFormatterAlias<DataObject> datum_formatter_alias(
    const std::string &_name,
    const BaseDatumFormatter<DataObject> &_inside,
    const std::string &_help = "") {
    return DatumFormatterAlias<DataObject>(_name, _inside, _help);
  }

  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  /// \brief Prints a string value specified at construction.  A header string can also be passed.
  ///
  /// \ingroup DataFormatterTypes
  ///
  template<typename ValueType, typename DataObject>
  class ConstantValueFormatter: public BaseDatumFormatter<DataObject> {
  public:
    using BaseDatumFormatter<DataObject>::name;

    ConstantValueFormatter(const std::string &_header, const ValueType &_value, bool _print_json = false) :
      BaseDatumFormatter<DataObject>(_header, "Constant value that will be printed on each line"), m_value(_value), m_print_json(_print_json) {}

    std::unique_ptr<ConstantValueFormatter> clone() const {
      return std::unique_ptr<ConstantValueFormatter>(this->_clone());
    }

    void inject(const DataObject &_data_obj, DataStream &_stream, Index pass_index) const override {
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

    ConstantValueFormatter *_clone() const override {
      return new ConstantValueFormatter(*this);
    }


    ValueType m_value;
    bool m_print_json;
  };

  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


  /// \brief Base class for creating scalar DatumFormatter
  ///
  /// A scalar DatumFormatter accepts an object of type DataObject and returns a
  /// scalar object of ValueType.
  ///
  /// To implement MyDatumFormatter, you must:
  /// - Write the MyDatumFormatter class:
  /// \code
  /// class MyDatumFormatter : public BaseValueFormatter<ValueType, DataObject> {
  ///   ... include implementations here ...
  /// };
  /// \endcode
  ///
  /// The MyDatumFormatter class is required to provide an implementation for:
  /// - ValueType evaluate(const DataObject& obj) const;
  /// - BaseDatumFormatter<DataObject> *clone() const;
  ///
  /// See base BaseDatumFormatter for default implementations which may be specialized.
  ///
  /// \seealso GenericDatumFormatter to specialize a scalar DatumFormatter at runtime
  ///
  /// \ingroup DataFormatterTypes
  ///
  template<typename ValueType, typename DataObject>
  class BaseValueFormatter :
    public BaseDatumFormatter<DataObject> {

  public:

    /// \brief Constructor
    BaseValueFormatter(const std::string &_name, const std::string &_desc) :
      BaseDatumFormatter<DataObject>(_name, _desc) {}

    /// \brief Destructor
    virtual ~BaseValueFormatter() {}

    std::unique_ptr<BaseValueFormatter> clone() const {
      return std::unique_ptr<BaseValueFormatter>(this->_clone());
    }

    /// \brief Return requested data from obj, throwing std::runtime_error if not valid
    virtual ValueType operator()(const DataObject &obj) const {
      if(!this->validate(obj)) {
        throw std::runtime_error(std::string("Invalid DataObject in ") + this->name());
      }
      return evaluate(obj);
    }


    // --- Derived classes require an evaluate member ----

    virtual ValueType evaluate(const DataObject &obj) const = 0;


    // --- These methods specialize virtual BaseDatumFormatter<DataObject> methods ----

    /// \brief Default implementation injects each element, via operator<<
    ///
    /// - sets DataStream::failbit if validation fails
    virtual void inject(const DataObject &_data_obj, DataStream &_stream, Index pass_index = 0) const override {
      if(!this->validate(_data_obj))
        _stream << DataStream::failbit << ValueType();
      else
        _stream << this->evaluate(_data_obj);
    }

    /// \brief Default implementation prints each element in a column, via operator<<
    ///
    /// - Prints "unknown" if validation fails
    virtual void print(const DataObject &_data_obj, std::ostream &_stream, Index pass_index = 0) const override {
      _stream.flags(std::ios::showpoint | std::ios::fixed | std::ios::right);
      _stream.precision(8);
      if(this->validate(_data_obj))
        _stream << this->evaluate(_data_obj);
      else
        _stream << "unknown";
    }

    /// \brief Default implementation calls jsonParser& to_json(const ValueType&, jsonParser&)
    ///
    /// - Does nothing if validation fails
    virtual jsonParser &to_json(const DataObject &_data_obj, jsonParser &json) const override {
      if(this->validate(_data_obj))
        json = this->evaluate(_data_obj);
      return json;
    }

  private:

    /// \brief Clone
    virtual BaseValueFormatter *_clone() const override = 0;
  };


  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


  /// \brief A DatumFormatter that returns a value of specified type, via functions
  ///        that may be specified at runtime
  ///
  /// \ingroup DataFormatterTypes
  ///
  template<typename ValueType, typename DataObject>
  class GenericDatumFormatter : public BaseValueFormatter<ValueType, DataObject> {

  public:

    using Evaluator = std::function<ValueType(const DataObject &)>;
    using Validator = std::function<bool(const DataObject &)>;

    /// \brief Constructor
    ///
    /// \param _name Name of formatter
    /// \param _desc Description of the formatter
    /// \param _evaluator Returns a ValueType type piece of data about a DataObject
    /// \param _validator Returns boolean indicating if the _evaluator can be used successfully
    ///
    GenericDatumFormatter(const std::string &_init_name,
                          const std::string &_desc,
                          Evaluator _evaluator,
                          Validator _validator = always_true<DataObject>) :
      BaseValueFormatter<ValueType, DataObject>(_init_name, _desc),
      m_evaluate(_evaluator),
      m_validate(_validator) {}

    // /// \brief Destructor
    // virtual ~GenericDatumFormatter() {}

    // --- Required implementations -----------

    /// \brief Returns Container type piece of data about a DataObject, the result of the Evaluator
    ValueType evaluate(const DataObject &obj) const override {
      return m_evaluate(obj);
    }

    /// \brief Clone
    std::unique_ptr<GenericDatumFormatter> clone() const {
      return std::unique_ptr<GenericDatumFormatter>(this->_clone());
    }


    // --- Specialized implementation -----------

    /// \brief Returns boolean indicating if the _evaluator can be used successfully,
    ///        the result of the Validator
    bool validate(const DataObject &obj) const override {
      return m_validate(obj);
    }


  private:

    /// \brief Clone
    GenericDatumFormatter *_clone() const override {
      return new GenericDatumFormatter(*this);
    }

    Evaluator m_evaluate;
    Validator m_validate;

  };

  template<typename DataObject>
  GenericDatumFormatter<std::string, DataObject> name() {
    return GenericDatumFormatter<std::string, DataObject>(
             "name",
             traits<DataObject>::name + " name",
    [](const DataObject & obj)->std::string {
      return obj.name();
    });
  }

  template<typename DataObject>
  GenericDatumFormatter<std::string, DataObject> alias() {
    return GenericDatumFormatter<std::string, DataObject>(
             "alias",
             traits<DataObject>::name + " alias (if exists, else \"none\")",
    [](const DataObject & obj)->std::string {
      std::string alias = obj.alias();
      if(alias.empty()) {
        return "none";
      }
      return alias;
    });
  }

  template<typename DataObject>
  GenericDatumFormatter<std::string, DataObject> alias_or_name() {
    return GenericDatumFormatter<std::string, DataObject>(
             "alias_or_name",
             traits<DataObject>::name + " alias (if exists), else name",
    [](const DataObject & obj)->std::string {
      std::string alias = obj.alias();
      if(alias.empty()) {
        return obj.name();
      }
      return alias;
    });
  }

  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  /// \brief Abstract base class for creating 1D DatumFormatter
  ///
  /// A 1D DatumFormatter accepts an object of type DataObject and returns a
  /// one-dimensional Container object.
  ///
  /// To implement My1DDatumFormatter, you must:
  /// - Write the My1DDatumFormatter class:
  /// \code
  /// class My1DDatumFormatter : public Base1dDatumFormatter<DataObject, Container> {
  ///   ... include implementations here ...
  /// };
  /// \endcode
  /// - Specialize the ContainerTraits class if it does not already exist:
  /// \code
  /// template<>
  /// struct ContainerTraits<Container> {
  ///   ... include specializations here ...
  /// };
  /// \endcode
  ///
  /// The ContainerTraits specialization requires specifying:
  /// - typename Container
  /// - typename value_type
  /// - typename size_type
  /// - typename Access // such as CASM_TMP::BracketAccess or CASM_TMP::ParenthesesAccess
  /// - size_type size(const Container&) const;
  ///
  /// The My1DDatumFormatter class is required to provide an implementation for:
  /// - Container evaluate(const DataObject& obj) const;
  /// - BaseDatumFormatter<DataObject> *clone() const;
  ///
  /// See BaseDatumFormatter and BaseValueFormatter for default implementations
  /// which may be specialized.
  ///
  /// \seealso Generic1DDatumFormatter to specialize a 1D DatumFormatter at runtime.
  ///
  /// \ingroup DataFormatterTypes
  ///
  template<typename Container, typename DataObject>
  class Base1DDatumFormatter :
    public BaseValueFormatter<Container, DataObject>,
    public ContainerTraits<Container> {

  public:

    /// \brief Access methods for Container
    typedef typename ContainerTraits<Container>::Access Access;
    typedef typename ContainerTraits<Container>::value_type ValueType;


    /// \brief Constructor
    Base1DDatumFormatter(const std::string &_name, const std::string &_desc) :
      BaseValueFormatter<Container, DataObject>(_name, _desc) {}

    /// \brief Destructor
    virtual ~Base1DDatumFormatter() {}

    /// \brief Clone
    std::unique_ptr<Base1DDatumFormatter> clone() const {
      return std::unique_ptr<Base1DDatumFormatter>(this->_clone());
    }


    // --- These methods specialize virtual BaseValueFormatter<Container, DataObject> methods ----


    /// \brief Default initialization adds rules for each element
    virtual bool init(const DataObject &_template_obj) const override {
      if(_index_rules().size())
        return true;

      if(!this->validate(_template_obj))
        return false;

      Index size = ContainerTraits<Container>::size(this->evaluate(_template_obj));
      for(Index i = 0; i < size; i++) {
        _add_rule(std::vector<Index>({i}));
      }
      return true;
    }

    /// \brief Default col_header uses 'name(index)' for each column
    ///
    /// Ex: "corr(0)" "corr(1)" "corr(5)" "corr(6)"
    virtual std::vector<std::string> col_header(const DataObject &_template_obj) const override {
      std::vector<std::string> _col;
      auto it(_index_rules().cbegin()), end_it(_index_rules().cend());
      Index s = max(8 - int(this->name().size()), 0);
      for(; it != end_it; ++it) {
        std::stringstream t_ss;
        t_ss << std::string(s, ' ') << this->name() << '(' << (*it)[0] << ')';
        _col.push_back(t_ss.str());
      }
      return _col;
    }

    /// \brief Default implementation calls _parse_index_expression
    virtual bool parse_args(const std::string &args) override {
      _parse_index_expression(args);
      return true;
    }

    /// \brief Default implementation injects each element
    ///
    /// - sets DataStream::failbit if validation fails
    virtual void inject(const DataObject &_data_obj, DataStream &_stream, Index pass_index = 0) const override {

      // add_rules to print all elements if not set yet
      if(!_index_rules().size()) {
        init(_data_obj);
      }

      Container val = this->evaluate(_data_obj);
      auto it(_index_rules().cbegin()), end_it(_index_rules().cend());
      if(!this->validate(_data_obj)) {
        _stream << DataStream::failbit;
        for(; it != end_it; ++it) {
          _stream << ValueType();
        }
      }
      else {
        for(; it != end_it; ++it) {
          _stream << Access::at(val, (*it)[0]);
        }
      }
    }

    /// \brief Default implementation prints each element in a column
    ///
    /// - Prints "unknown" if validation fails
    virtual void print(const DataObject &_data_obj, std::ostream &_stream, Index pass_index = 0) const override {

      // add_rules to print all elements if not set yet
      if(!_index_rules().size()) {
        init(_data_obj);
      }

      _stream.flags(std::ios::showpoint | std::ios::fixed | std::ios::right);
      _stream.precision(8);
      bool known = this->validate(_data_obj);
      Container val;
      if(known)
        val = this->evaluate(_data_obj);
      auto it(_index_rules().cbegin()), end_it(_index_rules().cend());
      for(; it != end_it; ++it) {
        if(known)
          _stream << "  " << Access::at(val, (*it)[0]);
        else
          _stream << "  unknown";

      }
    }

  protected:

    using BaseDatumFormatter<DataObject>::_add_rule;
    using BaseDatumFormatter<DataObject>::_index_rules;
    using BaseDatumFormatter<DataObject>::_parse_index_expression;

  private:

    /// \brief Clone
    virtual Base1DDatumFormatter *_clone() const override = 0;

  };


  /// \brief A DatumFormatter that returns a 1D value of specified type, via functions
  ///        that may be specified at runtime
  ///
  /// \ingroup DataFormatterTypes
  ///
  template<typename Container, typename DataObject>
  class Generic1DDatumFormatter : public Base1DDatumFormatter<Container, DataObject> {

  public:

    typedef std::function<Container(const DataObject &)> Evaluator;
    typedef std::function<bool (const DataObject &)> Validator;

    /// \brief Constructor
    ///
    /// \param _name Name of formatter
    /// \param _desc Description of the formatter
    /// \param _evaluator Returns a Container type piece of data about a DataObject
    /// \param _validator Returns boolean indicating if the _evaluator can be used successfully
    ///
    Generic1DDatumFormatter(const std::string &_name,
                            const std::string &_desc,
                            Evaluator _evaluator,
                            Validator _validator = always_true<DataObject>) :
      Base1DDatumFormatter<Container, DataObject>(_name, _desc),
      m_evaluate(_evaluator),
      m_validate(_validator) {}


    // --- Required implementations -----------

    /// \brief Returns Container type piece of data about a DataObject, the result of the Evaluator
    Container evaluate(const DataObject &obj) const override {
      return m_evaluate(obj);
    }

    /// \brief Clone using copy constructor
    std::unique_ptr<Generic1DDatumFormatter> clone() const {
      return std::unique_ptr<Generic1DDatumFormatter>(this->_clone());
    }


    // --- Specialized implementation -----------

    /// \brief Returns boolean indicating if the _evaluator can be used successfully,
    ///        the result of the Validator
    bool validate(const DataObject &obj) const override {
      return m_validate(obj);
    }

  private:

    /// \brief Clone using copy constructor
    Generic1DDatumFormatter *_clone() const override {
      return new Generic1DDatumFormatter(*this);
    }

    Evaluator m_evaluate;
    Validator m_validate;
  };

  /// \brief Dictionary of all Attribute
  ///
  /// Default includes String, Boolean, Integer, Scalar, and VectorXd attributes
  template<typename DataObject>
  DataFormatterDictionary<DataObject> make_attribute_dictionary() {
    DataFormatterDictionary<DataObject> dict;

    dict.insert(
      make_string_dictionary<DataObject>(),
      make_boolean_dictionary<DataObject>(),
      make_integer_dictionary<DataObject>(),
      make_scalar_dictionary<DataObject>(),
      make_vectorxi_dictionary<DataObject>(),
      make_vectorxd_dictionary<DataObject>(),
      make_matrixxd_dictionary<DataObject>()
    );

    return dict;
  }


  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


  /// \brief A DatumFormatter that returns a value of specified 2d container
  ///
  /// Functionality is specified via template parameters at compile time
  ///
  /// \ingroup DataFormatter
  ///
  template<typename Container, typename DataObject>
  class Base2DDatumFormatter :
    public BaseValueFormatter<Container, DataObject> {
  public:

    using typename BaseDatumFormatter<DataObject>::difference_type;
    /// \brief Access methods for Container
    typedef typename ContainerTraits<Container>::Access Access;
    typedef typename ContainerTraits<Container>::value_type2D ValueType;


    /// \brief Constructor
    Base2DDatumFormatter(const std::string &_name, const std::string &_desc) :
      BaseValueFormatter<Container, DataObject>(_name, _desc),
      m_current_ptr(nullptr),
      m_known(false) {}

    /// \brief Destructor
    virtual ~Base2DDatumFormatter() {}

    /// \brief Clone
    std::unique_ptr<Base2DDatumFormatter> clone() const {
      return std::unique_ptr<Base2DDatumFormatter>(this->_clone());
    }


    // --- These methods specialize virtual BaseValueFormatter<Container, DataObject> methods ----


    /// \brief Default initialization adds rules for each element
    virtual bool init(const DataObject &_template_obj) const override {
      _prepare(_template_obj);
      return m_known;
    }


    virtual Index num_passes(const DataObject &_template_obj) const override {
      _prepare(_template_obj);
      Index result(0);
      if(!m_known)
        result = 1;
      else if(_index_rules().size() && _index_rules()[0].size()) {
        if(valid_index(_index_rules()[0][0].first))
          result = _index_rules().size();
        else
          result = ContainerTraits<Container>::rows(m_cache);
      }

      //std::cout << "Requesting " << result << " passes.\n";
      return result;

    }

    /// \brief Default col_header uses 'name(index)' for each column
    ///
    /// Ex: "corr(i:j,0)" "corr(i:j,1)" "corr(i:j,5)" "corr(i:j,6)"
    virtual std::vector<std::string> col_header(const DataObject &_template_obj) const override {
      std::vector<std::string> _col;
      init(_template_obj);
      if(!_index_rules().size())
        return _col;

      Index s = max(8 - int(this->name().size()), 0);

      for(Index j = 0; j < _index_rules()[0].size(); ++j) {
        Index a(_index_rules()[0][j].first), b(_index_rules().back()[j].first);
        std::stringstream t_ss;
        t_ss << std::string(s, ' ') << this->name() << '(';

        if(!valid_index(a)) {
          t_ss << ":";
        }
        else if(a == b)
          t_ss << a;
        else
          t_ss << a << ":" << b;
        t_ss << ", " << _index_rules()[0][j].second << ')';
        _col.push_back(t_ss.str());
      }
      return _col;
    }

    /// \brief Default implementation calls _parse_index_expression
    virtual bool parse_args(const std::string &args) override {
      _parse_index_expression(args);
      return true;
    }

    /// \brief Default implementation injects each element
    ///
    /// - sets DataStream::failbit if validation fails
    virtual void inject(const DataObject &_data_obj, DataStream &_stream, Index pass_index = 0) const override {

      // add_rules to print all elements if not set yet
      _prepare(_data_obj);

      if(!m_known)
        _stream << DataStream::failbit;

      Index i = pass_index;
      Index row = 0;
      Index *row_ptr = &row;

      Index rows = _index_rules().size();
      Index cols = 0;
      if(rows && _index_rules()[0].size()) {
        cols = _index_rules()[0].size();
        if(m_known) {
          if(!valid_index(_index_rules()[0][0].first)) {
            rows = ContainerTraits<Container>::rows(m_cache);
            row_ptr = &pass_index;
            i = 0;
          }
        }
      }
      if(i <= _index_rules().size()) {
        for(Index j = 0; j < cols; ++j) {
          row = _index_rules()[i][j].first;
          _stream << (m_known ? Access::at(m_cache, *row_ptr, _index_rules()[i][j].second) : ValueType());
        }
      }
    }


    /// \brief Default implementation prints each element in a column
    ///
    /// - Prints "unknown" if validation fails
    virtual void print(const DataObject &_data_obj, std::ostream &_stream, Index pass_index = 0) const override {

      // add_rules to print all elements if not set yet
      _prepare(_data_obj);

      _stream.flags(std::ios::showpoint | std::ios::fixed | std::ios::right);
      _stream.precision(8);

      Index i = pass_index;
      Index row = 0;
      Index *row_ptr = &row;

      Index rows = _index_rules().size();
      Index cols = 0;
      if(rows && _index_rules()[0].size()) {
        cols = _index_rules()[0].size();
        if(m_known) {
          if(!valid_index(_index_rules()[0][0].first)) {
            rows = ContainerTraits<Container>::rows(m_cache);
            row_ptr = &pass_index;
            i = 0;
          }
        }
      }

      if(i <= _index_rules().size()) {
        for(Index j = 0; j < cols; ++j) {
          row = _index_rules()[i][j].first;
          if(m_known)
            _stream << "  " << Access::at(m_cache, *row_ptr, _index_rules()[i][j].second);
          else
            _stream << "  unknown";
        }
      }
    }

  protected:

    typedef multivector<std::pair<Index, Index> >::X<2> IndexContainer;

    mutable IndexContainer m_2D_index_rules;


    /// Derived DatumFormatters have some optional functionality for parsing index
    /// expressions in order to make it easy to handle ranges such as:
    /// \code
    ///       formatter_name(3,4:8)
    /// \endcode
    /// in which case, DerivedDatumFormatter::parse_args() is called with the string "5:6,4:8"
    /// by dispatching that string to BaseDatumFormatter::_parse_index_expression(),
    /// m_2D_index_rules will be populated with
    ///    {{{5,4},{5,5},{5,6},{5,7},{5,8}},
    ///     {{6,4},{6,5},{6,6},{6,7},{6,8}}
    void _parse_index_expression(const std::string &_expr) {
      //std::cout << "Parsing index expression: " << _expr << "\n";
      auto bounds = index_expression_to_bounds(_expr);
      //std::vector<difference_type> ind_begin(bounds.first.rbegin(), bounds.first.rend());
      //std::vector<difference_type> ind_end(bounds.second.rbegin(), bounds.second.rend());
      if(bounds.first.empty() && bounds.second.empty())
        return;

      if(bounds.first.size() != 2 || bounds.second.size() != 2) {
        throw std::runtime_error("Attempted to initialize 2D DatumFormatter with incompatible index expression: " + _expr);
      }

      Index r = 0;
      for(difference_type i = bounds.first[0]; i < bounds.second[0]; ++i, ++r) {
        for(difference_type j = bounds.first[1]; j < bounds.second[1]; ++j) {
          _add_rule(r, std::make_pair(i, j));
        }
      }

    }

    void _add_rule(Index row, std::pair<Index, Index> const &new_rule) const {
      while(m_2D_index_rules.size() < row + 1)
        m_2D_index_rules.push_back({});

      m_2D_index_rules[row].push_back(new_rule);
    }

    const IndexContainer &_index_rules() const {
      return m_2D_index_rules;
    }

    void _prepare(DataObject const &_data_obj) const {
      if(m_current_ptr != &_data_obj) {
        m_current_ptr = &_data_obj;
        m_known = this->validate(_data_obj);
        if(m_known) {
          m_cache = this->evaluate(_data_obj);
        }
      }

      if(m_known && _index_rules().empty()) {

        //std::cout << "creating idices, cache size: " << ContainerTraits<Container>::rows(m_cache) << ", " << ContainerTraits<Container>::cols(m_cache) << "\n"
        //        << "cache matrix:\n" << m_cache << "\n";
        Index cols = ContainerTraits<Container>::cols(m_cache);
        for(Index j = 0; j < cols; j++) {
          _add_rule(0, std::make_pair(-1, j));
        }
      }

    }

  private:
    mutable DataObject const *m_current_ptr;
    mutable bool m_known;
    mutable Container m_cache;

    /// \brief Clone
    virtual Base2DDatumFormatter *_clone() const override = 0;

  };

  template<typename DataObject>
  MatrixXdAttributeDictionary<DataObject> make_matrixxd_dictionary() {
    return MatrixXdAttributeDictionary<DataObject>();
  }

}

#endif
