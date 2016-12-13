#ifndef DATAFORMATTERTOOLS_HH
#define DATAFORMATTERTOOLS_HH
#include <numeric>
#include <iterator>
#include "casm/external/boost.hh"
#include "casm/casm_io/DataFormatter.hh"
#include "casm/casm_io/json_io/container.hh"
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

    void init(const DataObject &_template_obj) const  override {
      m_formatter->init(_template_obj);
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
    virtual void init(const DataObject &_template_obj) const override {
      if(_index_rules().size())
        return;

      Index size = ContainerTraits<Container>::size(this->evaluate(_template_obj));
      for(Index i = 0; i < size; i++) {
        _add_rule(std::vector<Index>({i}));
      }
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
          _stream << "    " << Access::at(val, (*it)[0]);
        else
          _stream << "    unknown";

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
      Base1DDatumFormatter<DataObject, Container>(_name, _desc),
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
    Generic1DDatumFormatter *_clone() const {
      return new Generic1DDatumFormatter(*this);
    }

    Evaluator m_evaluate;
    Validator m_validate;
  };

  /// \brief Template alias for BaseValueFormatter returning std::string
  ///
  /// \ingroup DataFormatterTypes
  ///
  template<typename DataObject>
  using StringAttribute = BaseValueFormatter<std::string, DataObject>;

  /// \brief Template to be specialized for constructing dictionaries for particular DataObject
  ///
  /// \ingroup DataFormatter
  ///
  template<typename DataObject>
  using StringAttributeDictionary = DataFormatterDictionary<DataObject, StringAttribute<DataObject> >;

  /// \brief Template to be specialized for constructing dictionaries for particular DataObject
  ///
  /// \ingroup DataFormatter
  ///
  template<typename DataObject>
  StringAttributeDictionary<DataObject> make_string_dictionary();


  /// \brief Template alias for BaseValueFormatter returning bool
  ///
  /// \ingroup DataFormatterTypes
  ///
  template<typename DataObject>
  using BooleanAttribute = BaseValueFormatter<bool, DataObject>;

  /// \brief Template to be specialized for constructing dictionaries for particular DataObject
  ///
  /// \ingroup DataFormatter
  ///
  template<typename DataObject>
  using BooleanAttributeDictionary = DataFormatterDictionary<DataObject, BooleanAttribute<DataObject> >;

  /// \brief Template to be specialized for constructing dictionaries for particular DataObject
  ///
  /// \ingroup DataFormatter
  ///
  template<typename DataObject>
  BooleanAttributeDictionary<DataObject> make_boolean_dictionary();


  /// \brief Template alias for BaseValueFormatter returning Index
  ///
  /// \ingroup DataFormatterTypes
  ///
  template<typename DataObject>
  using IntegerAttribute = BaseValueFormatter<Index, DataObject>;

  /// \brief Template to be specialized for constructing dictionaries for particular DataObject
  ///
  /// \ingroup DataFormatter
  ///
  template<typename DataObject>
  using IntegerAttributeDictionary = DataFormatterDictionary<DataObject, IntegerAttribute<DataObject> >;

  /// \brief Template to be specialized for constructing dictionaries for particular DataObject
  ///
  /// \ingroup DataFormatter
  ///
  template<typename DataObject>
  IntegerAttributeDictionary<DataObject> make_integer_dictionary();


  /// \brief Template alias for BaseValueFormatter returning double
  ///
  /// \ingroup DataFormatterTypes
  ///
  template<typename DataObject>
  using ScalarAttribute = BaseValueFormatter<double, DataObject>;

  /// \brief Template to be specialized for constructing dictionaries for particular DataObject
  ///
  /// \ingroup DataFormatter
  ///
  template<typename DataObject>
  using ScalarAttributeDictionary = DataFormatterDictionary<DataObject, ScalarAttribute<DataObject> >;

  /// \brief Template to be specialized for constructing dictionaries for particular DataObject
  ///
  /// \ingroup DataFormatter
  ///
  template<typename DataObject>
  ScalarAttributeDictionary<DataObject> make_scalar_dictionary();


  /// \brief Template alias for BaseValueFormatter returning Eigen::VectorXd
  ///
  /// \ingroup DataFormatterTypes
  ///
  template<typename DataObject>
  using VectorXdAttribute = Base1DDatumFormatter<Eigen::VectorXd, DataObject>;

  /// \brief Template to be specialized for constructing dictionaries for particular DataObject
  ///
  /// \ingroup DataFormatter
  ///
  template<typename DataObject>
  using VectorXdAttributeDictionary = DataFormatterDictionary<DataObject, VectorXdAttribute<DataObject> >;

  /// \brief Template to be specialized for constructing dictionaries for particular DataObject
  ///
  /// \ingroup DataFormatter
  ///
  template<typename DataObject>
  VectorXdAttributeDictionary<DataObject> make_vectorxd_dictionary();


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
      make_vectorxd_dictionary<DataObject>()
    );

    return dict;
  }

  /*
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    template<typename Container>
    std::vector<Index> container_size_2D(const Container &cont) {
      std::vector<Index> tsize(2, 0);
      tsize[0] = cont.size();
      if(tsize[0] > 0)
        tsize[1] = cont[0].size();
      return tsize;
    }

    /// \brief A DatumFormatter that returns a value of specified 2d container
    ///
    /// Functionality is specified via template parameters at compile time
    ///
    /// \ingroup DataFormatter
    ///
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
                              Sizer sizer = container_size_2D<_Container>) :
        BaseDatumFormatter<DataObject>(_init_name, _desc), m_evaluate(evaluator), m_validate(validator), m_size(sizer) {}

      BaseDatumFormatter<DataObject> *clone() const override {
        return new Generic2DDatumFormatter(*this);
      }

      void init(const DataObject &_template_obj) const override;

      std::string col_header(const DataObject &_template_obj) const override;
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

  */

}

#include "casm/casm_io/DataFormatterTools_impl.hh"
#endif
