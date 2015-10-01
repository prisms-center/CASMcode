#ifndef DATAFORMATTERTOOLS_HH
#define DATAFORMATTERTOOLS_HH
#include "casm/casm_io/DataFormatter.hh"
#include "casm/misc/CASM_TMP.hh"

namespace CASM {
  /* ConstantValueFormatter prints a string value specified at construction.  A header string can also be passed.
   */
  template<typename ValueType, typename DataObject>
  class ConstantValueFormatter: public BaseDatumFormatter<DataObject> {
  public:
    using BaseDatumFormatter<DataObject>::name;

    ConstantValueFormatter(const std::string &_header, const ValueType &_value, bool _print_json = false) :
      BaseDatumFormatter<DataObject>(_header, "Constant value, passed as string, that will be printed on each line"), m_value(_value), m_print_json(_print_json) {}

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

  template<typename Container>
  Index container_size_1D(const Container &cont) {
    return cont.size();
  }

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
                            Sizer sizer = container_size_1D<Container>) :
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

  template<typename Container>
  std::vector<Index> container_size_2D(const Container &cont) {
    std::vector<Index> tsize(2, 0);
    tsize[0] = cont.size();
    if(tsize[0] > 0)
      tsize[1] = cont[0].size();
    return tsize;
  }

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
                            Sizer sizer = container_size_2D<Container>) :
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
