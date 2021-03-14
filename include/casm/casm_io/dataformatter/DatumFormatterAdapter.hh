#ifndef CASM_DatumFormatterAdapter
#define CASM_DatumFormatterAdapter

#include "casm/casm_io/dataformatter/DataFormatter.hh"

namespace CASM {

namespace adapter {

template <typename ToType, typename FromType>
struct Adapter;

}  // namespace adapter

/// Adapt DatumFormatter to reuse with different types
///
template <typename ExternalType, typename InternalType>
class DatumFormatterAdapter : public BaseDatumFormatter<ExternalType> {
 public:
  using DataObject = ExternalType;
  using difference_type =
      typename BaseDatumFormatter<ExternalType>::difference_type;
  using DictType = typename BaseDatumFormatter<ExternalType>::DictType;

  /// Construct with apapter function and base formatter. Use different name and
  /// description than base formatter.
  DatumFormatterAdapter(const std::string &_name, const std::string &_desc,
                        adapter::Adapter<InternalType, ExternalType> _adapt,
                        BaseDatumFormatter<InternalType> const &_base_formatter)
      : BaseDatumFormatter<ExternalType>(_name, _desc),
        m_adapt(_adapt),
        m_base_ptr(notstd::clone(_base_formatter)) {}

  /// Construct with apapter function and base formatter. Use same name and
  /// description as base formatter.
  DatumFormatterAdapter(adapter::Adapter<InternalType, ExternalType> _adapt,
                        BaseDatumFormatter<InternalType> const &_base_formatter)
      : BaseDatumFormatter<ExternalType>(_base_formatter.name(),
                                         _base_formatter.description()),
        m_adapt(_adapt),
        m_base_ptr(notstd::clone(_base_formatter)) {}

  /// Allow polymorphic deletion
  virtual ~DatumFormatterAdapter(){};

  DatumFormatterClass type() const override { return m_base_ptr->type(); }

  /// \brief const Access the dictionary containing this formatter, set during
  /// DictType::lookup
  const DictType &home() const { return this->home(); }

  /// \brief Set the dictionary containing this formatter, set during
  /// DictType::lookup
  void set_home(const DictType &home) const { this->set_home(home); }

  std::unique_ptr<DatumFormatterAdapter<ExternalType, InternalType>> clone()
      const {
    return std::unique_ptr<DatumFormatterAdapter<ExternalType, InternalType>>(
        this->_clone());
  }

  bool init(ExternalType const &_template_obj) const override {
    return m_base_ptr->init(m_adapt(_template_obj));
  }

  bool validate(ExternalType const &_data_obj) const override {
    return m_base_ptr->validate(m_adapt(_data_obj));
  }

  std::vector<std::string> col_header(
      ExternalType const &_template_obj) const override {
    return m_base_ptr->col_header(m_adapt(_template_obj));
  }

  std::string long_header(ExternalType const &_template_obj) const override {
    return m_base_ptr->long_header(m_adapt(_template_obj));
  }

  std::string short_header(ExternalType const &_template_obj) const override {
    return m_base_ptr->short_header(m_adapt(_template_obj));
  }

  Index num_passes(ExternalType const &_data_obj) const override {
    return m_base_ptr->num_passes(m_adapt(_data_obj));
  }

  /// Print formatted data from _data_obj to _stream, while specifying which
  /// output pass is requested If implementation does not depend on pass_index,
  /// it may safely be ignored
  void print(ExternalType const &_data_obj, std::ostream &_stream,
             Index pass_index = 0) const override {
    m_base_ptr->print(m_adapt(_data_obj), _stream, pass_index);
  }

  void inject(ExternalType const &_data_obj, DataStream &_stream,
              Index pass_index = 0) const override {
    m_base_ptr->inject(m_adapt(_data_obj), _stream, pass_index);
  }

  jsonParser &to_json(ExternalType const &_data_obj,
                      jsonParser &json) const override {
    return m_base_ptr->to_json(m_adapt(_data_obj), json);
  }

  bool parse_args(const std::string &args) override {
    return m_base_ptr->parse_args(args);
  }

 private:
  /// \brief Make an exact copy of the formatter (including any initialized
  /// members)
  ///
  BaseDatumFormatter<ExternalType> *_clone() const override {
    return new DatumFormatterAdapter(*this);
  }

  adapter::Adapter<InternalType, ExternalType> m_adapt;
  notstd::cloneable_ptr<BaseDatumFormatter<InternalType>> m_base_ptr;
};

template <typename ExternalType, typename InternalType>
DatumFormatterAdapter<ExternalType, InternalType> make_datum_formatter_adapter(
    BaseDatumFormatter<InternalType> const &_base_formatter) {
  return DatumFormatterAdapter<ExternalType, InternalType>{
      adapter::Adapter<InternalType, ExternalType>(), _base_formatter};
}

}  // namespace CASM

#endif
