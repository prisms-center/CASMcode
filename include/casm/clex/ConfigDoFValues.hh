#ifndef CASM_ConfigDoFValues
#define CASM_ConfigDoFValues

#include "casm/basis_set/DoFTraits.hh"
namespace CASM {
  class jsonParser;

  class ConfigDoFValues {
  public:
    ConfigDoFValues() : m_n_basis(0), m_n_vol(0)
    {}
    ConfigDoFValues(DoFType::BasicTraits const &_traits, Index _n_basis, Index _n_vol) :
      m_type(_traits.type_name()),
      m_n_basis(_n_basis),
      m_n_vol(_n_vol) {
    }

    std::string const &type_name() const {
      return m_type;
    }

    Index n_vol() const {
      return m_n_vol;
    }

    Index n_basis() const {
      return m_n_basis;
    }

  private:
    DoFKey m_type;
    Index m_n_basis;
    Index m_n_vol;
  };

  class LocalDiscreteConfigDoFValues : public ConfigDoFValues {
  public:

    typedef Eigen::VectorXi ValueType;
    typedef Eigen::VectorXi &Reference;
    typedef const Eigen::VectorXi &ConstReference;

    typedef typename ValueType::Scalar SiteValueType;
    typedef int &SiteReference;
    typedef const int &ConstSiteReference;

    typedef ValueType SublatValueType;
    typedef typename ValueType::SegmentReturnType SublatReference;
    typedef typename ValueType::ConstSegmentReturnType ConstSublatReference;

    LocalDiscreteConfigDoFValues() {}

    LocalDiscreteConfigDoFValues(DoFType::BasicTraits const &_traits, Index _n_basis, Index _n_vol, Eigen::Ref< const ValueType > const &_vals) :
      ConfigDoFValues(_traits, _n_basis, _n_vol),
      m_vals(_vals) {

    }

    Reference values() {
      return m_vals;
    }

    ConstReference values() const {
      return m_vals;
    }

    SublatReference sublat(Index b) {
      return m_vals.segment((b - 1) * n_vol(), n_vol());
    }

    ConstSublatReference sublat(Index b) const {
      return m_vals.segment((b - 1) * n_vol(), n_vol());
    }

  private:
    ValueType m_vals;
  };

  jsonParser &to_json(LocalDiscreteConfigDoFValues const &_values, jsonParser &_json);
  void from_json(LocalDiscreteConfigDoFValues &_values, jsonParser const &_json);

  class LocalContinuousConfigDoFValues : public ConfigDoFValues {
  public:

    typedef Eigen::MatrixXd ValueType;
    typedef Eigen::MatrixXd &Reference;
    typedef const Eigen::MatrixXd &ConstReference;

    typedef Eigen::VectorXd SiteValueType;
    typedef typename ValueType::ColXpr SiteReference;
    typedef const typename ValueType::ConstColXpr ConstSiteReference;

    typedef ValueType SublatValueType;
    typedef typename Eigen::Block<ValueType> SublatReference;
    typedef const typename Eigen::Block<const ValueType> ConstSublatReference;

    LocalContinuousConfigDoFValues() {}

    LocalContinuousConfigDoFValues(DoFType::BasicTraits const &_traits, Index _n_basis, Index _n_vol,  Eigen::Ref< const ValueType > const &_vals) :
      ConfigDoFValues(_traits, _n_basis, _n_vol),
      m_vals(_vals) {

    }

    Reference values() {
      return m_vals;
    }

    ConstReference values() const {
      return m_vals;
    }

    SiteReference site_value(Index l) {
      return m_vals.col(l);
    }

    ConstSiteReference site_value(Index l) const {
      return m_vals.col(l);
    }

    SublatReference sublat(Index b) {
      return m_vals.block(0, (b - 1) * n_vol(), m_vals.rows(), n_vol());
    }

    ConstSublatReference sublat(Index b) const {
      return m_vals.block(0, (b - 1) * n_vol(), m_vals.rows(), n_vol());
    }

  private:
    ValueType m_vals;
  };

  jsonParser &to_json(LocalContinuousConfigDoFValues const &_values, jsonParser &_json);
  void from_json(LocalContinuousConfigDoFValues &_values, jsonParser const &_json);

  class GlobalContinuousConfigDoFValues : public ConfigDoFValues {
  public:

    typedef Eigen::VectorXd ValueType;
    typedef Eigen::VectorXd &Reference;
    typedef const Eigen::VectorXd &ConstReference;

    typedef typename ValueType::Scalar SiteValueType;
    typedef int &SiteReference;
    typedef const int &ConstSiteReference;

    GlobalContinuousConfigDoFValues() {}

    GlobalContinuousConfigDoFValues(DoFType::BasicTraits const &_traits, Index _n_basis, Index _n_vol, Eigen::Ref< const ValueType > const &_vals) :
      ConfigDoFValues(_traits, _n_basis, _n_vol),
      m_vals(_vals) {

    }

    Reference values() {
      return m_vals;
    }

    ConstReference values() const {
      return m_vals;
    }

  private:
    ValueType m_vals;
  };

  jsonParser &to_json(GlobalContinuousConfigDoFValues const &_values, jsonParser &_json);
  void from_json(GlobalContinuousConfigDoFValues &_values, jsonParser const &_json);


}

#endif
