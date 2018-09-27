#ifndef CASM_ConfigDoFValues
#define CASM_ConfigDoFValues

#include "casm/basis_set/DoFTraits.hh"
namespace CASM {
  class jsonParser;

  class ConfigDoFValues {
  public:
    ConfigDoFValues(DoFType::BasicTraits const &_traits);

    std::string const &type_name() const {
      return m_type;
    }


  private:
    std::string m_type;
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
    typedef typename ValueType::SegmentReturnType ConstSublatReference;


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

    typedef typename ValueType::Scalar SiteValueType;
    typedef int &SiteReference;
    typedef const int &ConstSiteReference;

    typedef ValueType SublatValueType;
    typedef typename ValueType::SegmentReturnType SublatReference;
    typedef typename ValueType::SegmentReturnType ConstSublatReference;

    LocalContinuousConfigDoFValues(DoFType::BasicTraits const &_traits, Eigen::Ref< const ValueType > const &_vals) :
      ConfigDoFValues(_traits),
      m_vals(_vals) {

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

    typedef ValueType SublatValueType;
    typedef typename ValueType::SegmentReturnType SublatReference;
    typedef typename ValueType::SegmentReturnType ConstSublatReference;

    GlobalContinuousConfigDoFValues(DoFType::BasicTraits const &_traits, Eigen::Ref< const ValueType > const &_vals) :
      ConfigDoFValues(_traits),
      m_vals(_vals) {

    }

  private:
    ValueType m_vals;
  };

  jsonParser &to_json(GlobalContinuousConfigDoFValues const &_values, jsonParser &_json);
  void from_json(GlobalContinuousConfigDoFValues &_values, jsonParser const &_json);


}

#endif
