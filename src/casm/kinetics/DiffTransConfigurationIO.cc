#include "casm/kinetics/DiffTransConfigurationIO.hh"

namespace CASM {

  template class BaseDatumFormatter<Kinetics::DiffTransConfiguration>;
  template class DataFormatterOperator<bool, std::string, Kinetics::DiffTransConfiguration>;
  template class DataFormatterOperator<bool, bool, Kinetics::DiffTransConfiguration>;
  template class DataFormatterOperator<bool, double, Kinetics::DiffTransConfiguration>;
  template class DataFormatterOperator<double, double, Kinetics::DiffTransConfiguration>;
  template class DataFormatterOperator<Index, double, Kinetics::DiffTransConfiguration>;
  template class DataFormatter<Kinetics::DiffTransConfiguration>;
  template bool DataFormatter<Kinetics::DiffTransConfiguration>::evaluate_as_scalar<bool>(Kinetics::DiffTransConfiguration const &) const;
  template double DataFormatter<Kinetics::DiffTransConfiguration>::evaluate_as_scalar<double>(Kinetics::DiffTransConfiguration const &) const;
  template class DataFormatterDictionary<Kinetics::DiffTransConfiguration>;

  template<>
  StringAttributeDictionary<Kinetics::DiffTransConfiguration> make_string_dictionary<Kinetics::DiffTransConfiguration>() {

    using namespace Kinetics::DiffTransConfigIO;
    StringAttributeDictionary<Kinetics::DiffTransConfiguration> dict;
    return dict;
  }

  template<>
  BooleanAttributeDictionary<Kinetics::DiffTransConfiguration> make_boolean_dictionary<Kinetics::DiffTransConfiguration>() {

    using namespace Kinetics::DiffTransConfigIO;
    BooleanAttributeDictionary<Kinetics::DiffTransConfiguration> dict;
    return dict;
  }

  template<>
  IntegerAttributeDictionary<Kinetics::DiffTransConfiguration> make_integer_dictionary<Kinetics::DiffTransConfiguration>() {

    using namespace Kinetics::DiffTransConfigIO;
    IntegerAttributeDictionary<Kinetics::DiffTransConfiguration> dict;
    return dict;
  }

  template<>
  ScalarAttributeDictionary<Kinetics::DiffTransConfiguration> make_scalar_dictionary<Kinetics::DiffTransConfiguration>() {

    using namespace Kinetics::DiffTransConfigIO;
    ScalarAttributeDictionary<Kinetics::DiffTransConfiguration> dict;
    return dict;
  }

  template<>
  VectorXiAttributeDictionary<Kinetics::DiffTransConfiguration> make_vectorxi_dictionary<Kinetics::DiffTransConfiguration>() {
    using namespace Kinetics::DiffTransConfigIO;
    VectorXiAttributeDictionary<Kinetics::DiffTransConfiguration> dict;
    return dict;
  }

  template<>
  VectorXdAttributeDictionary<Kinetics::DiffTransConfiguration> make_vectorxd_dictionary<Kinetics::DiffTransConfiguration>() {

    using namespace Kinetics::DiffTransConfigIO;
    VectorXdAttributeDictionary<Kinetics::DiffTransConfiguration> dict;
    return dict;
  }

}
