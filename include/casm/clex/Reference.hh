#ifndef CASM_Reference
#define CASM_Reference

#include <map>
#include <functional>
#include <memory>
#include <string>

#include "casm/external/Eigen/Dense"
#include "casm/misc/cloneable_ptr.hh"
#include "casm/casm_io/DataFormatter.hh"
#include "casm/casm_io/DataFormatterTools.hh"

namespace CASM {

  class Configuration;

  typedef ScalarAttribute<Configuration> Reference;

  /** \defgroup Reference
   *
   *  \ingroup Clex
   *  \brief Specify a reference for a cluster expanded property
   *  @{
   */

  /// \brief Maps all Configurations to the same value
  class ConstantReference : public Reference {

  public:

    static const std::string Name;
    static const std::string Desc;

    explicit ConstantReference(double _value = 0.0) :
      Reference(Name, Desc), m_value(_value) {}

    /// \brief Return the constant value used for all Configuration
    double value() const {
      return m_value;
    }

    /// \brief Returns the same constant value for all Configuration
    double evaluate(const Configuration &config) const override {
      return value();
    }

    std::unique_ptr<ConstantReference> clone() const {
      return notstd::make_unique<ConstantReference>(*this->_clone());
    }

  private:

    /// \brief Clone
    virtual ConstantReference *_clone() const override {
      return new ConstantReference(*this);
    }

    double m_value;

  };



  /// \brief Maps a Configuration to a scalar value via a hyperplane
  ///
  /// - HyperPlaneReferenceBase protects mutators
  ///
  /// A hyperplane reference, R, maps input vector coordinates, x, to output scalar value, y:
  /// - \code y = R.dot(x) \endcode
  ///
  /// A HyperPlaneReference is constructed with a global reference, but may be
  /// specialized to give a different R for a particular Supercell or Configuration.
  class HyperPlaneReferenceBase : public Reference {

  public:

    typedef std::map<std::string, Eigen::VectorXd> SpecializedRef;
    typedef std::function<Eigen::VectorXd(const Configuration &)> InputFunction;

    /// \brief Constructor
    ///
    /// \param _global_ref An Eigen::VectorXd giving the intercepts of the
    ///                    hyperplane used for the global reference
    /// \param _input The 'input' function maps a Configuration to an Eigen::VectorXd,
    ///        typically the composition of a Configuration
    /// \param _supercell_ref A map of scelname to Eigen::VectorXd specializing the
    ///                       the reference value by Supercell
    /// \param _config_ref A map of configname to Eigen::VectorXd specializing the
    ///                       the reference value by Configuration
    ///
    /// A hyperplane reference, R, maps input vector coordinates, x, to output scalar value, y:
    /// - \code y = R.dot(x) \endcode
    ///
    /// The global reference, '_global_ref', is required, but may be specialized
    /// to give a different R for a particular Supercell or Configuration via
    /// optional '_supercell_ref' and '_config_ref'.
    ///
    HyperPlaneReferenceBase(std::string _name,
                            std::string _desc,
                            const Eigen::VectorXd &_global_ref,
                            InputFunction _input,
                            SpecializedRef _supercell_ref = SpecializedRef(),
                            SpecializedRef _config_ref = SpecializedRef()) :
      Reference(_name, _desc),
      m_input(_input),
      m_config_ref(_config_ref),
      m_supercell_ref(_supercell_ref),
      m_global_ref(_global_ref) {}

    virtual ~HyperPlaneReferenceBase() {}

    /// \brief Clone
    std::unique_ptr<HyperPlaneReferenceBase> clone() const {
      return notstd::make_unique<HyperPlaneReferenceBase>(*this->_clone());
    }

    /// \brief const Access the global reference
    ///
    const Eigen::VectorXd &global() const {
      return m_global_ref;
    }


    /// \brief const Access a map of scelname to reference for Supercell
    ///        specialized references
    ///
    const std::map<std::string, Eigen::VectorXd> &supercell() const {
      return m_supercell_ref;
    }

    /// \brief const Access a map of configname to reference for Configuration
    ///        specialized references
    ///
    const std::map<std::string, Eigen::VectorXd> &config() const {
      return m_config_ref;
    }

    /// \brief Return the 'input' function that maps a Configuration to coordinates
    std::function<Eigen::VectorXd(const Configuration &)> input() const {
      return m_input;
    }

    /// \brief Return the 'input' coordinates that a Configuration is mapped to
    Eigen::VectorXd input(const Configuration &config) const {
      return m_input(config);
    }

    /// \brief Return the reference hyperplane used for a particular configuration
    ///
    /// Returns the Configuration specific hyperplane if it exists, else the
    /// Supercell specific hyperplane if it exists, else the global hyperplane.
    Eigen::VectorXd hyperplane(const Configuration &config) const;


    /// \brief Return the reference for a particular configuration
    ///
    /// \returns
    /// \code
    /// hyperplane(config).dot(input(config));
    /// \endcode
    double evaluate(const Configuration &config) const override {

      return hyperplane(config).dot(input(config));
    }

  protected:

    /// \brief Access the global reference
    ///
    Eigen::VectorXd &global() {
      return m_global_ref;
    }

    /// \brief Access a map of scelname to reference for Supercell
    ///        specialized references
    ///
    std::map<std::string, Eigen::VectorXd> &supercell() {
      return m_supercell_ref;
    }

    /// \brief Access a map of configname to reference for Configuration
    ///        specialized references
    ///
    std::map<std::string, Eigen::VectorXd> &config() {
      return m_config_ref;
    }


  private:

    /// \brief Clone
    virtual HyperPlaneReferenceBase *_clone() const override {
      return new HyperPlaneReferenceBase(*this);
    }


    InputFunction m_input;

    std::map<std::string, Eigen::VectorXd> m_config_ref;

    std::map<std::string, Eigen::VectorXd> m_supercell_ref;

    Eigen::VectorXd m_global_ref;

  };


  /// \brief Maps a Configuration to a scalar value via a hyperplane
  ///
  /// A hyperplane reference, R, maps input vector coordinates, x, to output scalar value, y:
  /// - \code y = R.dot(x) \endcode
  ///
  /// A HyperPlaneReference is constructed with a global reference, but may be
  /// specialized to give a different R for a particular Supercell or Configuration.
  class HyperPlaneReference : public HyperPlaneReferenceBase {

  public:

    static const std::string Name;
    static const std::string Desc;

    /// \brief Constructor
    ///
    /// \param _global_ref An Eigen::VectorXd giving the intercepts of the
    ///                    hyperplane used for the global reference
    /// \param _input The 'input' function maps a Configuration to an Eigen::VectorXd,
    ///        typically the composition of a Configuration
    /// \param _supercell_ref A map of scelname to Eigen::VectorXd specializing the
    ///                       the reference value by Supercell
    /// \param _config_ref A map of configname to Eigen::VectorXd specializing the
    ///                       the reference value by Configuration
    ///
    /// A hyperplane reference, R, maps input vector coordinates, x, to output scalar value, y:
    /// - \code y = R.dot(x) \endcode
    ///
    /// The global reference, '_global_ref', is required, but may be specialized
    /// to give a different R for a particular Supercell or Configuration via
    /// optional '_supercell_ref' and '_config_ref'.
    ///
    HyperPlaneReference(const Eigen::VectorXd &_global_ref,
                        InputFunction _input,
                        SpecializedRef _supercell_ref = SpecializedRef(),
                        SpecializedRef _config_ref = SpecializedRef()) :
      HyperPlaneReferenceBase(Name, Desc, _global_ref, _input, _config_ref, _supercell_ref) {}


    /// \brief Clone
    std::unique_ptr<HyperPlaneReference> clone() const {
      return notstd::make_unique<HyperPlaneReference>(*this->_clone());
    }


    // --- Prevent hiding of const accessors ---

    using HyperPlaneReferenceBase::global;
    using HyperPlaneReferenceBase::supercell;
    using HyperPlaneReferenceBase::config;


    // --- Make protected inherited members public ---

    /// \brief Access the global reference
    ///
    Eigen::VectorXd &global() {
      return HyperPlaneReferenceBase::global();
    }

    /// \brief Access a map of scelname to reference for Supercell
    ///        specialized references
    ///
    std::map<std::string, Eigen::VectorXd> &supercell() {
      return HyperPlaneReferenceBase::supercell();
    }

    /// \brief Access a map of configname to reference for Configuration
    ///        specialized references
    ///
    std::map<std::string, Eigen::VectorXd> &config() {
      return HyperPlaneReferenceBase::config();
    }

  private:

    /// \brief Clone
    HyperPlaneReference *_clone() const override {
      return new HyperPlaneReference(*this);
    }

  };

  /** @} */

}

#endif
