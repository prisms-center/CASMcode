#ifndef CLEXPARAMPACK_HH
#define CLEXPARAMPACK_HH
#include <cstddef>
#include <map>
#include <vector>
#include "casm/misc/cloneable_ptr.hh"

namespace CASM {

  namespace ClexParamPack_impl {
    /// \brief BaseKey class that hides implementation-specific access attributes
    class BaseKey {
    public:
      BaseKey(std::string const &_name) : m_name(_name) {}

      virtual ~BaseKey() {};

      std::string const &name() const {
        return m_name;
      }

      /// \brief Clone the ClexParamKey
      std::unique_ptr<BaseKey> clone() const {
        return std::unique_ptr<BaseKey>(_clone());
      }


    protected:

      /// \brief Clone the ClexParamKey
      virtual BaseKey *_clone() const = 0;

    private:
      std::string m_name;
    };
  }

  class ClexParamPack;

  namespace Clexulator_impl {
    class Base;
  }

  /// \brief Key for indexing clexulator parameters
  class ClexParamKey {
  public:
    ClexParamKey() {
    }

    ClexParamKey(ClexParamPack_impl::BaseKey const &_key) :
      m_key_ptr(_key.clone()) {
    }

    ~ClexParamKey() {

    }

    //friend class Clexulator_impl::Base;
    std::string const &name() const {
      return m_key_ptr->name();
    }

    ClexParamPack_impl::BaseKey const *ptr() const {
      return m_key_ptr.unique().get();
    }
  private:
    /// \brief  ptr to BaseKey class that hides implementation-specific access attributes
    notstd::cloneable_ptr<ClexParamPack_impl::BaseKey> m_key_ptr;
  };


  /// \brief Abstract base class for reading/writing clexulator parameters
  class ClexParamPack {
  public:

    typedef unsigned int size_type;

    virtual ~ClexParamPack() {
    }

    std::map<std::string, ClexParamKey> const &keys() const {
      return m_keys;
    }

    ClexParamKey const &key(std::string const &_name) const {
      auto it = keys().find(_name);
      if(it == keys().end())
        throw std::runtime_error("In ClexParamPack::key(), ClexParamPack does not contain parameters corresponding to name " + _name + ".");
      return it->second;
    }

    virtual size_type size(ClexParamKey  const &_key) const = 0;

    virtual size_type dim(ClexParamKey const &_key) const = 0;

    virtual std::string eval_mode(ClexParamKey const &_key) const = 0;

    virtual Eigen::MatrixXd const &read(ClexParamKey  const &_key) const = 0;

    virtual double const &read(ClexParamKey  const &_key, size_type _i) const = 0;

    virtual double const &read(ClexParamKey  const &_key, size_type _i, size_type _j) const = 0;

    virtual void set_eval_mode(ClexParamKey const &_key, std::string const &_mode) = 0;

    virtual void write(ClexParamKey const &_key, Eigen::Ref<const Eigen::MatrixXd> const &_val) = 0;
    virtual void write(ClexParamKey const &_key, size_type _i, double val) = 0;
    virtual void write(ClexParamKey const &_key, size_type _i, size_type _j, double val) = 0;

  protected:
    std::map<std::string, ClexParamKey> m_keys;
  private:
    //possible implementation:
    //std::vector<std::vector<double> m_data;
  };

}

#endif

