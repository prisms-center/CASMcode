#ifndef CLEXPARAMPACK_HH
#define CLEXPARAMPACK_HH
#include <cstddef>

namespace CASM {

  namespace ClexParamPack_impl {
    /// \brief BaseKey class that hides implementation-specific access attributes
    class BaseKey {
    public:
      /// \brief Clone the ClexParamKey
      std::unique_ptr<BaseKey> clone() const {
        return std::unique_ptr<BaseKey>(_clone());
      }

    private:

      /// \brief Clone the ClexParamKey
      virtual BaseKey *_clone() const = 0;

    }
  }

  class ClexParamPack;

  /// \brief Key for indexing clexulator parameters
  class ClexParamKey {
    friend class Clexulator_impl::Base;
  private:
    /// \brief  ptr to BaseKey class that hides implementation-specific access attributes
    std::cloneable_ptr<ClexParamPack_impl::BaseKey> m_key_ptr;
  };


  /// \brief Abstract base class for reading/writing clexulator parameters
  class ClexParamPack {

  public:

    typedef unsigned int size_type;

    size_type size(ClexParamKey  const &_key, size_type _ind) const;

    double read(ClexParamKey  const &_key, size_type _ind) const;
    double read(ClexParamKey  const &_key, size_type _ind) const;


    void write(ClexParamKey const &_key, std::vector<double> const &_val);
    void write(ClexParamKey const &_key, _ind, double val);
  };



}

#endif

