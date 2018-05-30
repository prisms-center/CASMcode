#ifndef CLEXPARAMPACK_HH
#define CLEXPARAMPACK_HH
#include <cstddef>
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

    virtual size_type size(ClexParamKey  const &_key) const = 0;

    virtual std::vector<double> const &read(ClexParamKey  const &_key) const = 0;
    virtual double read(ClexParamKey  const &_key, size_type _ind) const = 0;

    virtual void write(ClexParamKey const &_key, std::vector<double> const &_val) = 0;
    virtual void write(ClexParamKey const &_key, size_type _ind, double val) = 0;

  private:
    //possible implementation:
    //std::vector<std::vector<double> m_data;
  };

  class BasicClexParamKey : public ClexParamPack_impl::BaseKey {
  public:
    typedef ClexParamPack::size_type size_type;

    BasicClexParamKey(std::string const &_name, size_type _ind) :
      ClexParamPack_impl::BaseKey(_name),
      m_index(_ind) {}

    size_type index() const {
      return m_index;
    }
  protected:

    ClexParamPack_impl::BaseKey *_clone() const override {
      return new BasicClexParamKey(*this);
    }

    size_type m_index;
  };

  /// \brief Abstract base class for reading/writing clexulator parameters
  class BasicClexParamPack : public ClexParamPack {
  public:

    BasicClexParamPack(Index Nkey);

    size_type size(ClexParamKey  const &_key) const override {
      return size(*static_cast<BasicClexParamKey const *>(_key.ptr()));
    }

    size_type size(BasicClexParamKey  const &_key) const {
      return m_data[_key.index()].size();
    }

    std::vector<double> const &read(ClexParamKey const &_key) const override {
      return read(*static_cast<BasicClexParamKey const *>(_key.ptr()));
    }

    std::vector<double> const &read(BasicClexParamKey const &_key) const {
      return m_data[_key.index()];
    }

    double read(ClexParamKey const &_key, size_type _ind) const override {
      return read(*static_cast<BasicClexParamKey const *>(_key.ptr()), _ind);
    }

    double read(BasicClexParamKey const &_key, size_type _ind) const {
      return m_data[_key.index()][_ind];
    }

    void write(ClexParamKey const &_key, std::vector<double> const &_val) override {
      write(*static_cast<BasicClexParamKey const *>(_key.ptr()), _val);
    }

    void write(BasicClexParamKey const &_key, std::vector<double> const &_val) {
      m_data[_key.index()] = _val;
    }

    void write(ClexParamKey const &_key, size_type _ind, double _val) override {
      write(*static_cast<BasicClexParamKey const *>(_key.ptr()), _ind, _val);
    }

    void write(BasicClexParamKey const &_key, size_type _ind, double _val) {
      m_data[_key.index()][_ind] = _val;
    }

    void resize(BasicClexParamKey  const &_key, Index _size) {
      m_data[_key.index()].resize(_size);
    }


  private:
    std::vector<std::vector<double> > m_data;
  };
}

#endif

