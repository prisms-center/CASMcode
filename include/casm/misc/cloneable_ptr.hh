#ifndef CASM_cloneable_ptr_HH
#define CASM_cloneable_ptr_HH

#include <memory>
#include <iostream>
#include <utility>

#include "casm/misc/CASM_TMP.hh"

/// \brief Non-std smart pointer classes and functions
namespace notstd {

  /// \brief c++11 does not include 'make_unique'
  template<typename T, typename ...Args>
  std::unique_ptr<T> make_unique(Args &&...args) {
    return std::unique_ptr<T>(new T(std::forward<Args>(args)...));
  }

  class Cloneable {
  public:
    virtual ~Cloneable() {}

    std::unique_ptr<Cloneable> clone() const {
      return std::unique_ptr<Cloneable>(this->_clone());
    }
  private:
    virtual Cloneable *_clone() const = 0;
  };

  template<typename Type>
  class cloneable_ptr;

  /// \brief make a cloneable_ptr<T> via T(Args... args)
  template<typename T, typename ...Args>
  cloneable_ptr<T> make_cloneable(Args &&...args) {
    return cloneable_ptr<T>(new T(std::forward<Args>(args)...));
  }

  /// \brief Base type inherits from std::false_type if T does not have clone method
  template <typename T, typename = void>
  struct has_clone : std::false_type {
  };

  /// \brief Specialized case inherits from std::true_type if T does have clone method
  template <typename T>
struct has_clone<T, CASM::CASM_TMP::void_t<decltype(&T::clone)> > : std::true_type {
  };

  template < typename T, typename std::enable_if < !has_clone<T>::value, void >::type * = nullptr >
  std::unique_ptr<T> clone(const T &obj);

  template<typename T, typename std::enable_if<has_clone<T>::value, void>::type * = nullptr>
  std::unique_ptr<T> clone(const T &obj);


  /// \brief A 'cloneable_ptr' can be used in place of 'unique_ptr'
  ///
  /// If you are creating a class 'MyNewClass' that will have 'std::unique_ptr<Type> m_data'
  /// as a class member, you have to explicitly write a copy constructor.
  ///
  /// If Type has method 'virtual std::unique_ptr<Type> clone() const'
  /// you can use 'notstd::cloneable_ptr<Type> m_data' in place of 'std::unique_ptr<Type> m_data'
  /// as a class member and not have to explicity write copy constructors.
  ///
  template<typename Type>
  class cloneable_ptr {

  public:

    typedef Type element_type;
    typedef Type *pointer;
    typedef Type &reference;

    /// \brief Default constructor
    cloneable_ptr() {}

    /// \brief Construct by taking ownership of ptr
    explicit cloneable_ptr(pointer ptr);

    /// \brief Construct by cloning other
    cloneable_ptr(const cloneable_ptr &other);

    /// \brief Construct by cloning other
    template<typename U>
    cloneable_ptr(const cloneable_ptr<U> &other);

    /// \brief Construct by moving other
    cloneable_ptr(cloneable_ptr &&other);

    /// \brief Construct by moving other
    template<typename U>
    cloneable_ptr(cloneable_ptr<U> &&other);


    /// \brief Construct by cloning other
    template<typename U>
    cloneable_ptr(const std::unique_ptr<U> &other);

    /// \brief Construct by moving
    template<typename U>
    cloneable_ptr(std::unique_ptr<U> &&other);

    ~cloneable_ptr();


    /// \brief Assignment via copy-swap
    cloneable_ptr &operator=(cloneable_ptr other);

    /// \brief Assignment via move
    template<typename U>
    cloneable_ptr &operator=(cloneable_ptr<U> &&other);


    /// \brief Assignment via clone
    template<typename U>
    cloneable_ptr &operator=(const std::unique_ptr<U> &other);

    /// \brief Assignment via move
    template<typename U>
    cloneable_ptr &operator=(std::unique_ptr<U> &&other);


    reference operator*() const;

    pointer operator->() const;

    /// \brief Reset contained unique_ptr
    void reset();

    /// \brief Access contained unique_ptr
    std::unique_ptr<Type> &unique();

    /// \brief const Access contained unique_ptr
    const std::unique_ptr<Type> &unique() const;

    /// \brief Checks whether *this owns an object
    explicit operator bool() const;

  private:

    std::unique_ptr<Type> m_unique;

  };

  template<typename Type>
  void swap(cloneable_ptr<Type> &A, cloneable_ptr<Type> &B);

  template<typename Type>
  bool operator<(const cloneable_ptr<Type> &A, const cloneable_ptr<Type> &B);

  template<typename Type>
  bool operator==(const cloneable_ptr<Type> &A, const cloneable_ptr<Type> &B);

  template<typename Type>
  bool operator!=(const cloneable_ptr<Type> &A, const cloneable_ptr<Type> &B);
}

#include "casm/misc/cloneable_ptr_impl.hh"

#endif
