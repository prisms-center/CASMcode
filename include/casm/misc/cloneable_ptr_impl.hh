#ifndef CASM_cloneable_ptr_impl
#define CASM_cloneable_ptr_impl

#include <memory>
#include <iostream>
#include <utility>

#include "casm/misc/cloneable_ptr.hh"

/// \brief Non-std smart pointer classes and functions
namespace notstd {

  template < typename T, typename std::enable_if < !has_clone<T>::value, void >::type * >
  std::unique_ptr<T> clone(const T &obj) {
    return std::unique_ptr<T>(new T(obj));
  }

  template<typename T, typename std::enable_if<has_clone<T>::value, void>::type * >
  std::unique_ptr<T> clone(const T &obj) {
    return obj.clone();
  }

  /// \brief Construct by taking ownership of ptr
  template<typename Type>
  cloneable_ptr<Type>::cloneable_ptr(pointer ptr) :
    m_unique(ptr) {}

  /// \brief Construct by cloning other
  template<typename Type>
  cloneable_ptr<Type>::cloneable_ptr(const cloneable_ptr &other) :
    m_unique() {
    if(other) {
      m_unique = std::move(clone(*other));
    }
  }

  /// \brief Construct by cloning other
  template<typename Type>
  template<typename U>
  cloneable_ptr<Type>::cloneable_ptr(const cloneable_ptr<U> &other) :
    m_unique() {
    if(other) {
      m_unique = std::move(clone(*other));
    }
  }

  /// \brief Construct by moving other
  template<typename Type>
  cloneable_ptr<Type>::cloneable_ptr(cloneable_ptr &&other) :
    m_unique(std::move(other.unique())) {}

  /// \brief Construct by moving other
  template<typename Type>
  template<typename U>
  cloneable_ptr<Type>::cloneable_ptr(cloneable_ptr<U> &&other) :
    m_unique(std::move(other.unique())) {}


  /// \brief Construct by cloning other
  template<typename Type>
  template<typename U>
  cloneable_ptr<Type>::cloneable_ptr(const std::unique_ptr<U> &other) :
    m_unique() {
    if(other) {
      m_unique = std::move(clone(*other));
    }
  }

  /// \brief Construct by moving
  template<typename Type>
  template<typename U>
  cloneable_ptr<Type>::cloneable_ptr(std::unique_ptr<U> &&other) :
    m_unique(std::move(other)) {}

  template<typename Type>
  cloneable_ptr<Type>::~cloneable_ptr() {}

  /// \brief Assignment via copy-swap
  template<typename Type>
  cloneable_ptr<Type> &cloneable_ptr<Type>::operator=(cloneable_ptr other) {
    swap(*this, other);
    return *this;
  }

  /// \brief Assignment via move
  template<typename Type>
  template<typename U>
  cloneable_ptr<Type> &cloneable_ptr<Type>::operator=(cloneable_ptr<U> &&other) {
    unique() = std::move(other.unique());
    return *this;
  }


  /// \brief Assignment via clone
  template<typename Type>
  template<typename U>
  cloneable_ptr<Type> &cloneable_ptr<Type>::operator=(const std::unique_ptr<U> &other) {
    if(other) {
      unique() = clone(*other);
    }
    else {
      unique().reset();
    }
    return *this;
  }

  /// \brief Assignment via move
  template<typename Type>
  template<typename U>
  cloneable_ptr<Type> &cloneable_ptr<Type>::operator=(std::unique_ptr<U> &&other) {
    unique() = std::move(other);
    return *this;
  }


  template<typename Type>
  typename cloneable_ptr<Type>::reference cloneable_ptr<Type>::operator*() const {
    return *m_unique;
  }

  template<typename Type>
  typename cloneable_ptr<Type>::pointer cloneable_ptr<Type>::operator->() const {
    return &(this->operator*());
  }

  /// \brief Reset contained unique_ptr
  template<typename Type>
  void cloneable_ptr<Type>::reset() {
    return m_unique.reset();
  }

  /// \brief Access contained unique_ptr
  template<typename Type>
  std::unique_ptr<Type> &cloneable_ptr<Type>::unique() {
    return m_unique;
  }

  /// \brief const Access contained unique_ptr
  template<typename Type>
  const std::unique_ptr<Type> &cloneable_ptr<Type>::unique() const {
    return m_unique;
  }

  /// \brief Checks whether *this owns an object
  template<typename Type>
  cloneable_ptr<Type>::operator bool() const {
    return static_cast<bool>(m_unique);
  }

  template<typename Type>
  void swap(cloneable_ptr<Type> &A, cloneable_ptr<Type> &B) {
    A.unique().swap(B.unique());
  }

  template<typename Type>
  bool operator<(const cloneable_ptr<Type> &A, const cloneable_ptr<Type> &B) {
    return A.unique() < B.unique();
  }

  template<typename Type>
  bool operator==(const cloneable_ptr<Type> &A, const cloneable_ptr<Type> &B) {
    return A.unique() == B.unique();
  }

  template<typename Type>
  bool operator!=(const cloneable_ptr<Type> &A, const cloneable_ptr<Type> &B) {
    return A.unique() != B.unique();
  }
}

#endif
