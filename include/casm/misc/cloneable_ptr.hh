#ifndef CASM_cloneable_ptr_HH
#define CASM_cloneable_ptr_HH

#include <memory>
#include <iostream>
#include <utility>

#include "casm/misc/type_traits.hh"

/// Include in a class/struct definition that inherits from Cloneable
/// to make it also cloneable
#define ABSTRACT_CLONEABLE(T) \
public: \
  virtual ~T() {}\
  \
  std::unique_ptr<T> clone() const {\
    return std::unique_ptr<T>(this->_clone());\
  }\
  std::unique_ptr<T> move() {\
    return std::unique_ptr<T>(this->_move());\
  }\
private:\
  virtual T *_clone() const = 0;\
  virtual T *_move() = 0;\
public:\

/// Include in a class/struct definition that inherits from Cloneable
/// to make it also cloneable
#define ABSTRACT_CLONEABLE_DERIVED(T) \
public: \
  virtual ~T() {}\
  \
  std::unique_ptr<T> clone() const {\
    return std::unique_ptr<T>(this->_clone());\
  }\
  std::unique_ptr<T> move() {\
    return std::unique_ptr<T>(this->_move());\
  }\
private:\
  virtual T *_clone() const override = 0;\
  virtual T *_move() override = 0;\
public:\

/// Include in a class/struct definition that inherits from Cloneable
/// to make it also cloneable
#define CLONEABLE(T) \
public: \
  virtual ~T() {}\
  \
  std::unique_ptr<T> clone() const {\
    return std::unique_ptr<T>(this->_clone());\
  }\
  std::unique_ptr<T> move() {\
    return std::unique_ptr<T>(this->_move());\
  }\
private:\
  virtual T *_clone() const override {\
    return new T(*this);\
  }\
  virtual T *_move() override {\
    return new T(std::move(*this));\
  }\
public:\

/// Include in a class/struct definition that inherits from Cloneable
/// to make it also cloneable, but do not include destructor defintion
#define CLONEABLE_NEEDS_DESTRUCTOR_DEF(T) \
public: \
  virtual ~T();\
  \
  std::unique_ptr<T> clone() const {\
    return std::unique_ptr<T>(this->_clone());\
  }\
  std::unique_ptr<T> move() {\
    return std::unique_ptr<T>(this->_move());\
  }\
private:\
  virtual T *_clone() const override {\
    return new T(*this);\
  }\
  virtual T *_move() override {\
    return new T(std::move(*this));\
  }\
public:\


/// \brief Non-std smart pointer classes and functions
namespace notstd {

  /// \brief c++11 does not include 'make_unique'
  template<typename T, typename ...Args>
  std::unique_ptr<T> make_unique(Args &&...args) {
    return std::unique_ptr<T>(new T(std::forward<Args>(args)...));
  }

  /// \brief Base class for cloning
  class Cloneable {
  public:
    virtual ~Cloneable() {}

    std::unique_ptr<Cloneable> clone() const {
      return std::unique_ptr<Cloneable>(this->_clone());
    }
    std::unique_ptr<Cloneable> move() {
      return std::unique_ptr<Cloneable>(this->_move());
    }
  private:
    virtual Cloneable *_clone() const = 0;
    virtual Cloneable *_move() = 0;
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
struct has_clone<T, void_t<decltype(&T::clone)> > : std::true_type {
  };

  template < typename T, typename std::enable_if < !has_clone<T>::value, void >::type * = nullptr >
  std::unique_ptr<T> clone(const T &obj);

  template<typename T, typename std::enable_if<has_clone<T>::value, void>::type * = nullptr>
  std::unique_ptr<T> clone(const T &obj);


  /// \brief Base type inherits from std::false_type if T does not have move method
  template <typename T, typename = void>
  struct has_move : std::false_type {
  };

  /// \brief Specialized case inherits from std::true_type if T does have move method
  template <typename T>
struct has_move<T, void_t<decltype(&T::move)> > : std::true_type {
  };

  /// \brief Construct std::unique_ptr<T> from rvalue reference
  template < typename T, typename std::enable_if < !has_move<T>::value, void >::type * = nullptr >
  std::unique_ptr<T> clone_move(T && obj);

  /// \brief Construct std::unique_ptr<T> from rvalue reference
  template<typename T, typename std::enable_if<has_move<T>::value, void>::type * = nullptr>
  std::unique_ptr<T> clone_move(T && obj);


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

  template<typename Type>
  bool operator==(std::nullptr_t, const cloneable_ptr<Type> &B);

  template<typename Type>
  bool operator!=(std::nullptr_t, const cloneable_ptr<Type> &B);

  template<typename Type>
  bool operator==(const cloneable_ptr<Type> &A, std::nullptr_t);

  template<typename Type>
  bool operator!=(const cloneable_ptr<Type> &A, std::nullptr_t);


  // --- Immplementation ---

  template < typename T, typename std::enable_if < !has_clone<T>::value, void >::type * >
  std::unique_ptr<T> clone(const T &obj) {
    return std::unique_ptr<T>(new T(obj));
  }

  template<typename T, typename std::enable_if<has_clone<T>::value, void>::type * >
  std::unique_ptr<T> clone(const T &obj) {
    return obj.clone();
  }

  /// \brief Construct std::unique_ptr<T> from rvalue reference
  ///
  /// - If obj.move() is valid, use that;
  /// - else if obj.clone() is valid use that;
  /// - else copy-construct (fails if T is abstract)
  template < typename T, typename std::enable_if < !has_move<T>::value, void >::type * >
  std::unique_ptr<T> clone_move(T &&obj) {
    return clone(std::move(obj));
  }

  /// \brief Construct std::unique_ptr<T> from rvalue reference
  ///
  /// - If obj.move() is valid, use that;
  /// - else if obj.clone() is valid use that;
  /// - else copy-construct (fails if T is abstract)
  template<typename T, typename std::enable_if<has_move<T>::value, void>::type * >
  std::unique_ptr<T> clone_move(T &&obj) {
    return obj.move();
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

  template<typename Type>
  bool operator==(std::nullptr_t, const cloneable_ptr<Type> &B) {
    return B.unique() == nullptr;
  }

  template<typename Type>
  bool operator!=(std::nullptr_t, const cloneable_ptr<Type> &B) {
    return B.unique() != nullptr;
  }

  template<typename Type>
  bool operator==(const cloneable_ptr<Type> &A, std::nullptr_t) {
    return A.unique() == nullptr;
  }

  template<typename Type>
  bool operator!=(const cloneable_ptr<Type> &A, std::nullptr_t) {
    return A.unique() != nullptr;
  }

}

#endif
