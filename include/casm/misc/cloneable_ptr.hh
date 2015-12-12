#ifndef cloneable_ptr_HH
#define cloneable_ptr_HH

#include <memory>

namespace notstd {
  
  /// \brief c++11 does not include 'make_unique'
  template<typename T, typename ...Args>
  std::unique_ptr<T> make_unique( Args&& ...args )
  {
    return std::unique_ptr<T>( new T( std::forward<Args>(args)... ) );
  }
  
  
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
    
    typedef Type* pointer;
    typedef Type& reference;
    
    cloneable_ptr(pointer ptr) : 
      m_unique(ptr) {}
    
    cloneable_ptr(const Type& obj) : 
      m_unique(obj.clone()) {}
    
    cloneable_ptr(const cloneable_ptr& other) : 
      m_unique(other->clone()) {}
      
    cloneable_ptr& operator=(cloneable_ptr other) {
      swap(*this, other);
      return *this;
    }
    
    reference operator*() const {
      return *m_unique;
    }
    
    pointer operator->() const {
      return m_unique.operator->();
    }
    
    std::unique_ptr<Type>& unique() {
      return m_unique;
    }
    
    const std::unique_ptr<Type>& unique() const {
      return m_unique;
    }
    
    private:
    
    std::unique_ptr<Type> m_unique;
    
  };
  
  template<typename Type>
  void swap(cloneable_ptr<Type>& A, cloneable_ptr<Type>& B) {
    A.unique().swap(B.unique());
  }
}

#endif
