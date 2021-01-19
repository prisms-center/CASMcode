#ifndef IDHIERARCHY_HH
#define IDHIERARCHY_HH

#include <iostream>

namespace CASM {

template <typename Base>
class HierarchyID;
template <typename Derived, typename Base>
class DerivedID;

template <typename Base>
class HierarchyID {
 protected:
  static int new_class_ID() {
    static int Nclass = 0;
    Base::extend_hierarchy();
    // std::cout << "Nclass is " << Nclass << '\n';
    return Nclass++;
  };
};

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

template <typename Derived, typename Base>
class DerivedID : protected HierarchyID<Base> {
 protected:
  // static bool table_empty;
  static int get_class_ID() {
    static int ID_val = HierarchyID<Base>::new_class_ID();
    static bool table_empty = true;
    if (table_empty) {
      table_empty = false;
      Derived::fill_dispatch_table();
    }
    return ID_val;
  };
};

/*  template< typename Derived, typename Base >
bool DerivedID<Derived, Base> :: table_empty(true);
*/
}  // namespace CASM
#endif
