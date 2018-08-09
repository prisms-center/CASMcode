#include "casm/clex/ScelEnum_impl.hh"

extern "C" {
  CASM::EnumInterfaceBase *make_ScelEnum_interface() {
    return new CASM::EnumInterface<CASM::ScelEnum>();
  }
}

namespace CASM {

  namespace {
    typedef std::vector<std::string>::iterator str_vec_it;
  }

  template class ScelEnumByNameT<true>;
  template class ScelEnumByNameT<false>;

  template ScelEnumByNameT<true>::ScelEnumByNameT(PrimClex &, str_vec_it, str_vec_it);
  template ScelEnumByNameT<false>::ScelEnumByNameT(PrimClex &, str_vec_it, str_vec_it);

  template class ScelEnumByPropsT<true>;
  template class ScelEnumByPropsT<false>;

  template struct CASM_TMP::traits<ScelEnumT<true> >;
  template struct CASM_TMP::traits<ScelEnumT<false> >;
  template class EnumInterface<ScelEnumT<true> >;
  template class EnumInterface<ScelEnumT<false> >;
  template class ScelEnumT<true>;
  template class ScelEnumT<false>;

}
