#ifndef CASM_MappedPropertiesTools
#define CASM_MappedPropertiesTools
namespace CASM {
struct MappedProperties;
class PermuteIterator;
///
MappedProperties copy_apply(PermuteIterator const &op,
                            MappedProperties const &props);

}  // namespace CASM

#endif
