#ifndef PERMUTATIONIO_HH
#define PERMUTATIONIO_HH

namespace CASM {
class jsonParser;
class Permutation;

jsonParser &to_json(const Permutation &value, jsonParser &json);
void from_json(Permutation &value, const jsonParser &json);

}  // namespace CASM

#endif
