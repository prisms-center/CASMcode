#ifndef CASM_jsonIO_clex
#define CASM_jsonIO_clex

#include "casm/casm_io/jsonParser.hh"
#include "casm/CASM_global_definitions.hh"

#include "casm/clex/ChemicalReference.hh"

namespace CASM {

  class Structure;

  /**
   * \ingroup ProjectIO
   *
   * @{
   */

  // --- ChemicalReferenceState -------------------

  /// \brief Write ChemicalReferenceState to: '{"A" : X, "B" : X, ..., "energy_per_species" : X }'
  jsonParser &to_json(const ChemicalReferenceState &ref_state, jsonParser &json);

  /// \brief Read ChemicalReferenceState from: '{"A" : X, "B" : X, ..., "energy_per_species" : X }'
  template<>
  ChemicalReferenceState from_json<ChemicalReferenceState>(const jsonParser &json);

  /// \brief Read ChemicalReferenceState from: '{"A" : X, "B" : X, ..., "energy_per_species" : X }'
  void from_json(ChemicalReferenceState &ref_state, const jsonParser &json);


  // --- HyperPlaneReference -------------------

  jsonParser &to_json(const HyperPlaneReference &ref, jsonParser &json);

  template<>
  struct jsonConstructor<HyperPlaneReference> {

    static HyperPlaneReference from_json(const jsonParser &json,
                                         HyperPlaneReference::InputFunction f);
  };

  void from_json(HyperPlaneReference &ref,
                 const jsonParser &json,
                 HyperPlaneReference::InputFunction f);


  // --- ChemicalReference -------------------

  /// \brief Write chemical reference
  jsonParser &to_json(const ChemicalReference &ref, jsonParser &json);

  /// \brief Read chemical reference from one of 3 alternative forms
  std::pair<Eigen::VectorXd, std::vector<ChemicalReferenceState> >
  one_chemical_reference_from_json(const Structure &prim,
                                   const jsonParser &json);

  /// \brief Read chemical reference from JSON
  template<>
  struct jsonConstructor<ChemicalReference> {

    static ChemicalReference from_json(const jsonParser &json,
                                       const Structure &prim,
                                       double tol = 1e-14);
  };

  /// \brief Read chemical reference from JSON
  void from_json(ChemicalReference &ref,
                 const jsonParser &json,
                 const Structure &prim,
                 double tol = 1e-14);

  /** @} */

}

#endif
