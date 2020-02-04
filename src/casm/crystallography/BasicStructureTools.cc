#include "casm/crystallography/BasicStructureTools.hh"
#include "casm/crystallography/BasicStructure.hh"
#include "casm/crystallography/Coordinate.hh"
#include "casm/crystallography/LatticePointWithin.hh"
#include "casm/crystallography/Niggli.hh"
#include "casm/crystallography/UnitCellCoord.hh"
#include "casm/crystallography/Site.hh"
#include "casm/crystallography/SuperlatticeEnumerator.hh"
#include "casm/crystallography/SymTools.hh"
#include "casm/crystallography/SymType.hh"

// TODO: This has to get gobbled into crystallography
#include "casm/basis_set/DoFIsEquivalent.hh"
#include "casm/external/Eigen/Core"
#include "casm/misc/CASM_Eigen_math.hh"
#include <algorithm>
#include <atomic>
#include <unordered_set>
#include <utility>
#include <vector>

namespace {
  using namespace CASM;

  /// Returns false if the given operation does not create an equivalent degree of freedom.
  /// For example, a structure with only zz strain must have uniaxial symmetry, but this won't
  /// be taken into account when generating the point group of the lattice. This method is used
  /// to discard operations from the factor group that aren't compatible with degrees of freedom
  /// that weren't taken into account during its creation.
  bool is_compatible_with_operation(const xtal::SymOp &operation, const std::map<DoFKey, DoFSet> &global_dof_map) {
    for(const auto &dof : global_dof_map) {
      if(!DoFIsEquivalent(dof.second)(operation)) {
        return false;
      }
    }
    return true;
  }

  /// Doubles the number of symmetry operations by multiplying each one with a time reversal flip
  void expand_with_time_reversal(xtal::SymOpVector *growing_point_group) {
    int size = growing_point_group->size();
    for(int ix = 0; ix < size; ++ix) {
      growing_point_group->emplace_back(growing_point_group->at(ix) * xtal::SymOp::time_reversal());
    }
    return;
  }

  /// Adds the translation to every basis site
  std::vector<xtal::Site> apply_translation(std::vector<xtal::Site> &translatable_basis, const xtal::Coordinate &translation) {
    std::vector<xtal::Site> translated_basis;
    for(const xtal::Site &s : translatable_basis) {
      translated_basis.emplace_back(s + translation);
    }

    return translated_basis;
  }

  /// Returns the indexes of basis where you can find sites of the translated basis. If the site can't be found, then the size
  /// of the basis is pushed back
  std::vector<Index> map_translated_basis(const std::vector<xtal::Site> &basis,
                                          const std::vector<xtal::Site> &translatable_basis,
                                          const xtal::Coordinate &translation) {
    assert(basis.size() == translatable_basis.size());
    std::vector<Index> basis_map;

    for(const xtal::Site &s_tb : translatable_basis) {
      Index ix = xtal::find_index_with_translation(basis, s_tb, translation);
      basis_map.push_back(ix);
    }

    return basis_map;
  }

  /// Returns the average distance between the original basis and the transformed basis.
  /// The basis map describes the relationship between basis and transformed basis.
  /// For site at index i in basis, the corresponding index in transformed basis is given
  /// by basis_map[i]. So basis[i] maps to translated_transformed_basis[basis_map[i]].
  /// The lattice is given to ensure that all the sites are consistent with the structure you're
  /// working on, since all sites should have the same lattice.
  /// The mapping translation error is calculated using translated_transformed_basis-basis
  xtal::Coordinate average_mapping_translation_error(const std::vector<xtal::Site> &basis,
                                                     const std::vector<xtal::Site> &translated_transformed_basis,
                                                     const std::vector<Index> &basis_map,
                                                     const xtal::Lattice &structure_lattice) {
    assert(basis.size() == translated_transformed_basis.size());
    assert(basis.size() == basis_map.size());
    xtal::Coordinate average_mapping_translation(structure_lattice);

    for(Index ix = 0; ix < basis_map.size(); ++ix) {
      const xtal::Site &mapped_site = translated_transformed_basis[basis_map[ix]];
      const xtal::Site &mapped_to_site = basis[ix];
      assert(mapped_site.lattice() == mapped_to_site.lattice());
      assert(mapped_site.lattice() == structure_lattice);

      xtal::Coordinate single_mapping_translation = mapped_site.min_translation(mapped_to_site);
      average_mapping_translation += single_mapping_translation;
    }

    average_mapping_translation.cart() *= (1.0 / basis.size());
    return average_mapping_translation;
  }

  /// Generates the factor group by applying every possible translation to the provided
  /// point group, and checking if the structure maps onto itself.
  /// The structure is considered to map onto itself by just looking at the position of the sites,
  /// so it does not take into account orientation of molecules, or global degrees of freedom.
  /// This routine is slow for non primitive structures!
  xtal::SymOpVector make_factor_group_from_point_group_translations(const xtal::BasicStructure &struc, const xtal::SymOpVector &point_group) {
    if(struc.basis().size() == 0) {
      return point_group;
    }

    xtal::SymOpVector factor_group;

    for(const xtal::SymOp &point_group_operation : point_group) {
      // apply the point group operation to tall the sites
      std::vector<xtal::Site> transformed_basis;
      for(const xtal::Site &s : struc.basis()) {
        transformed_basis.emplace_back(point_group_operation * s);
      }

      // Using the symmetrically transformed basis, find all possible
      // translations that MIGHT map the symmetrically transformed basis onto
      // the original basis
      const xtal::Site &reference_site = struc.basis()[0];
      for(const xtal::Site &transformed_site : transformed_basis) {
        // If the types don't match don't even bother with anything else
        if(reference_site.compare_type(transformed_site)) {
          continue;
        }

        xtal::Coordinate translation(struc.lattice());
        translation = reference_site - transformed_site;
        translation.within();

        std::vector<xtal::Site> translated_transformed_basis = apply_translation(transformed_basis, translation);

        // By construction, the current transformed_site matches the first
        // basis site, do the rest of them match too? If not, continue to the
        // next site for a new translation
        std::vector<Index> basis_map = map_translated_basis(struc.basis(), transformed_basis, translation);
        std::unordered_set<Index> unique_mappings(basis_map.begin(), basis_map.end());
        if(unique_mappings.size() != struc.basis().size()) {
          continue;
        }

        // You found a valid translation! BUT, before you construct the
        // symmetry operation, remove some bias from the translation. We want
        // the average mapping error to be zero, but we arbitrarily
        // constructed the translation by taking the first basis site of the
        // structure as a reference. Here we use the average mapping error to
        // correct this.
        xtal::Coordinate mapping_translation_error =
          average_mapping_translation_error(struc.basis(), translated_transformed_basis, basis_map, struc.lattice());
        translation -= mapping_translation_error;

        // Now that the translation has been adjusted, create the symmetry
        // operation and add it if we don't have an equivalent one already
        xtal::SymOp translation_operation = xtal::SymOp::translation_operation(translation.const_cart());
        xtal::SymOp factor_group_operation(translation_operation * point_group_operation);
        xtal::SymOpPeriodicCompare_f fg_op_compare(factor_group_operation, struc.lattice(), CASM::TOL);
        if(std::find_if(factor_group.begin(), factor_group.end(), fg_op_compare) == factor_group.end()) {
          factor_group.push_back(factor_group_operation);
        }
      }
    }

    xtal::close_group<xtal::SymOpPeriodicCompare_f>(&factor_group, struc.lattice(), TOL);
    //!!TODO!! Bring the operations within
    return factor_group;
  }

  xtal::SymOpVector make_translation_group(const xtal::BasicStructure &struc) {
    xtal::SymOpVector identity_group{xtal::SymOp::identity()};
    xtal::SymOpVector translation_group = make_factor_group_from_point_group_translations(struc, identity_group);
    return translation_group;
  }

  /// Given a structure, make it primitive and calculate its factor group. Return the primitive structure and the
  /// factor group of the primitive structure
  std::pair<xtal::BasicStructure, xtal::SymOpVector> make_primitive_factor_group(const xtal::BasicStructure &non_primitive_struc) {
    xtal::BasicStructure primitive_struc = xtal::make_primitive(non_primitive_struc);

    xtal::SymOpVector primitive_point_group = xtal::make_point_group(primitive_struc.lattice());
    if(primitive_struc.is_time_reversal_active()) {
      // Duplicate each symmetry operation so that the second version has time reversal enabled
      ::expand_with_time_reversal(&primitive_point_group);
    }

    xtal::SymOpVector primitive_factor_group = ::make_factor_group_from_point_group_translations(primitive_struc, primitive_point_group);
    return std::make_pair(primitive_struc, primitive_factor_group);
  }

} // namespace

namespace CASM {
  namespace xtal {
    Index find_index(const std::vector<Site> &basis, const Site &test_site) {
      for(Index i = 0; i < basis.size(); i++) {
        if(basis[i].compare(test_site)) {
          return i;
        }
      }
      return basis.size();
    }

    Index find_index_with_translation(const std::vector<Site> &basis, const Site &test_site, const Coordinate &shift) {
      for(Index i = 0; i < basis.size(); i++) {
        if(basis[i].compare(test_site, shift)) {
          return i;
        }
      }
      return basis.size();
    }

    bool is_primitive(const BasicStructure &struc) {
      SymOpVector translation_group = ::make_translation_group(struc);
      // For a primitive structure, the only possible translation is no translation
      return translation_group.size() == 1;
    }

    BasicStructure make_primitive(const BasicStructure &non_primitive_struc) {
      SymOpVector translation_group = ::make_translation_group(non_primitive_struc);
      double minimum_possible_primitive_volume = std::abs(0.5 * non_primitive_struc.lattice().volume() / non_primitive_struc.basis().size());

      // The candidate lattice vectors are the original lattice vectors, plus all the possible translations that map the basis
      // of the non primitive structure onto itself
      std::vector<Eigen::Vector3d> possible_lattice_vectors{non_primitive_struc.lattice()[0], non_primitive_struc.lattice()[1],
                                                            non_primitive_struc.lattice()[2]};
      for(const SymOp &trans_op : translation_group) {
        possible_lattice_vectors.push_back(trans_op.translation);
      }

      // Attempt every combination of vectors, picking one that doesn't have colinearity (minimum volume check), but results in the smallest
      // lattice possible (running lattice volume check)
      double minimum_volume = 2 * non_primitive_struc.lattice().volume();
      Eigen::Vector3d a_vector_primitive, b_vector_primitive, c_vector_primitive;
      for(const Eigen::Vector3d a_vector_candidate : possible_lattice_vectors) {
        for(const Eigen::Vector3d b_vector_candidate : possible_lattice_vectors) {
          for(const Eigen::Vector3d c_vector_candidate : possible_lattice_vectors) {
            double possible_volume = std::abs(triple_product(a_vector_candidate, b_vector_candidate, c_vector_candidate));
            if(possible_volume < minimum_volume && possible_volume > TOL) {
              minimum_volume = possible_volume;
              a_vector_primitive = a_vector_candidate;
              b_vector_primitive = b_vector_candidate;
              c_vector_primitive = c_vector_candidate;
            }
          }
        }
      }
      Lattice non_reduced_form_primitive_lattice(a_vector_primitive, b_vector_primitive, c_vector_primitive);
      Lattice primitive_lattice = niggli(non_reduced_form_primitive_lattice, non_primitive_struc.lattice().tol());

      // The primitive lattice could be noisy, so we smoothen it out to match an integer transformation to the
      // original lattice exactly
      Superlattice prim_to_original = Superlattice::smooth_prim(primitive_lattice, non_primitive_struc.lattice());
      primitive_lattice = prim_to_original.prim_lattice();

      // Fill up the basis
      BasicStructure primitive_struc(primitive_lattice);
      for(Site site_for_prim : non_primitive_struc.basis()) {
        site_for_prim.set_lattice(primitive_struc.lattice(), CART);
        if(find_index(primitive_struc.basis(), site_for_prim) == primitive_struc.basis().size()) {
          site_for_prim.within();
          primitive_struc.set_basis().emplace_back(std::move(site_for_prim));
        }
      }

      //TODO: Do we want this?
      primitive_struc.set_title(non_primitive_struc.title());
      return primitive_struc;
    }

    std::vector<SymOp> make_factor_group(const BasicStructure &struc) {
      auto prim_factor_group_pair =::make_primitive_factor_group(struc);
      const BasicStructure &primitive_struc = prim_factor_group_pair.first;
      const std::vector<SymOp> &primitive_factor_group = prim_factor_group_pair.second;

      auto all_lattice_points = make_lattice_points(primitive_struc.lattice(), struc.lattice(), struc.lattice().tol());
      std::vector<SymOp> point_group = make_point_group(struc.lattice());
      std::vector<SymOp> factor_group;

      for(const SymOp &prim_op : primitive_factor_group) {
        //If the primitive factor group operation with translations removed can't map the original structure's
        //lattice onto itself, then ditch that operation.
        SymOpMatrixCompare_f match_ignoring_translations(prim_op, TOL);
        if(std::find_if(point_group.begin(), point_group.end(), match_ignoring_translations) == point_group.end()) {
          continue;
        }

        //Otherwise take that factor operation, and expand it by adding additional translations within the structure
        for(const UnitCell &lattice_point : all_lattice_points) {
          Coordinate lattice_point_coordinate = make_superlattice_coordinate(lattice_point, primitive_struc.lattice(), struc.lattice());
          factor_group.emplace_back(SymOp::translation_operation(lattice_point_coordinate.cart())*prim_op);
        }
      }

      return factor_group;
    }

  } // namespace xtal
} // namespace CASM
