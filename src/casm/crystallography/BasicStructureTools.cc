#include "casm/crystallography/BasicStructureTools.hh"
#include "casm/crystallography/BasicStructure.hh"
#include "casm/crystallography/Coordinate.hh"
#include "casm/crystallography/Site.hh"
#include "casm/crystallography/SuperlatticeEnumerator.hh"
#include "casm/crystallography/SymType.hh"

// TODO: This has to get gobbled into crystallography
#include "casm/basis_set/DoFIsEquivalent.hh"
#include <algorithm>
#include <unordered_set>
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
      Index ix = xtal::find_index(basis, s_tb, translation);
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

    xtal::close_group(&factor_group, struc.lattice());
    return factor_group;
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

    Index find_index(const std::vector<Site> &basis, const Site &test_site, const Coordinate &shift) {
      for(Index i = 0; i < basis.size(); i++) {
        if(basis[i].compare(test_site, shift)) {
          return i;
        }
      }
      return basis.size();
    }

    /* bool is_primitive(const BasicStructure& struc){ */
    /*   SymGroup valid_translations, identity_group; */
    /*   identity_group.push_back(CASM::SymOp()); */
    /*   _generate_factor_group_slow(valid_translations, identity_group, false); */
    /*   if(valid_translations.size() == 1) */
    /*     return true; */

    /*   return false; */
    /* } */

  } // namespace xtal
} // namespace CASM
