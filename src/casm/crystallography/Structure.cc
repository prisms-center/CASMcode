#include "casm/crystallography/Structure.hh"

#include <sys/stat.h>
#include <sys/types.h>

#include <boost/filesystem.hpp>
#include <boost/filesystem/fstream.hpp>
#include <exception>
#include <memory>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

#include "casm/basis_set/Adapter.hh"
#include "casm/basis_set/DoF.hh"
#include "casm/basis_set/DoFIsEquivalent_impl.hh"
#include "casm/basis_set/DoFTraits.hh"
#include "casm/casm_io/Log.hh"
#include "casm/container/algorithm.hh"
#include "casm/crystallography/Adapter.hh"
#include "casm/crystallography/BasicStructure.hh"
#include "casm/crystallography/BasicStructureTools.hh"
#include "casm/crystallography/DoFSet.hh"
#include "casm/crystallography/IntegralCoordinateWithin.hh"
#include "casm/crystallography/SymTools.hh"
#include "casm/crystallography/SymType.hh"
#include "casm/external/Eigen/src/Core/Matrix.h"
#include "casm/misc/algorithm.hh"
#include "casm/symmetry/SupercellSymInfo.hh"
#include "casm/symmetry/SymBasisPermute.hh"
#include "casm/symmetry/SymGroupRep.hh"
#include "casm/symmetry/SymMatrixXd.hh"
#include "casm/symmetry/SymOp.hh"
#include "casm/symmetry/SymPermutation.hh"

namespace CASM {
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Structure::Structure(const fs::path &filepath) {
  if (!fs::exists(filepath)) {
    err_log() << "Error in Structure::Structure(const fs::path &filepath)."
              << std::endl;
    err_log() << "  File does not exist at: " << filepath << std::endl;
    exit(1);
  }
  fs::ifstream infile(filepath);

  BasicStructure read_struc = BasicStructure::from_poscar_stream(infile);
  m_structure_ptr = std::make_shared<BasicStructure>(read_struc);
  this->generate_factor_group();
}

// TODO: This is awful
Structure::Structure()
    : m_structure_ptr(std::make_shared<BasicStructure>(BasicStructure())) {}

Structure::Structure(const Structure &RHS)
    : m_structure_ptr(RHS.m_structure_ptr) {
  copy_attributes_from(RHS);
};

Structure::Structure(const BasicStructure &base)
    : m_structure_ptr(std::make_shared<BasicStructure>(base)) {
  this->generate_factor_group();
}

Structure &Structure::operator=(const Structure &RHS) {
  m_structure_ptr = RHS.m_structure_ptr;
  copy_attributes_from(RHS);
  return *this;
}

Structure::operator const BasicStructure &() const { return this->structure(); }

//***********************************************************

void Structure::copy_attributes_from(const Structure &RHS) {
  m_basis_perm_rep_ID = RHS.m_basis_perm_rep_ID;      // this *should* work
  m_site_dof_symrep_IDs = RHS.m_site_dof_symrep_IDs;  // this *should* work
  m_factor_group = RHS.m_factor_group;
  m_factor_group.set_lattice(this->lattice());
}

//***********************************************************

void Structure::generate_factor_group() {
  m_factor_group.clear();
  m_factor_group.set_lattice(this->lattice());

  // Don't copy a MasterSymGroup or you'll have bad luck
  xtal::SymOpVector factor_group_operations =
      xtal::make_factor_group(*this, this->lattice().tol());
  for (const xtal::SymOp &op : factor_group_operations) {
    m_factor_group.push_back(adapter::Adapter<CASM::SymOp, xtal::SymOp>()(op));
  }

  m_factor_group.sort();

  _generate_basis_symreps();
  _generate_global_symreps();

  return;
}

//************************************************************
const MasterSymGroup &Structure::factor_group() const { return m_factor_group; }

//************************************************************
const SymGroup &Structure::point_group() const {
  return factor_group().point_group();
}

//***********************************************************

SymGroupRep const *Structure::basis_permutation_symrep() const {
  return &(factor_group().representation(basis_permutation_symrep_ID()));
}

//***********************************************************

SymGroupRepID Structure::basis_permutation_symrep_ID() const {
  return m_basis_perm_rep_ID;
}

//************************************************************

std::vector<SymGroupRepID> Structure::occupant_symrep_IDs() const {
  return this->m_occupant_symrep_IDs;
}

std::vector<std::map<DoFKey, SymGroupRepID>> Structure::site_dof_symrep_IDs()
    const {
  return this->m_site_dof_symrep_IDs;
}

SymGroupRepID Structure::global_dof_symrep_ID(
    const std::string dof_name) const {
  return this->m_global_dof_symrep_IDs.at(dof_name);
}

void Structure::_reset_occupant_symrep_IDs() {
  this->m_occupant_symrep_IDs.clear();
  for (const Site &s : this->basis()) {
    this->m_occupant_symrep_IDs.emplace_back(
        SymGroupRepID::identity(s.allowed_occupants().size()));
  }
  return;
}

void Structure::_reset_site_dof_symrep_IDs() {
  this->m_site_dof_symrep_IDs =
      std::vector<std::map<DoFKey, SymGroupRepID>>(this->basis().size());
}

// This function gets the permutation representation of the
// factor group operations of the structure. It first applies
// the factor group operation to the structure, and then tries
// to map the new position of the basis atom to the various positions
// before symmetry was applied. It only checks the positions after
// it brings the basis within the crystal.
void Structure::_generate_basis_symreps() {
  std::string clr(100, ' ');
  if (factor_group().size() <= 0) {
    err_log() << "ERROR in generate_basis_permutation_representation"
              << std::endl;
    err_log() << "Factor group is empty." << std::endl;
    ;
    exit(1);
  }

  std::vector<UnitCellCoord> sitemap;

  m_basis_perm_rep_ID = m_factor_group.allocate_representation();

  this->_reset_site_dof_symrep_IDs();
  for (Index b = 0; b < basis().size(); ++b) {
    for (auto const &dof : basis()[b].dofs()) {
      this->m_site_dof_symrep_IDs[b][dof.first] =
          this->m_factor_group.allocate_representation();
    }
  }

  // The sitemap specifies that op*basis(b) -> sitemap[b] (which is a
  // UnitCellCoord) The new dofs of site specified by UCC sitemap[b] will be
  // transformations of the dofs that previously resided at basis(b). As such,
  // for dofs, we use the inverse permutation and
  //   basis()[b].symrep(doftype.name()) =
  //   basis()[b].dof(doftype.name()).basis().transpose()
  //                                       * doftype.symop_to_matrix(op)
  //                                       * basis()[sitemap[b].sublattice()].dof(doftype.name().basis())
  this->_reset_occupant_symrep_IDs();
  for (Index s = 0; s < m_factor_group.size(); ++s) {
    auto const &op = m_factor_group[s];

    // It seems that ideally xtal::make_factor_group would provide the mappings
    // that are being re-determined here. Since the factor group has already
    // been determined at this point using xtal::make_factor_group, no mapping
    // failures should be tolerated here.
    // Fundamentally, both xtal::make_factor_group and xtal::symop_site_map use
    // something equivalent to `basis[i].min_dist(test_site) < lattice().tol()`
    // to find the basis site mapping and therefore should agree. This should
    // not re-check that `cart2frac` is integer.
    sitemap = xtal::symop_site_map(op, *this);
    Eigen::Matrix3l point_mat = lround(cart2frac(op.matrix(), lattice()));
    op.set_rep(m_basis_perm_rep_ID, SymBasisPermute(point_mat, sitemap));

    for (Index b = 0; b < basis().size(); ++b) {
      // copy_aply(symop,dofref_from) = P.permute(dofref_to);
      auto const &dofref_to = basis()[sitemap[b].sublattice()].occupant_dof();
      auto const &dofref_from = basis()[b].occupant_dof();

      auto &symrep_from = this->m_occupant_symrep_IDs[b];
      OccupantDoFIsEquivalent<xtal::Molecule> eq(dofref_from);

      if (eq(adapter::Adapter<xtal::SymOp, CASM::SymOp>()(op), dofref_to)) {
        if (symrep_from.is_identity()) {
          if (!eq.perm().is_identity()) {
            symrep_from = m_factor_group.allocate_representation();
            Index s2;
            for (s2 = 0; s2 < s; ++s2) {
              m_factor_group[s2].set_rep(
                  symrep_from,
                  SymPermutation(sequence<Index>(0, dofref_from.size() - 1)));
            }
            m_factor_group[s2].set_rep(symrep_from,
                                       SymPermutation(eq.perm().inverse()));
          }
        } else {
          op.set_rep(symrep_from, SymPermutation(eq.perm().inverse()));
        }
      } else
        throw std::runtime_error(
            "In Structure::_generate_basis_symreps(), Sites originally "
            "identified as equivalent cannot be mapped by symmetry.");
    }

    for (auto const &dof_key : xtal::continuous_local_dof_types(*this)) {
      for (Index from_b = 0; from_b < basis().size(); ++from_b) {
        if (!basis()[from_b].has_dof(dof_key)) continue;

        xtal::DoFSet const &_dofref_from = basis()[from_b].dof(dof_key);

        Index to_b = sitemap[from_b].sublattice();
        xtal::DoFSet const &_dofref_to = basis()[to_b].dof(dof_key);

        // Transform the xtal::SiteDoFSet to the CASM::DoFSet version
        CASM::DoFSet dofref_from =
            adapter::Adapter<CASM::DoFSet, xtal::SiteDoFSet>()(
                _dofref_from, SymGroupRepID(), from_b);

        CASM::DoFSet dofref_to =
            adapter::Adapter<CASM::DoFSet, xtal::SiteDoFSet>()(
                _dofref_to, SymGroupRepID(), to_b);

        DoFIsEquivalent eq(dofref_from);
        // TODO
        // Calling the adapter here, because we said we don't want anything
        // outside of crystallography to invoke crystallography/Adapter.hh
        if (!eq(adapter::Adapter<xtal::SymOp, CASM::SymOp>()(op), dofref_to)) {
          throw std::runtime_error(
              "While generating symmetry representation for local DoF \"" +
              dof_key +
              "\", a symmetry operation was identified that invalidates the "
              "degree of freedom. " +
              "Degrees of freedom must be fully specified before performing "
              "symmetry analyses.");
        }

        SymGroupRepID from_symrep_ID =
            this->site_dof_symrep_IDs()[from_b][dof_key];
        op.set_rep(from_symrep_ID, SymMatrixXd(eq.U()));
      }
    }
  }

  return;
}

//***********************************************************

// TODO: Simplify the DoF equivalence checks. You can probably completely erase
// the comparators in basis_set in favor of the ones in the xtal namespace
void Structure::_generate_global_symreps() {
  if (factor_group().size() <= 0) {
    err_log() << "ERROR in generate_global_dof_representations" << std::endl;
    err_log() << "Factor group is empty." << std::endl;
    exit(1);
  }
  for (auto const &name_dof_pr : this->structure().global_dofs()) {
    std::string dof_name = name_dof_pr.first;
    const xtal::DoFSet &dof = name_dof_pr.second;

    this->m_global_dof_symrep_IDs[dof_name] =
        this->factor_group().allocate_representation();

    for (auto const &op : m_factor_group) {
      /* DoFIsEquivalent eq(dof.second); */
      xtal::DoFSetIsEquivalent_f dof_equals(dof, TOL);

      xtal::DoFSet transformed_dof = sym::copy_apply(
          adapter::Adapter<xtal::SymOp, CASM::SymOp>()(op), dof);
      if (!dof_equals(transformed_dof)) {
        throw std::runtime_error(
            "While generating symmetry representation for global DoF \"" +
            dof_name +
            "\", a symmetry operation was identified that invalidates the "
            "degree of freedom. " +
            "Degrees of freedom must be fully specified before performing "
            "symmetry analyses.");
      }
      Eigen::MatrixXd basis_change_representation;
      try {
        basis_change_representation = xtal::dofset_transformation_matrix(
            dof.basis(), transformed_dof.basis(), TOL);
      } catch (std::runtime_error &e) {
        throw std::runtime_error(std::string(e.what()) +
                                 " Attempted to make representation for " +
                                 dof_name + ".");
      }
      op.set_rep(this->m_global_dof_symrep_IDs.at(dof_name),
                 SymMatrixXd(basis_change_representation));
    }
  }
}
//***********************************************************

/// Returns 'converter' which converts Site::site_occupant indices to 'mol_list'
/// indices:
///   mol_list_index = converter[basis_site][site_occupant_index]
std::vector<std::vector<Index>> make_index_converter(
    const Structure &struc, std::vector<xtal::Molecule> mol_list) {
  std::vector<std::vector<Index>> converter(struc.basis().size());

  for (Index i = 0; i < struc.basis().size(); i++) {
    converter[i].resize(struc.basis()[i].occupant_dof().size());

    for (Index j = 0; j < struc.basis()[i].occupant_dof().size(); j++) {
      converter[i][j] =
          CASM::find_index(mol_list, struc.basis()[i].occupant_dof()[j]);
    }
  }

  return converter;
}

/// Returns 'converter' which converts Site::site_occupant indices to
/// 'mol_name_list' indices:
///   mol_name_list_index = converter[basis_site][site_occupant_index]
std::vector<std::vector<Index>> make_index_converter(
    const Structure &struc, std::vector<std::string> mol_name_list) {
  std::vector<std::vector<Index>> converter(struc.basis().size());

  for (Index i = 0; i < struc.basis().size(); i++) {
    converter[i].resize(struc.basis()[i].occupant_dof().size());

    for (Index j = 0; j < struc.basis()[i].occupant_dof().size(); j++) {
      converter[i][j] = CASM::find_index(
          mol_name_list, struc.basis()[i].occupant_dof()[j].name());
    }
  }

  return converter;
}

/// Returns 'converter_inverse' which converts 'mol_name_list' indices to
/// Site::site_occupant indices:
///  site_occupant_index = converter_inverse[basis_site][mol_name_list_index]
///
/// If mol is not allowed on basis_site, return
/// struc.basis()[basis_site].occupant_dof().size()
std::vector<std::vector<Index>> make_index_converter_inverse(
    const Structure &struc, std::vector<std::string> mol_name_list) {
  std::vector<std::vector<Index>> converter_inv(struc.basis().size());

  for (Index i = 0; i < struc.basis().size(); i++) {
    converter_inv[i].resize(mol_name_list.size());

    std::vector<std::string> site_occ_name_list;
    for (Index j = 0; j < struc.basis()[i].occupant_dof().size(); j++) {
      site_occ_name_list.push_back(struc.basis()[i].occupant_dof()[j].name());
    }

    for (Index j = 0; j < mol_name_list.size(); j++) {
      converter_inv[i][j] =
          CASM::find_index(site_occ_name_list, mol_name_list[j]);
    }
  }

  return converter_inv;
}

std::map<DoFKey, DoFSetInfo> global_dof_info(Structure const &_struc) {
  std::map<DoFKey, DoFSetInfo> result;
  for (auto const &dof : _struc.structure().global_dofs()) {
    result.emplace(dof.first,
                   adapter::Adapter<CASM::DoFSet, xtal::DoFSet>()(
                       dof.second, _struc.global_dof_symrep_ID(dof.first))
                       .info());
  }
  return result;
}

std::map<DoFKey, std::vector<DoFSetInfo>> local_dof_info(
    Structure const &_struc) {
  std::map<DoFKey, std::vector<DoFSetInfo>> result;

  for (DoFKey const &type : xtal::continuous_local_dof_types(_struc)) {
    std::vector<CASM::DoFSetInfo> tresult(
        _struc.basis().size(),
        CASM::DoFSetInfo(
            SymGroupRepID(),
            Eigen::MatrixXd::Zero(DoF::BasicTraits(type).dim(), 0)));

    for (Index b = 0; b < _struc.basis().size(); ++b) {
      if (_struc.basis()[b].has_dof(type)) {
        const auto &dofset = _struc.basis()[b].dof(type);
        tresult[b] = adapter::Adapter<CASM::DoFSet, xtal::SiteDoFSet>()(
                         dofset, _struc.site_dof_symrep_IDs()[b].at(type), b)
                         .info();
      }
    }
    result.emplace(type, std::move(tresult));
  }
  return result;
}

SupercellSymInfo make_supercell_sym_info(Structure const &prim,
                                         Lattice const &super_lattice) {
  // Structure data needs to be reorganized for SupercellSymInfo construction

  // map of global DoFKey -> SymGroupRepID
  std::map<DoFKey, SymGroupRepID> global_dof_symrep_IDs;
  for (auto const &key : xtal::global_dof_types(prim)) {
    global_dof_symrep_IDs.emplace(
        std::make_pair(key, prim.global_dof_symrep_ID(key)));
  }

  // map of site DoFKey -> std::vector<SymGroupRepID>
  std::map<DoFKey, std::vector<SymGroupRepID>> local_dof_symrep_IDs;
  for (auto const &key : xtal::continuous_local_dof_types(prim)) {
    std::vector<SymGroupRepID> treps(prim.basis().size());
    for (Index b = 0; b < prim.basis().size(); ++b) {
      if (prim.basis()[b].has_dof(key))
        treps[b] = prim.site_dof_symrep_IDs()[b][key];
    }
    local_dof_symrep_IDs.emplace(std::make_pair(key, std::move(treps)));
  }

  return SupercellSymInfo(
      prim.lattice(), super_lattice, prim.basis().size(), prim.factor_group(),
      prim.basis_permutation_symrep_ID(), global_dof_symrep_IDs,
      prim.occupant_symrep_IDs(), local_dof_symrep_IDs);
}
}  // namespace CASM
