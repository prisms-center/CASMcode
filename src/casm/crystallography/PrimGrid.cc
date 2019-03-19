#include "casm/crystallography/PrimGrid.hh"

#include <iostream>
#include <cmath>
#include <cassert>

#include "casm/misc/CASM_Eigen_math.hh"
#include "casm/container/Counter.hh"
#include "casm/casm_io/Log.hh"

#include "casm/crystallography/Lattice.hh"
#include "casm/crystallography/Coordinate.hh"
#include "casm/crystallography/UnitCellCoord.hh"

#include "casm/symmetry/SymGroupRep.hh"
#include "casm/symmetry/SymPermutation.hh"
#include "casm/symmetry/SymBasisPermute.hh"

namespace CASM {
  PrimGrid::PrimGrid(const Lattice &p_lat, const Lattice &s_lat, Index NB) {
    m_lat[static_cast<int>(PRIM)] = &p_lat;
    m_lat[static_cast<int>(SCEL)] = &s_lat;

    m_NB = NB;

    Eigen::Matrix3d dtrans_mat(p_lat.lat_column_mat().inverse()*s_lat.lat_column_mat());

    for(int i = 0; i < 3; i++) {
      for(int j = 0; j < 3; j++) {
        m_trans_mat(i, j) = round(dtrans_mat(i, j));
      }
    }
    if(m_trans_mat.determinant() == 0) {
      default_err_log() << "CRITICAL ERROR:  Attempting to construct a PrimGrid for a superlattice that is smaller than lattice passed as its prim.\n"
                        << "floating-point transformation matrix:\n" << dtrans_mat << "\n was rounded to integer matrix:\n" << m_trans_mat << "\n\n"
                        << "This usage of PrimGrid is not supported. Exiting...\n" << "Prim Lat\n" << p_lat.lat_column_mat() << "Scel Lat\n" << s_lat.lat_column_mat();
      assert(0);
      exit(1);
    }

    matrix_type Smat, V;
    smith_normal_form(m_trans_mat, m_U, Smat, V);

    m_invU = inverse(m_U);

    m_S = Smat.diagonal();

    m_stride[0] = m_S[0];
    m_stride[1] = m_S[0] * m_S[1];
    m_N_vol = m_S[0] * m_S[1] * m_S[2];

    m_trans_mat = m_U * Smat * V;

    Smat(0, 0) = m_N_vol / Smat(0, 0);
    Smat(1, 1) = m_N_vol / Smat(1, 1);
    Smat(2, 2) = m_N_vol / Smat(2, 2);
    m_plane_mat = inverse(V) * Smat * inverse(m_U);


  }

  //**********************************************************************************************
  // Constructor for forcing specific choice of 'U' matrix.  Use only in very specific cases (such as applying symmetry to a PrimGrid
  PrimGrid::PrimGrid(const Lattice &p_lat,
                     const Lattice &s_lat,
                     const Eigen::Ref<const PrimGrid::matrix_type> &U,
                     const Eigen::Ref<const PrimGrid::matrix_type> &Smat,
                     Index NB) : m_U(U) {
    m_lat[static_cast<int>(PRIM)] = &p_lat;
    m_lat[static_cast<int>(SCEL)] = &s_lat;

    m_NB = NB;

    Eigen::Matrix3d dtrans_mat(p_lat.lat_column_mat().inverse()*s_lat.lat_column_mat());

    for(int i = 0; i < 3; i++) {
      for(int j = 0; j < 3; j++) {
        m_trans_mat(i, j) = round(dtrans_mat(i, j));
      }
    }

    if(m_trans_mat.determinant() == 0) {
      default_err_log() << "CRITICAL ERROR:  Attempting to construct a PrimGrid for a superlattice that is smaller than lattice passed as its prim.\n"
                        << "floating-point transformation matrix:\n" << dtrans_mat << "\n was rounded to integer matrix:\n" << m_trans_mat << "\n\n"
                        << "This usage of PrimGrid is not supported. Exiting...\n";
      assert(0);
      exit(1);
    }

    m_invU = inverse(m_U);

    for(int i = 0; i < 3; i++) {
      m_S[i] = Smat(i, i);
    }

    //default_err_log() << "Smith decomposition is:\n";
    m_stride[0] = m_S[0];
    m_stride[1] = m_S[0] * m_S[1];
    m_N_vol = m_S[0] * m_S[1] * m_S[2];

    m_plane_mat = adjugate(m_trans_mat);


  }

  //**********************************************************************************************

  const Lattice &PrimGrid::prim_lattice() const {
    return *m_lat[static_cast<int>(PRIM)];
  }

  //**********************************************************************************************

  const Lattice &PrimGrid::scel_lattice() const {
    return *m_lat[static_cast<int>(SCEL)];
  }

  //**********************************************************************************************

  const Lattice &PrimGrid::lattice(CELL_TYPE lat_mode) const {
    return *m_lat[static_cast<int>(lat_mode)];
  }


  //**********************************************************************************************

  const Eigen::DiagonalWrapper<const PrimGrid::vector_type> PrimGrid::matrixS()const {
    return m_S.asDiagonal();
  }
  //**********************************************************************************************
  Index PrimGrid::find(const Coordinate &_coord) const {

    auto lambda = [](double val) {
      return floor(val);
    };
    auto frac((prim_lattice().inv_lat_column_mat()*_coord.cart()).array() + _coord.lattice().tol());
    UnitCell ijk(frac.unaryExpr(lambda).matrix().cast<long>());

    return find(ijk);
  }

  //**********************************************************************************************

  Index PrimGrid::find(const UnitCell &_unitcell) const {
    UnitCell mnp = to_canonical(within(_unitcell));
    return mnp[0] + mnp[1] * m_stride[0] + mnp[2] * m_stride[1];
  }

  //**********************************************************************************************

  Index PrimGrid::find_cart(const Eigen::Ref<const Eigen::Vector3d> &_cart_coord) const {
    return find(Coordinate(_cart_coord, prim_lattice(), CART));
  }

  //**********************************************************************************************

  /*
    /// The m_plane_mat works because of the following:
    ///
    ///     prim_frac_coord = m_trans_mat * scel_frac_coord
    ///
    /// m_trans_mat is integer, and (if scel_frac_coord is a lattice translation) then
    /// prim_frac_coord is also integer.
    ///
    /// transf_mat.inverse() does the inverse mapping, but it is NOT integer.  However,
    ///
    ///      m_plane_mat = transf_mat.determinant()*transf_mat.inverse()
    ///
    /// is integer. This is because it is simply the adjugate matrix
    /// (i.e., the transpose of the cofactor matrix, which can be obtained
    ///  without using division).
    ///
    /// This means that
    ///
    ///       m_plane_mat*prim_frac_coord = m_N_vol*scel_frac_coord
    ///
    /// again, for a lattice translation, the left hand side is integer, so the
    /// right hand side must be also.  Moreover, the elements of the RHS should
    /// be on the interval [0,m_N_vol-1] if it is within the supercell.
  */

  /// map a UnitCell inside the supercell
  UnitCell PrimGrid::within(const UnitCell &ijk)const {

    vector_type vec2 = m_plane_mat * ijk;

    vec2[0] = ((vec2[0] % m_N_vol) + m_N_vol) % m_N_vol;
    vec2[1] = ((vec2[1] % m_N_vol) + m_N_vol) % m_N_vol;
    vec2[2] = ((vec2[2] % m_N_vol) + m_N_vol) % m_N_vol;

    return (m_trans_mat * vec2) / m_N_vol;
  }

  /// map a UnitCellCoord inside the supercell
  UnitCellCoord PrimGrid::within(const UnitCellCoord &bijk)const {

    return UnitCellCoord(bijk.unit(), bijk.sublat(), within(bijk.unitcell()));
  }

  //**********************************************************************************************

  Coordinate PrimGrid::coord(const UnitCell &ijk, CELL_TYPE lat_mode)const {

    Coordinate tcoord(ijk.cast<double>(), prim_lattice(), FRAC);

    tcoord.set_lattice(lattice(lat_mode), CART);
    return tcoord;
  }

  //**********************************************************************************************

  Coordinate PrimGrid::coord(Index l, CELL_TYPE lat_mode)const {
    return coord(unitcell(l), lat_mode);
  }

  //**********************************************************************************************

  UnitCell PrimGrid::unitcell(Index i) const {
    assert(i >= 0 && i < m_N_vol && "PrimGrid::uccoord(Index i) -> index 'i' out of range!");

    UnitCell mnp((i % m_stride[1]) % m_stride[0],
                 (i % m_stride[1]) / m_stride[0],
                 i / m_stride[1]);
    return from_canonical(mnp);
  }

  //**********************************************************************************************

  SymGroupRepID PrimGrid::make_permutation_representation(const SymGroup &group, SymGroupRepID basis_permute_ID)const {

    SymGroupRepID perm_rep_ID = group.allocate_representation();
    PrimGrid::matrix_type mat_mnp;
    UnitCell mnp_shift;
    Index old_l, new_l;
    for(Index ng = 0; ng < group.size(); ng++) {
      SymBasisPermute const *rep = group[ng].get_basis_permute_rep(basis_permute_ID);
      if(!rep) {
        default_err_log() << "CRITICAL ERROR: In PrimGrid::make_permutation_representation, BasisPermute representation is incorrectly initialized!\n"
                          << "                basis_permute_ID is " << basis_permute_ID << " and op index is " << group[ng].index() << "\n"
                          << "                Exiting...\n";

        exit(1);
      }

      mat_mnp = m_invU * rep->matrix() * m_U;

      std::vector<UnitCellCoord> const &b_permute = rep->data();

      std::vector<Index> ipermute(b_permute.size()*size());
      //default_err_log() << "PRINTING b_permute array for op " << ng << ":\n";
      //begin loop over sites
      for(Index nb = 0; nb < b_permute.size(); nb++) {
        //default_err_log() << b_permute.at(nb) << '\n';
        mnp_shift = to_canonical(b_permute.at(nb).unitcell());

        vector_type old_mnp, new_mnp;
        EigenCounter<vector_type> mnp_count(vector_type::Zero(), m_S - vector_type::Ones(), vector_type::Ones());
        for(; mnp_count.valid(); ++mnp_count) {
          new_mnp = mat_mnp * mnp_count() + mnp_shift;

          //map within bounds
          for(int i = 0; i < 3; i++) {
            new_mnp[i] = ((new_mnp[i] % m_S[i]) + m_S[i]) % m_S[i];
          }

          old_l = mnp_count[0] + mnp_count[1] * m_stride[0] + mnp_count[2] * m_stride[1] + nb * size();
          new_l = new_mnp[0] + new_mnp[1] * m_stride[0] + new_mnp[2] * m_stride[1] + b_permute.at(nb).sublat() * size();
          assert(old_l < b_permute.size()*size() && new_l < b_permute.size()*size());
          // We have found uccoord(new_l) = symop*uccoord(old_l) -- this describes how indexing of the uccoordinates change
          // However, the indexing of the uccoords remains fixed, and we want to describe the permutation of something *at* the sites,
          // like an occupation bit. So if uccord(old_l) transforms into uccoord(new_l), we know that the object originally at 'new_l'
          // was effectively transformed into the object that was originally at 'old_l'. (in other words, the objects permute inversely to the labels)
          //             i.e., new_occ(new_l) = old_occ(old_l)  --> this matches our permutation convention, which is specified as
          //                   new_occ(new_l) = old_occ(ipermute[new_l])
          // and thus:
          ipermute[new_l] = old_l;
        }
      }//end loop over sites
      //default_err_log() << "\\end " << ng << '\n';
      group[ng].set_rep(perm_rep_ID, SymPermutation(ipermute));
    }
    return perm_rep_ID;
  }

  //**********************************************************************************************
  // for std::vector<std::vector<int> > perms=prim_grid.make_translation_permutations(NB);
  // perms[l] is the Supercell permutation that results when sites are translated by
  // prim_grid.uccoord(l).  In other words, site 'l' is translated so that (i,j,k)=(0,0,0)

  std::vector<Permutation> PrimGrid::make_translation_permutations(Index NB)const {
    std::vector<Permutation > perms;
    perms.reserve(size());
    std::vector<Index> ipermute(NB * size(), 0);
    UnitCell shift_mnp, mnp;
    Index new_l;
    for(Index shift_l = 0; shift_l < size(); shift_l++) {
      //shift_mnp describes translation by uccoord(shift_l)
      shift_mnp[0] = (shift_l % m_stride[1]) % m_stride[0];
      shift_mnp[1] = (shift_l % m_stride[1]) / m_stride[0];
      shift_mnp[2] =   shift_l / m_stride[1];
      for(Index old_l = 0; old_l < size(); old_l++) {
        mnp[0] = (old_l % m_stride[1]) % m_stride[0];
        mnp[1] = (old_l % m_stride[1]) / m_stride[0];
        mnp[2] = old_l / m_stride[1];

        mnp[0] = (mnp[0] + shift_mnp[0]) % m_S[0];
        mnp[1] = (mnp[1] + shift_mnp[1]) % m_S[1];
        mnp[2] = (mnp[2] + shift_mnp[2]) % m_S[2];

        for(Index nb = 0; nb < NB; nb++) {
          new_l = mnp[0] + mnp[1] * m_stride[0] + mnp[2] * m_stride[1];

          // See comments in PrimGrid::make_permutation_representation() above
          ipermute[new_l + nb * size()] = old_l + nb * size();
        }
      }
      perms.push_back(Permutation(ipermute));
    }
    return perms;
  }
  // private functions:

  //**********************************************************************************************
  /// Convert UnitCell (ijk) to canonical UnitCell (mnp)
  /// mnp = m_invU * ijk
  UnitCell PrimGrid::to_canonical(const UnitCell &ijk) const {
    UnitCell mnp(0, 0, 0);
    for(int i = 0; i < 3; i++) {
      for(int j = 0; j < 3; j++) {
        mnp[i] += m_invU(i, j) * ijk[j];
      }
      //Map within bounds
      mnp[i] = ((mnp[i] % m_S[i]) + m_S[i]) % m_S[i];
    }

    return mnp;
  }

  //**********************************************************************************************

  /// Convert canonical UnitCell (mnp) to UnitCell (ijk)
  /// U*mnp = ijk
  UnitCell PrimGrid::from_canonical(const UnitCell &mnp) const {
    UnitCell ijk(0, 0, 0);
    for(int i = 0; i < 3; i++) {
      for(int j = 0; j < 3; j++) {
        ijk[i] += m_U(i, j) * mnp[j];
      }
    }
    return within(ijk);
  };

  //==============================================================================================
  SymOp PrimGrid::sym_op(Index l) const {
    return SymOp::translation(coord(l, PRIM).cart());
  }
}

