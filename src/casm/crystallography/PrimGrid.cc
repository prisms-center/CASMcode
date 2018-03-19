#include "casm/crystallography/PrimGrid.hh"

#include <iostream>
#include <cmath>
#include <cassert>

#include "casm/crystallography/Lattice.hh"
#include "casm/crystallography/Coordinate.hh"
#include "casm/crystallography/UnitCellCoord.hh"

#include "casm/symmetry/SymGroupRep.hh"
#include "casm/symmetry/SymPermutation.hh"
#include "casm/symmetry/SymBasisPermute.hh"

namespace CASM {
  PrimGrid::PrimGrid(const Lattice &p_lat, const Lattice &s_lat, Index NB) {
    m_lat[PRIM] = &p_lat;
    m_lat[SCEL] = &s_lat;

    m_NB = NB;

    Eigen::Matrix3d dtrans_mat(p_lat.lat_column_mat().inverse()*s_lat.lat_column_mat());

    for(int i = 0; i < 3; i++) {
      for(int j = 0; j < 3; j++) {
        m_trans_mat(i, j) = round(dtrans_mat(i, j));
      }
    }
    if(m_trans_mat.determinant() == 0) {
      std::cerr << "CRITICAL ERROR:  Attempting to construct a PrimGrid for a superlattice that is smaller than lattice passed as its prim.\n"
                << "floating-point transformation matrix:\n" << dtrans_mat << "\n was rounded to integer matrix:\n" << m_trans_mat << "\n\n"
                << "This usage of PrimGrid is not supported. Exiting...\n";
      assert(0);
      exit(1);
    }
    /*
    std::cerr << "prim_lat is:\n";
    m_lat[PRIM]->print(std::cerr);
    std::cerr << '\n';

    std::cerr << "scel_lat is:\n";
    m_lat[SCEL]->print(std::cerr);
    std::cerr << '\n';

    std::cerr << "dtrans_mat is\n" << dtrans_mat << '\n';
    */
    matrix_type Smat, V;
    smith_normal_form(m_trans_mat, m_U, Smat, V);

    m_invU = inverse(m_U);

    m_S = Smat.diagonal();

    //std::cerr << "Smith decomposition is:\n";
    //std::cerr << m_U << "\n\n" << Smat << "\n\n" << V << "\n\n";
    m_stride[0] = m_S[0];
    m_stride[1] = m_S[0] * m_S[1];
    m_N_vol = m_S[0] * m_S[1] * m_S[2];

    m_trans_mat = m_U * Smat * V;

    Smat(0, 0) = m_N_vol / Smat(0, 0);
    Smat(1, 1) = m_N_vol / Smat(1, 1);
    Smat(2, 2) = m_N_vol / Smat(2, 2);
    m_plane_mat = inverse(V) * Smat * inverse(m_U);

    /*
    std::cerr << "trans_mat is:\n" << m_trans_mat << "\n\nplane_mat is:\n" << m_plane_mat << "\n\n";

    std::cerr << "Testing PrimGrid:\n"
            << "UnitCellCoords:\n";
    for(int i = 0; i < size(); i++)
    std::cerr << uccoord(i);

    std::cerr << "\n\nSupercell Coordinates:\n";
    for(int i = 0; i < size(); i++)
      std::cerr << coord(i, SCEL)(FRAC) << '\n';
    std::cerr << "\n";
    */
  }

  //**********************************************************************************************
  // Constructor for forcing specific choice of 'U' matrix.  Use only in very specific cases (such as applying symmetry to a PrimGrid
  PrimGrid::PrimGrid(const Lattice &p_lat,
                     const Lattice &s_lat,
                     const Eigen::Ref<const PrimGrid::matrix_type> &U,
                     const Eigen::Ref<const PrimGrid::matrix_type> &Smat,
                     Index NB) : m_U(U) {
    m_lat[PRIM] = &p_lat;
    m_lat[SCEL] = &s_lat;

    m_NB = NB;

    Eigen::Matrix3d dtrans_mat(p_lat.lat_column_mat().inverse()*s_lat.lat_column_mat());

    for(int i = 0; i < 3; i++) {
      for(int j = 0; j < 3; j++) {
        m_trans_mat(i, j) = round(dtrans_mat(i, j));
      }
    }

    if(m_trans_mat.determinant() == 0) {
      std::cerr << "CRITICAL ERROR:  Attempting to construct a PrimGrid for a superlattice that is smaller than lattice passed as its prim.\n"
                << "floating-point transformation matrix:\n" << dtrans_mat << "\n was rounded to integer matrix:\n" << m_trans_mat << "\n\n"
                << "This usage of PrimGrid is not supported. Exiting...\n";
      assert(0);
      exit(1);
    }

    m_invU = inverse(m_U);

    for(int i = 0; i < 3; i++) {
      m_S[i] = Smat(i, i);
    }

    //std::cerr << "Smith decomposition is:\n";
    m_stride[0] = m_S[0];
    m_stride[1] = m_S[0] * m_S[1];
    m_N_vol = m_S[0] * m_S[1] * m_S[2];

    m_plane_mat = adjugate(m_trans_mat);

    /*
    std::cerr << "trans_mat is:\n" << m_trans_mat << "\n\nplane_mat is:\n" << m_plane_mat << "\n\n";

    std::cerr << "Testing PrimGrid:\n"
            << "UnitCellCoords:\n";
    for(int i = 0; i < size(); i++)
    std::cerr << uccoord(i);

    std::cerr << "\n\nSupercell Coordinates:\n";
    for(int i = 0; i < size(); i++)
      std::cerr << coord(i, SCEL)(FRAC) << '\n';
    std::cerr << "\n";
    */
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
    UnitCellCoord bijk(-1, ((((m_lat[PRIM]->inv_lat_column_mat())*_coord.cart()).array() + TOL).unaryExpr(lambda)).matrix().cast<long>());

    return find(bijk);
  }

  //**********************************************************************************************

  Index PrimGrid::find(const UnitCell &_unitcell) const {
    return find(UnitCellCoord(-1, _unitcell));
  }

  //**********************************************************************************************

  Index PrimGrid::find(const UnitCellCoord &bijk) const {

    UnitCellCoord bmnp = to_canonical(get_within(bijk));

    return bmnp[1] + bmnp[2] * m_stride[0] + bmnp[3] * m_stride[1];
  }

  //**********************************************************************************************

  Index PrimGrid::find_cart(const Eigen::Ref<const Eigen::Vector3d> &_cart_coord) const {
    return find(Coordinate(_cart_coord, *m_lat[PRIM], CART));
  }

  //**********************************************************************************************

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

  UnitCellCoord PrimGrid::get_within(const UnitCellCoord &bijk)const {

    vector_type vec2 = m_plane_mat * bijk.unitcell();

    vec2[0] = ((vec2[0] % m_N_vol) + m_N_vol) % m_N_vol;
    vec2[1] = ((vec2[1] % m_N_vol) + m_N_vol) % m_N_vol;
    vec2[2] = ((vec2[2] % m_N_vol) + m_N_vol) % m_N_vol;

    return UnitCellCoord(bijk[0], (m_trans_mat * vec2) / m_N_vol);
  }

  //**********************************************************************************************

  Coordinate PrimGrid::coord(const UnitCellCoord &bijk, CELL_TYPE lat_mode)const {

    Coordinate tcoord(bijk.unitcell().cast<double>(), *(m_lat[PRIM]), FRAC);

    tcoord.set_lattice(*(m_lat[lat_mode]), CART);
    return tcoord;
  }

  //**********************************************************************************************

  Coordinate PrimGrid::coord(Index l, CELL_TYPE lat_mode)const {
    return coord(uccoord(l), lat_mode);
  }

  //**********************************************************************************************

  UnitCell PrimGrid::unitcell(Index i) const {
    return uccoord(i).unitcell();
  }

  //**********************************************************************************************

  UnitCellCoord PrimGrid::uccoord(Index i)const {
    assert(i >= 0 && i < m_N_vol && "PrimGrid::uccoord(Index i) -> index 'i' out of range!");

    UnitCellCoord bmnp(-1,
                       (i % m_stride[1]) % m_stride[0],
                       (i % m_stride[1]) / m_stride[0],
                       i / m_stride[1]);
    return from_canonical(bmnp);
  }

  //**********************************************************************************************

  SymGroupRepID PrimGrid::make_permutation_representation(const SymGroup &group, SymGroupRepID basis_permute_ID)const {

    SymGroupRepID perm_rep_ID = group.add_empty_representation();
    PrimGrid::matrix_type mat_mnp;
    UnitCellCoord mnp_shift;
    Index old_l, new_l;
    for(Index ng = 0; ng < group.size(); ng++) {
      SymBasisPermute const *rep = group[ng].get_basis_permute_rep(basis_permute_ID);
      if(!rep) {
        std::cerr << "CRITICAL ERROR: In PrimGrid::make_permutation_representation, BasisPermute representation is incorrectly initialized!\n"
                  << "                basis_permute_ID is " << basis_permute_ID << " and op index is " << group[ng].index() << "\n"
                  << "                Exiting...\n";

        exit(1);
      }

      mat_mnp = m_invU * rep->matrix() * m_U;

      std::vector<UnitCellCoord> const &b_permute = rep->data();

      Array<Index> ipermute(b_permute.size()*size());
      //std::cerr << "PRINTING b_permute array for op " << ng << ":\n";
      //begin loop over sites
      for(Index nb = 0; nb < b_permute.size(); nb++) {
        //std::cerr << b_permute.at(nb) << '\n';
        mnp_shift = to_canonical(b_permute.at(nb));

        vector_type old_mnp, new_mnp;
        EigenCounter<vector_type> mnp_count(vector_type::Zero(), m_S - vector_type::Ones(), vector_type::Ones());
        for(; mnp_count.valid(); ++mnp_count) {
          new_mnp = mat_mnp * mnp_count() + mnp_shift.unitcell();

          //map within bounds
          for(int i = 0; i < 3; i++) {
            new_mnp[i] = ((new_mnp[i] % m_S[i]) + m_S[i]) % m_S[i];
          }

          old_l = mnp_count[0] + mnp_count[1] * m_stride[0] + mnp_count[2] * m_stride[1] + nb * size();
          new_l = new_mnp[0] + new_mnp[1] * m_stride[0] + new_mnp[2] * m_stride[1] + mnp_shift[0] * size();
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
      //std::cerr << "\\end " << ng << '\n';
      group[ng].set_rep(perm_rep_ID, SymPermutation(ipermute));
    }
    return perm_rep_ID;
  }

  //**********************************************************************************************
  // for Array<Array<int> > perms=prim_grid.make_translation_permutations(NB);
  // perms[l] is the Supercell permutation that results when sites are translated by
  // prim_grid.uccoord(l).  In other words, site 'l' is translated so that (i,j,k)=(0,0,0)

  ReturnArray<Permutation> PrimGrid::make_translation_permutations(Index NB)const {
    Array<Permutation > perms;
    perms.reserve(size());
    Array<Index> ipermute(NB * size(), 0);
    UnitCellCoord shift_bmnp, bmnp;
    //std::cerr << "In PrimGrid::make_translation_permutations:\n";
    Index new_l;
    for(Index shift_l = 0; shift_l < size(); shift_l++) {
      //shift_bmnp describes translation by uccoord(shift_l)
      shift_bmnp[1] = (shift_l % m_stride[1]) % m_stride[0];
      shift_bmnp[2] = (shift_l % m_stride[1]) / m_stride[0];
      shift_bmnp[3] =   shift_l / m_stride[1];
      //std::cerr << "shift_bmnp " << shift_l << " is " << shift_bmnp;
      for(Index old_l = 0; old_l < size(); old_l++) {
        bmnp[1] = (old_l % m_stride[1]) % m_stride[0];
        bmnp[2] = (old_l % m_stride[1]) / m_stride[0];
        bmnp[3] = old_l / m_stride[1];
        //std::cerr << "old_bmnp " << old_l << " is " << bmnp;

        bmnp[1] = (bmnp[1] + shift_bmnp[1]) % m_S[0];
        bmnp[2] = (bmnp[2] + shift_bmnp[2]) % m_S[1];
        bmnp[3] = (bmnp[3] + shift_bmnp[3]) % m_S[2];

        //std::cerr << "new_bmnp " << bmnp[1] + bmnp[2] * m_stride[0] + bmnp[3] * m_stride[1] << " is " << bmnp;
        for(Index nb = 0; nb < NB; nb++) {
          new_l = bmnp[1] + bmnp[2] * m_stride[0] + bmnp[3] * m_stride[1];

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
  /// Convert UnitCellCoord (bijk) to canonical UnitCellCoord (bmnp)
  /// mnp = m_invU * ijk
  UnitCellCoord PrimGrid::to_canonical(const UnitCellCoord &bijk) const {
    UnitCellCoord bmnp(bijk[0], 0, 0, 0);
    for(int i = 0; i < 3; i++) {
      for(int j = 0; j < 3; j++) {
        bmnp[i + 1] += m_invU(i, j) * bijk[j + 1];
      }
      //Map within bounds
      bmnp[i + 1] = ((bmnp[i + 1] % m_S[i]) + m_S[i]) % m_S[i];
    }

    return bmnp;
  }

  //**********************************************************************************************
  /// Convert canonical UnitCellCoord (bmnp) to UnitCellCoord (bijk)
  /// U*mnp = ijk
  UnitCellCoord PrimGrid::from_canonical(const UnitCellCoord &bmnp) const {
    UnitCellCoord bijk(bmnp[0], 0, 0, 0);
    for(int i = 0; i < 3; i++) {
      for(int j = 0; j < 3; j++) {
        bijk[i + 1] += m_U(i, j) * bmnp[j + 1];
      }
    }
    return get_within(bijk);
  };

  //==============================================================================================
  SymOp PrimGrid::sym_op(Index l) const {
    return SymOp::translation(coord(l, PRIM).cart());
  }
}

