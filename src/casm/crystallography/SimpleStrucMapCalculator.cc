#include "casm/crystallography/SimpleStrucMapCalculator.hh"
#include "casm/crystallography/StrucMapping.hh"
#include "casm/crystallography/Coordinate.hh"

namespace CASM {
  namespace xtal {

    std::vector<Eigen::Vector3d> SimpleStrucMapCalculator::translations(MappingNode const &_node,
                                                                        SimpleStructure const &ichild_struc) const {

      SimpleStructure::Info const &p_info(info(parent()));
      SimpleStructure::Info const &c_info(info(ichild_struc));

      std::vector<Eigen::Vector3d> result;


      // Find a site in child whose occupant has the lowest number of possible
      // mapped sites in the parent. This will minimize translations
      Index i_trans(c_info.size()), n_best(p_info.size() + 2);
      for(Index i = 0; i < c_info.size() && n_best > 1; ++i) {
        auto it = max_n_species().find(c_info.names[i]);
        if(it != max_n_species().end()) {
          if(it->second < n_best) {
            n_best = it->second;
            i_trans = i;
          }
        }
        else {
          //std::cout << "\n NO TRANSLATIONS AVAILABLE FOR: " << c_info.names[i] << "\n";
          //child species not known in parent -- incompatible
          return result;
        }
      }

      std::string sp = c_info.names[i_trans];
      //std::cout << "\ni_trans: " << i_trans << "  sp: " << sp << "  c_info.size(): " << c_info.size() << "  n_best: " << n_best << "\n";

      result.reserve(n_best);

      Coordinate translation(_node.lat_node.parent.scel_lattice());
      // Try translating child atom at i_trans onto each chemically compatible site of parent
      for(Index j = 0; j < _allowed_species().size(); ++j) {
        if(_allowed_species()[j].count(sp)) {
          translation.cart() = p_info.coords.col(j) - c_info.coords.col(i_trans);
          //std::cout << "trans before within: " << translation.const_cart().transpose() << "\n";
          translation.voronoi_within();
          //std::cout << "trans after within: " << translation.const_cart().transpose() << "\n";
          result.push_back(translation.const_cart());
        }
      }

      return result;
    }

    //*******************************************************************************************

    /// \brief Creates copy of _child_struc by applying isometry, lattice transformation, translation, and site permutation of _node
    SimpleStructure SimpleStrucMapCalculator::resolve_setting(MappingNode const &_node, SimpleStructure const &_child_struc) const {
      SimpleStructure::Info const &c_info(info(_child_struc));

      SimpleStructure result(_child_struc.prefix());
      result.mol_info.resize(_node.permutation.size());
      Eigen::MatrixXd coords = _node.lat_node.stretch * _node.lat_node.isometry * c_info.coords;
      coords.colwise() += _node.basis_node.translation;
      PrimGrid const &cgrid = _node.lat_node.child;
      Index nmol = _child_struc.mol_info.size();
      for(Index i = 0; i < _node.permutation.size(); ++i) {
        Index j = _node.permutation[i];
        if(j < (nmol * cgrid.size())) {
          result.mol_info.coords.col(i) = coords.col(j % nmol) + cgrid.scel_coord(j / nmol).const_cart();
          result.mol_info.names[i] = _child_struc.mol_info.names[j % nmol];
        }
      }

      std::cerr << "SimpleStrucMapCalculator::resolve_setting() not implemented!\n";
      exit(1);
      return result;
    }
    //***************************************************************************************************

    void SimpleStrucMapCalculator::finalize(MappingNode &_node,
                                            SimpleStructure const &child_struc) const {

      populate_displacement(_node, child_struc);
      //std::cout << "**Finalizing: cost: " << _node.cost;
      _node.cost = _node.basis_weight * StrucMapping::basis_cost(_node, info(child_struc).size()) + _node.strain_weight * _node.lat_node.cost;
      //std::cout << " -> " << _node.cost << "\n";
      _node.is_valid = true;
      return;
    }

    //****************************************************************************************************************
    //            Assignment Problem methods
    //****************************************************************************************************************

    void SimpleStrucMapCalculator::populate_displacement(MappingNode &_node,
                                                         SimpleStructure const &child_struc) const {

      PrimGrid const &pgrid(_node.lat_node.parent);
      PrimGrid const &cgrid(_node.lat_node.child);
      SimpleStructure::Info const &p_info(info(parent()));
      SimpleStructure::Info const &c_info(info(child_struc));


      _node.permutation = _node.basis_node.permutation();

      // initialize displacement matrix with all zeros
      _node.displacement.setZero(3, pgrid.size()*p_info.size());

      Coordinate disp_coord(pgrid.scel_lattice());
      Index cN = c_info.size() * cgrid.size();
      // Populate displacements given as the difference in the Coordinates
      // as described by node.permutation.
      for(Index i = 0; i < _node.permutation.size(); i++) {

        // If we are dealing with a vacancy, its displacment must be zero.
        //if(node.permutation(i) >= child_struc.n_mol()) {
        //  --DO NOTHING--
        //}


        // Using min_dist routine to calculate the displacement vector that corresponds
        // to the distance used in the Cost Matrix and Hungarian Algorithm
        // The method returns the displacement vector pointing from the
        // IDEAL coordinate to the RELAXED coordinate
        if(_node.permutation[i] < cN) {

          Coordinate child_coord(c_info.coords.col(_node.permutation[i] / cgrid.size())
                                 + cgrid.scel_coord(_node.permutation[i] % cgrid.size()).const_cart()
                                 + _node.basis_node.translation,
                                 pgrid.scel_lattice(), CART);

          Coordinate parent_coord = pgrid.scel_coord(i % pgrid.size());
          parent_coord.cart() += p_info.coords.col(i / pgrid.size());
          child_coord.min_dist(parent_coord, disp_coord);
          //std::cout << "\nMap " << _node.permutation[i] << "->" << i << "\n"
          //        << "Child coord (" << _node.permutation[i] % cgrid.size()<< ", "<< _node.permutation[i] / cgrid.size() << "): "
          //        << child_coord.const_cart().transpose() << "\n"
          //        << "Parent coord: (" << i % pgrid.size()<< ", "<< i / pgrid.size() << "): "
          //        << parent_coord.const_cart().transpose() << "\n"
          //        << "Disp coord: " << disp_coord.const_cart().transpose() << "\n";
          _node.displacement.col(i) = disp_coord.const_cart();

        }
      }

      Eigen::Vector3d avg_disp = _node.displacement.rowwise().sum() / max<double>(double(cN), 1.);
      //std::cout << "avg disp: " << avg_disp << "\n";
      for(Index i = 0; i < _node.permutation.size(); i++) {
        if(_node.permutation[i] < cN) {
          _node.displacement.col(i) -= avg_disp;
        }
      }
      _node.basis_node.translation -= avg_disp;
      // End of filling displacements
    }

    /*
     * Finding the cost_matrix given the relaxed structure
     * This will always return a square matrix with the extra elements
     * reflecting the vacancies specified in the ideal supercell.
     * Costs are calculated in context of the lattice.
     * cost_matrix(i,j) is cost of mapping child site 'j' onto parent site 'i'
     */
    //****************************************************************************************************************
    bool SimpleStrucMapCalculator::populate_cost_mat(MappingNode &_node,
                                                     SimpleStructure const &child_struc) const {
      PrimGrid const &pgrid(_node.lat_node.parent);
      PrimGrid const &cgrid(_node.lat_node.child);
      Eigen::Vector3d const &translation(_node.basis_node.translation);
      Eigen::MatrixXd &cost_matrix(_node.basis_node.cost_mat);
      Eigen::Matrix3d metric = (_node.lat_node.stretch * _node.lat_node.stretch).inverse();

      SimpleStructure::Info const &p_info(info(parent()));
      SimpleStructure::Info const &c_info(info(child_struc));

      Index pN = p_info.size() * pgrid.size();
      Index cN = c_info.size() * cgrid.size();

      cost_matrix = Eigen::MatrixXd::Constant(pN, pN, StrucMapping::small_inf());

      if(pN < cN) {
        //std::cout << "Insufficient number of parent sites to accept child atoms\n";
        return false;
      }
      Index residual = pN - cN;
      if(residual > m_va_allowed.size()*pgrid.size()) {
        //std::cout << "Parent cannot accommodate enough vacancies to make mapping possible.\n";
        return false;
      }

      // index of basis atom in child_struc
      Index bc = 0;
      // index of atom in child supercell
      Index ac = 0;

      // loop through all the sites of the child structure
      // As always, each sublattice is traversed contiguously
      for(; bc < c_info.size(); bc++) {
        auto it = max_n_species().find(c_info.names[bc]);
        if(it == max_n_species().end() || (it->second * pgrid.size()) < cgrid.size()) {
          //std::cout << "Parent structure cannot accommodate enough atoms of type '" << c_info.names[bc] << "'.\n"
          //          << "Must accommodate at least " <<cgrid.size()<<" but only accommodates " << (it->second*pgrid.size()) << "\n";
          return false;
        }

        // For each sublattice, loop over lattice points, 'n'. 'ac' tracks linear index of atoms in child supercell
        for(Index lc = 0; lc < cgrid.size(); ++lc, ++ac) {
          Coordinate child_coord(c_info.coords.col(bc) + cgrid.scel_coord(lc).const_cart() + translation, pgrid.scel_lattice(), CART);
          // loop through all the sites in the parent supercell
          Index ap = 0;
          for(Index bp = 0; bp < p_info.size(); ++bp) {
            if(!_allowed_species()[bp].count(c_info.names[bc])) {
              ap += pgrid.size();
              continue;
            }
            Coordinate parent_coord(p_info.coords.col(bp), pgrid.scel_lattice(), CART);

            for(Index lp = 0; lp < pgrid.size(); ++lp, ++ap) {
              cost_matrix(ap, ac) = (parent_coord + pgrid.scel_coord(lp)).min_dist2(child_coord, metric);
            }
          }
        }
      }

      // If there are unvisited columns of cost_mat (because fewer sites in the child structure than in parent),
      // we will treat them as a vacant species in the child struc and set them to zero cost for mapping onto
      // parent sites that allow vacancies

      if(residual) {
        // loop over parent sublats that allow Va, set corresponding block to zero cost
        for(Index bp : m_va_allowed) {
          cost_matrix.block(bp * pgrid.size(), cN - 1, pgrid.size(), residual).setZero();
        }
      }


      // JCT: I'm not sure if there's an easy way to check if the cost matrix is viable in all cases
      //      Some of the simpler checks I could think of failed for edge cases with vacancies.
      //      If we return an invalid cost matrix, the Hungarian routines will detect that it is invalid,
      //      so maybe there's no point in doing additional checks here.
      return true;
    }

  }
}
