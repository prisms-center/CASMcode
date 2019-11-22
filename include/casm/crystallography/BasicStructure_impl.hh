#ifndef CASM_BasicStructure_impl
#define CASM_BasicStructure_impl

#include <sstream>
#include <iostream>
#include <string>
#include <sys/stat.h>
#include <sys/types.h>
#include <vector>

#include <boost/filesystem.hpp>
#include <boost/filesystem/fstream.hpp>
#include "casm/crystallography/BasicStructure.hh"

#include "casm/misc/CASM_Eigen_math.hh"
#include "casm/crystallography/LatticeIsEquivalent.hh"
#include "casm/crystallography/Coordinate.hh"
#include "casm/crystallography/UnitCellCoord.hh"
#include "casm/crystallography/PrimGrid.hh"
#include "casm/crystallography/Molecule.hh"
#include "casm/basis_set/DoFIsEquivalent.hh"
#include "casm/basis_set/OccupationDoFTraits.hh"
#include "casm/symmetry/SymPermutation.hh"
#include "casm/symmetry/SymBasisPermute.hh"
#include "casm/symmetry/SymGroupRep.hh"
#include "casm/crystallography/Niggli.hh"

namespace CASM {
  namespace xtal {
    template<typename CoordType>
    BasicStructure<CoordType>::BasicStructure(const fs::path &filepath) : m_lattice() {
      if(!fs::exists(filepath)) {
        std::cerr << "Error in BasicStructure<CoordType>::BasicStructure<CoordType>(const fs::path &filepath)." << std::endl;
        std::cerr << "  File does not exist at: " << filepath << std::endl;
        exit(1);
      }
      fs::ifstream infile(filepath);
      read(infile);
    }

    //***********************************************************

    template<typename CoordType>
    BasicStructure<CoordType>::BasicStructure(const BasicStructure &RHS) :
      m_lattice(RHS.lattice()),
      m_title(RHS.title()),
      m_basis(RHS.basis()),
      m_dof_map(RHS.m_dof_map) {
      for(Index i = 0; i < basis().size(); i++) {
        m_basis[i].set_lattice(lattice(), CART);
      }
    }

    //***********************************************************

    template<typename CoordType>
    BasicStructure<CoordType>::~BasicStructure() {}

    //***********************************************************

    template<typename CoordType>
    BasicStructure<CoordType> &BasicStructure<CoordType>::operator=(const BasicStructure<CoordType> &RHS) {
      m_lattice = RHS.lattice();
      m_title = RHS.title();
      set_basis(RHS.basis());
      m_dof_map = RHS.m_dof_map;

      for(Index i = 0; i < basis().size(); i++)
        m_basis[i].set_lattice(lattice(), CART);

      return *this;
    }

    //************************************************************

    template<typename CoordType>
    DoFSet const &BasicStructure<CoordType>::global_dof(std::string const &_dof_type) const {
      auto it = m_dof_map.find(_dof_type);
      if(it != m_dof_map.end())
        return (it->second);
      else
        throw std::runtime_error(std::string("In BasicStructure::dof(), this structure does not contain any global DoF's of type " + _dof_type));

    }



    //***********************************************************

    template<typename CoordType>
    void BasicStructure<CoordType>::copy_attributes_from(const BasicStructure<CoordType> &RHS) {

    }

    //***********************************************************
    /*
      template<typename CoordType>
      BasicStructure<CoordType> &BasicStructure<CoordType>::apply_sym(const SymOp &op) {
      for(Index i = 0; i < basis().size(); i++) {
      m_basis[i].apply_sym(op);
      }
      return *this;
      }
    */
    //***********************************************************

    template<typename CoordType>
    void BasicStructure<CoordType>::reset() {
      set_site_internals();

      within();
    }

    //***********************************************************

    template<typename CoordType>
    void BasicStructure<CoordType>::update() {
      set_site_internals();
    }

    //***********************************************************

    template<typename CoordType>
    void BasicStructure<CoordType>::within() {
      for(Index i = 0; i < basis().size(); i++) {
        m_basis[i].within();
      }
      return;
    }


    //***********************************************************

    template<typename CoordType>
    void BasicStructure<CoordType>::_generate_factor_group_slow(SymGroup &factor_group, SymGroup const &super_group, bool time_reversal_enabled) const {
      std::vector<CoordType> trans_basis;
      Index pg, b0, b1, b2;
      Coordinate t_tau(lattice());
      Index num_suc_maps;

      SymGroup const &point_group = super_group;
      Index time_max = Index(_time_reversal_active() && time_reversal_enabled);
      //std::cout << "time_max is " << time_max << "\n";
      //std::cout << "time_reversal_enabled is " << time_reversal_enabled << "\n";
      CASM::SymOp test_op;
      if(factor_group.size() != 0) {
        std::cerr << "WARNING in BasicStructure<CoordType>::generate_factor_group_slow" << std::endl;
        std::cerr << "The factor group passed isn't empty and it's about to be rewritten!" << std::endl;
        factor_group.clear();
      }
      factor_group.set_lattice(lattice());
      //Loop over all point group ops of the lattice
      for(pg = 0; pg < point_group.size(); pg++) {
        test_op = point_group[pg];
        for(Index tr = 0; tr <= time_max; tr++) {
          if(tr) {
            test_op = test_op * CASM::SymOp::time_reversal_op();
          }

          if(!_is_lattice_pg_op(test_op)) {
            continue;
          }

          trans_basis.clear();
          //First, generate the symmetrically transformed basis sites
          //Loop over all sites in basis
          for(b0 = 0; b0 < basis().size(); b0++) {
            trans_basis.push_back(test_op * m_basis[b0]);
          }

          //Using the symmetrically transformed basis, find all possible translations
          //that MIGHT map the symmetrically transformed basis onto the original basis
          for(b0 = 0; b0 < trans_basis.size(); b0++) {

            if(!m_basis[0].compare_type(trans_basis[b0])) {
              continue;
            }

            t_tau = m_basis[0] - trans_basis[b0];

            t_tau.within();
            num_suc_maps = 0; //Keeps track of number of old->new basis site mappings that are found

            double tdist = 0.0;
            double max_error = 0.0;
            std::vector<Index> mappings(basis().size(), basis().size());
            for(b1 = 0; b1 < basis().size(); b1++) { //Loop over original basis sites
              for(b2 = 0; b2 < trans_basis.size(); b2++) { //Loop over symmetrically transformed basis sites

                //see if translation successfully maps the two sites uniquely

                if(basis()[b1].compare(trans_basis[b2], t_tau)) {// maybe nuke this with Hungarian
                  tdist = basis()[b1].min_dist(Coordinate(trans_basis[b2]) + t_tau);
                  if(tdist > max_error) {
                    max_error = tdist;
                  }
                  mappings[b1] = b2;
                  num_suc_maps++;
                  break;
                }
              }

              //break out of outer loop if inner loop finds no successful map
              if(b2 == trans_basis.size()) {
                break;
              }

            }
            std::set<Index> unique_mappings;
            for(auto &e : mappings)
              unique_mappings.insert(e);
            if(num_suc_maps == basis().size() && unique_mappings.size() == mappings.size()) {
              //If all atoms in the basis are mapped successfully, try to add the corresponding
              //symmetry operation to the factor_group
              Coordinate center_of_mass(lattice());
              for(Index b = 0; b < basis().size(); b++) {
                //for each basis site loop through all trans_basis to find the closest one
                Coordinate tshift(lattice());
                //in tshift is stored trans_basis - basis
                tshift.cart() *= (1.0 / basis().size());
                center_of_mass += tshift;
              }

              /*
                test_op operates on all of basis and save = trans_basis
                //t_shift is the vector from basis to op_basis magnitude of tshift=min_dist
                average all t_shifts and add/subtract to t_tau

                if t_shift is b-> ob then subtract
                if t_shift is ob -> b then add
                matrix * basis + tau = operbasis
              */
              t_tau -= center_of_mass;/**/

              CASM::SymOp tSym(CASM::SymOp::translation(t_tau.cart())*test_op);
              tSym.set_map_error(max_error);

              if(!factor_group.contains_periodic(tSym)) {
                factor_group.push_back(tSym);
              }
            }
          }
        }
      } //End loop over point_group operations
      factor_group.enforce_group(lattice().tol());
      factor_group.sort();
      factor_group.max_error();

      return;
    }

    //************************************************************
    //***********************************************************

    template<typename CoordType>
    void BasicStructure<CoordType>::generate_factor_group_slow(SymGroup &factor_group) const {
      SymGroup point_group((SymGroup::lattice_point_group(lattice())));

      _generate_factor_group_slow(factor_group, point_group);
      return;
    }

    //************************************************************

    template<typename CoordType>
    void BasicStructure<CoordType>::generate_factor_group(SymGroup &factor_group) const {
      BasicStructure<CoordType> tprim(lattice());
      factor_group.clear();
      factor_group.set_lattice(lattice());
      // CASE 1: Structure is primitive
      if(is_primitive(tprim)) {
        generate_factor_group_slow(factor_group);
        return;
      }


      // CASE 2: Structure is not primitive

      PrimGrid prim_grid(tprim.lattice(), lattice());
      SymGroup prim_fg;
      tprim.generate_factor_group_slow(prim_fg);

      SymGroup point_group((SymGroup::lattice_point_group(lattice())));

      for(Index i = 0; i < prim_fg.size(); i++) {
        if(point_group.find_no_trans(prim_fg[i]) == point_group.size()) {
          continue;
        }
        else {
          for(Index j = 0; j < prim_grid.size(); j++) {
            factor_group.push_back(CASM::SymOp::translation(prim_grid.scel_coord(j).cart())*prim_fg[i]);
            // set lattice, in case CASM::SymOp::operator* ever changes
          }
        }
      }

      return;
    }

    //***********************************************************
    /**
     * It is NOT wise to use this function unless you have already
     * initialized a superstructure with lattice vectors.
     *
     * It is more wise to use the two methods that call this method:
     * Either the overloaded * operator which does:
     *  SCEL_Lattice * Prim_Structrue = New_Superstructure
     *       --- or ---
     *  New_Superstructure=Prim_BasicStructure<CoordType>.create_superstruc(SCEL_Lattice);
     *
     *  Both of these will return NEW superstructures.
     */
    //***********************************************************

    template<typename CoordType>
    void BasicStructure<CoordType>::fill_supercell(const BasicStructure<CoordType> &prim) {
      Index i, j;

      copy_attributes_from(prim);

      PrimGrid prim_grid(prim.lattice(), lattice());

      m_basis.clear();

      //loop over basis sites of prim
      for(j = 0; j < prim.basis().size(); j++) {

        //loop over prim_grid points
        for(i = 0; i < prim_grid.size(); i++) {

          //push back translated basis site of prim onto superstructure basis
          push_back(prim.basis()[j] + prim_grid.scel_coord(i));

          m_basis.back().within();
        }
      }

      return;
    }

    //***********************************************************
    /**
     * Operates on the primitive structure and takes as an argument
     * the supercell lattice.  It then returns a new superstructure.
     *
     * This is similar to the Lattice*Primitive routine which returns a
     * new superstructure.  Unlike the fill_supercell routine which takes
     * the primitive structure, this WILL fill the sites.
     */
    //***********************************************************

    template<typename CoordType>
    BasicStructure<CoordType> BasicStructure<CoordType>::create_superstruc(const Lattice &scel_lat) const {
      BasicStructure<CoordType> tsuper(scel_lat);
      tsuper.fill_supercell(*this);
      return tsuper;
    }


    //***********************************************************
    /**
     * Determines if structure is primitive description of the crystal
     */
    //***********************************************************

    template<typename CoordType>
    bool BasicStructure<CoordType>::is_primitive() const {
      SymGroup valid_translations, identity_group;
      identity_group.push_back(CASM::SymOp());
      _generate_factor_group_slow(valid_translations, identity_group, false);
      if(valid_translations.size() == 1)
        return true;

      return false;
    }


    //***********************************************************
    /**
     * Determines if structure is primitive description of the crystal
     * If not, finds primitive cell and copies to new_prim
     */
    //***********************************************************

    template<typename CoordType>
    bool BasicStructure<CoordType>::is_primitive(BasicStructure<CoordType> &new_prim) const {
      Eigen::Vector3d prim_vec0(lattice()[0]), prim_vec1(lattice()[1]), prim_vec2(lattice()[2]);
      std::vector<Eigen::Vector3d > shift;
      double tvol, min_vol;
      bool prim_flag = true;
      double prim_vol_tol = std::abs(0.5 * lattice().vol() / double(basis().size())); //sets a hard lower bound for the minimum value of the volume of the primitive cell

      SymGroup valid_translations, identity_group;
      identity_group.push_back(CASM::SymOp());
      _generate_factor_group_slow(valid_translations, identity_group, false);
      if(valid_translations.size() > 1) {
        prim_flag = false;
        for(auto &trans : valid_translations) {
          shift.push_back(trans.tau());
        }
      }



      if(prim_flag) {
        new_prim = *this;
        return true;
      }


      shift.push_back(lattice()[0]);
      shift.push_back(lattice()[1]);
      shift.push_back(lattice()[2]);

      //We want to minimize the volume of the primitivized cell, but to make it not a weird shape
      //that leads to noise we also minimize the dot products like reduced cell would
      min_vol = std::abs(lattice().vol());
      for(Index sh = 0; sh < shift.size(); sh++) {
        for(Index sh1 = sh + 1; sh1 < shift.size(); sh1++) {
          for(Index sh2 = sh1 + 1; sh2 < shift.size(); sh2++) {
            tvol = std::abs(triple_prod(shift[sh], shift[sh1], shift[sh2]));
            if(tvol < min_vol && tvol > prim_vol_tol) {
              min_vol = tvol;
              prim_vec0 = shift[sh];
              prim_vec1 = shift[sh1];
              prim_vec2 = shift[sh2];
            }
          }
        }

      }


      Lattice new_lat(prim_vec0, prim_vec1, prim_vec2);
      Lattice reduced_new_lat = niggli(new_lat, lattice().tol());
      //The lattice so far is OK, but it's noisy enough to matter for large
      //superstructures. We eliminate the noise by reconstructing it now via
      //rounded to integer transformation matrix.

      Eigen::Matrix3d transmat, invtransmat;
      //lattice().lat_column_mat() = reduced_new_lat.lat_column_mat()*transmat
      transmat = reduced_new_lat.inv_lat_column_mat() * lattice().lat_column_mat();

      invtransmat = iround(transmat).cast<double>().inverse();

      //By using invtransmat, the new prim is guaranteed to perfectly tile the old prim
      Lattice reconstructed_reduced_new_lat(lattice().lat_column_mat() * invtransmat, lattice().tol());


      new_prim.set_lattice(reconstructed_reduced_new_lat, CART);
      CoordType tsite(new_prim.lattice());
      for(Index nb = 0; nb < basis().size(); nb++) {
        tsite = m_basis[nb];
        tsite.set_lattice(new_prim.lattice(), CART);
        if(new_prim.find(tsite) == new_prim.basis().size()) {
          tsite.within();
          new_prim.push_back(tsite);
        }
      }

      return false;
    }

    //***********************************************************

    template<typename CoordType>
    void BasicStructure<CoordType>::set_site_internals() {
      for(Index nb = 0; nb < basis().size(); nb++) {
        m_basis[nb].set_basis_ind(nb);
      }
    }

    //***********************************************************

    template<typename CoordType> template<typename CoordType2>
    Index BasicStructure<CoordType>::find(const CoordType2 &test_site) const {
      for(Index i = 0; i < basis().size(); i++) {
        if(m_basis[i].compare(test_site)) {
          return i;
        }
      }
      return basis().size();
    }

    //***********************************************************

    template<typename CoordType> template<typename CoordType2>
    Index BasicStructure<CoordType>::find(const CoordType2 &test_site, const Coordinate &shift) const {
      for(Index i = 0; i < basis().size(); i++) {
        if(m_basis[i].compare(test_site, shift)) {
          return i;
        }
      }
      return basis().size();
    }

    //***********************************************************

    template<typename CoordType>
    void BasicStructure<CoordType>::set_lattice(const Lattice &new_lat, COORD_TYPE mode) {

      m_lattice = new_lat;

      for(Index nb = 0; nb < basis().size(); nb++) {
        m_basis[nb].set_lattice(lattice(), mode);
      }
    }

    //***********************************************************


    template<typename CoordType>
    void BasicStructure<CoordType>::set_title(std::string const &_title) {
      m_title = _title;
    }

    //\Liz D 032514
    //***********************************************************
    /**
     * Allows for the basis elements of a basic structure to be
     * manually set, e.g. as in jsonParser.cc.
     */
    //***********************************************************


    template<typename CoordType>
    void BasicStructure<CoordType>::set_basis(std::vector<CoordType> const &_basis, COORD_TYPE mode) {
      m_basis.clear();
      m_basis.reserve(_basis.size());
      for(CoordType const &site : _basis)
        push_back(site, mode);

    }

    template<typename CoordType>
    void BasicStructure<CoordType>::clear_basis() {
      m_basis.clear();
      this->reset();

    }

    template<typename CoordType>
    void BasicStructure<CoordType>::set_occ(Index basis_ind, int _val) {
      m_basis[basis_ind].set_occ_value(_val);
    }

    template<typename CoordType>
    void BasicStructure<CoordType>::push_back(CoordType const &_site, COORD_TYPE mode) {
      m_basis.push_back(_site);
      m_basis.back().set_basis_ind(basis().size() - 1);
      m_basis.back().set_lattice(lattice(), mode);
    }

    //************************************************************
    /// Counts sites that allow vacancies
    template<typename CoordType>
    Index BasicStructure<CoordType>::max_possible_vacancies()const {
      Index result(0);
      for(Index i = 0; i < basis().size(); i++) {
        if(m_basis[i].contains("Va"))
          ++result;
      }
      return result;
    }

    //************************************************************
    //read a POSCAR like file and collect all the structure variables
    //modified to read PRIM file and determine which basis to use
    //Changed by Ivy to read new VASP POSCAR format

    template<typename CoordType>
    void BasicStructure<CoordType>::read(std::istream &stream) {
      int i, t_int;
      char ch;
      std::vector<double> num_elem;
      std::vector<std::string> elem_array;
      bool read_elem = false;
      std::string tstr;
      std::stringstream tstrstream;

      CoordType tsite(lattice());

      bool SD_flag = false;
      getline(stream, m_title);
      if(title().back() == '\r')
        throw std::runtime_error(std::string("Structure file is formatted for DOS. Please convert to Unix format. (This can be done with the dos2unix command.)"));

      m_lattice.read(stream);

      stream.ignore(100, '\n');

      //Search for Element Names
      ch = stream.peek();
      while(ch != '\n' && !stream.eof()) {
        if(isalpha(ch)) {
          read_elem = true;
          stream >> tstr;
          elem_array.push_back(tstr);
          ch = stream.peek();
        }
        else if(ch == ' ' || ch == '\t') {
          stream.ignore();
          ch = stream.peek();
        }
        else if(ch >= '0' && ch <= '9') {
          break;
        }
        else {
          throw std::runtime_error(std::string("Error attempting to read Structure. Error reading atom names."));
        }
      }

      if(read_elem == true) {
        stream.ignore(10, '\n');
        ch = stream.peek();
      }

      //Figure out how many species
      int num_sites = 0;
      while(ch != '\n' && !stream.eof()) {
        if(ch >= '0' && ch <= '9') {
          stream >> t_int;
          num_elem.push_back(t_int);
          num_sites += t_int;
          ch = stream.peek();
        }
        else if(ch == ' ' || ch == '\t') {
          stream.ignore();
          ch = stream.peek();
        }
        else {
          throw std::runtime_error(std::string("Error in line 6 of structure input file. Line 6 of structure input file should contain the number of sites."));
        }
      }
      stream.get(ch);

      // fractional coordinates or cartesian
      COORD_MODE input_mode(FRAC);

      stream.get(ch);
      while(ch == ' ' || ch == '\t') {
        stream.get(ch);
      }

      if(ch == 'S' || ch == 's') {
        SD_flag = true;
        stream.ignore(1000, '\n');
        while(ch == ' ' || ch == '\t') {
          stream.get(ch);
        }
        stream.get(ch);
      }

      if(ch == 'D' || ch == 'd') {
        input_mode.set(FRAC);
      }
      else if(ch == 'C' || ch == 'c') {
        input_mode.set(CART);
      }
      else if(!SD_flag) {
        throw std::runtime_error(std::string("Error in line 7 of structure input file. Line 7 of structure input file should specify Direct, Cartesian, or Selective Dynamics."));
      }
      else if(SD_flag) {
        throw std::runtime_error(std::string("Error in line 8 of structure input file. Line 8 of structure input file should specify Direct or Cartesian when Selective Dynamics is on."));
      }

      stream.ignore(1000, '\n');
      //Clear basis if it is not empty
      if(basis().size() != 0) {
        std::cerr << "The structure is going to be overwritten." << std::endl;
        m_basis.clear();
      }

      if(read_elem) {
        int j = -1;
        int sum_elem = 0;
        m_basis.reserve(num_sites);
        for(i = 0; i < num_sites; i++) {
          if(i == sum_elem) {
            j++;
            sum_elem += num_elem[j];
          }

          tsite.read(stream, elem_array[j], SD_flag);
          push_back(tsite);
        }
      }
      else {
        //read the site info
        m_basis.reserve(num_sites);
        for(i = 0; i < num_sites; i++) {
          tsite.read(stream, SD_flag);
          if((stream.rdstate() & std::ifstream::failbit) != 0) {
            std::cerr << "Error reading site " << i + 1 << " from structure input file." << std::endl;
            exit(1);
          }
          push_back(tsite);
        }
      }

      // Check whether there are additional sites listed in the input file
      std::string s;
      getline(stream, s);
      std::istringstream tmp_stream(s);
      Eigen::Vector3d coord;
      tmp_stream >> coord;
      if(tmp_stream.good()) {
        throw std::runtime_error(std::string("ERROR: too many sites listed in structure input file."));
      }

      update();
      return;

    }

    //***********************************************************

    template<typename CoordType>
    void BasicStructure<CoordType>::print_xyz(std::ostream &stream, bool frac) const {
      stream << basis().size() << '\n';
      stream << title() << '\n';
      stream.precision(7);
      stream.width(11);
      stream.flags(std::ios::showpoint | std::ios::fixed | std::ios::right);
      stream << "      a       b       c" << '\n';
      stream << lattice().lat_column_mat() << '\n';
      for(Index i = 0; i < basis().size(); i++) {
        stream << std::setw(2) << basis()[i].occ_name() << " ";
        if(frac) {
          stream << std::setw(12) << basis()[i].frac().transpose() << '\n';
        }
        else {
          stream << std::setw(12) << basis()[i].cart() << '\n';
        }
      }

    }

    //***********************************************************
    /*
      template<typename CoordType>
      void BasicStructure<CoordType>::print_cif(std::ostream &stream) const {
      const char quote = '\'';
      const char indent[] = "   ";

      //double amag, bmag, cmag;
      //double alpha, beta, gamma;

      // Copying format based on VESTA .cif output.

      // Heading text.

      stream << '#';
      for(int i = 0; i < 70; i++) {
      stream << '=';
      }
      stream << "\n\n";
      stream << "# CRYSTAL DATA\n\n";
      stream << '#';
      for(int i = 0; i < 70; i++) {
      stream << '-';
      }
      stream << "\n\n";
      stream << "data_CASM\n\n\n";

      stream.precision(5);
      stream.width(11);
      stream.flags(std::ios::showpoint | std::ios::fixed | std::ios::left);

      stream << std::setw(40) << "_pd_phase_name" << quote << title << quote << '\n';
      stream << std::setw(40) << "_cell_length_a" << lattice().lengths[0] << '\n';
      stream << std::setw(40) << "_cell_length_b" << lattice().lengths[1] << '\n';
      stream << std::setw(40) << "_cell_length_c" << lattice().lengths[2] << '\n';
      stream << std::setw(40) << "_cell_angle_alpha" << lattice().angles[0] << '\n';
      stream << std::setw(40) << "_cell_angle_beta" << lattice().angles[1] << '\n';
      stream << std::setw(40) << "_cell_angle_gamma" << lattice().angles[2] << '\n';
      stream << std::setw(40) << "_symmetry_space_group_name_H-M" << quote << "TBD" << quote << '\n';
      stream << std::setw(40) << "_symmetry_Int_Tables_number" << "TBD" << "\n\n";

      stream << "loop_\n";
      stream << "_symmetry_equiv_pos_as_xyz\n";

      // Equivalent atom positions here. Form: 'x, y, z', '-x, -y, -z', 'x+1/2, y+1/2, z', etc.
      // Use stream << indent << etc.

      stream << '\n';
      stream << "loop_\n";
      stream << indent << "_atom_site_label" << '\n';
      stream << indent << "_atom_site_occupancy" << '\n';
      stream << indent << "_atom_site_fract_x" << '\n';
      stream << indent << "_atom_site_fract_y" << '\n';
      stream << indent << "_atom_site_fract_z" << '\n';
      stream << indent << "_atom_site_adp_type" << '\n';
      stream << indent << "_atom_site_B_iso_or_equiv" << '\n';
      stream << indent << "_atom_site_type_symbol" << '\n';

      // Use stream << indent << etc.
      }
    */
    //***********************************************************

    template<typename CoordType>
    BasicStructure<CoordType> &BasicStructure<CoordType>::operator+=(const Coordinate &shift) {

      for(Index i = 0; i < basis().size(); i++) {
        m_basis[i] += shift;
      }

      return (*this);
    }


    //***********************************************************

    template<typename CoordType>
    BasicStructure<CoordType> &BasicStructure<CoordType>::operator-=(const Coordinate &shift) {

      for(Index i = 0; i < basis().size(); i++) {
        m_basis[i] -= shift;
      }
      //factor_group -= shift;
      //asym_unit -= shift;
      return (*this);
    }

    //***********************************************************

    template<typename CoordType>
    BasicStructure<CoordType> operator*(const Lattice &LHS, const BasicStructure<CoordType> &RHS) {
      BasicStructure<CoordType> tsuper(LHS);
      tsuper.fill_supercell(RHS);
      return tsuper;
    }

    //***********************************************************

    /// \brief Returns true if @param _op leaves lattice and global DoFs (if any) invariant
    template<typename CoordType>
    bool BasicStructure<CoordType>::_is_lattice_pg_op(CASM::SymOp const &evil_op) const {
      //TODO: Resolve where the Adapter should be called. I'm doing it here because
      //DoFIsEquivalent is outside of crystallography, and we decided that crystallography/Adapter.hh
      //should be limited to its own module.
      auto _op = adapter::Adapter<xtal::SymOp, CASM::SymOp>()(evil_op);
      //std::cout << "CHECKING OP: \n" << _op.matrix() << std::endl;
      if(!LatticeIsEquivalent(lattice())(_op)) {
        //std::cout << "FAILED LATTICE CHECK\n"<< std::endl;
        return false;

      }
      for(auto const &dof : m_dof_map) {
        if(!DoFIsEquivalent(dof.second)(_op)) {
          //std::cout << "FAILED DOF CHECK\n" << std::endl;
          return false;
        }
      }
      //std::cout << "PASSED\n"<< std::endl;
      return true;

    }

    //****************************************************

    /// \brief Returns true if structure has attributes affected by time reversal
    // private for now, expose if necessary
    template<typename CoordType>
    bool BasicStructure<CoordType>::_time_reversal_active() const {
      for(auto const &dof : m_dof_map)
        if(dof.second.traits().time_reversal_active())
          return true;
      for(CoordType const &site : basis())
        if(site.time_reversal_active())
          return true;
      return false;
    }

    //***********************************************************

    template<typename CoordType>
    std::vector<UnitCellCoord> symop_site_map(CASM::SymOp const &_op, BasicStructure<CoordType> const &_struc) {
      return symop_site_map(_op, _struc, _struc.lattice().tol());
    }

    //***********************************************************

    template<typename CoordType>
    std::vector<UnitCellCoord> symop_site_map(CASM::SymOp const &_op, BasicStructure<CoordType> const &_struc, double _tol) {
      std::vector<UnitCellCoord> result;
      // Determine how basis sites transform from the origin unit cell
      for(int b = 0; b < _struc.basis().size(); b++) {
        result.push_back(UnitCellCoord(_struc, CASM::copy_apply(_op, _struc.basis()[b]), _tol));
      }
      return result;
    }

    //************************************************************

    /// Returns an std::vector of each *possible* Specie in this Structure
    template<typename CoordType>
    std::vector<std::string> struc_species(BasicStructure<CoordType> const &_struc) {

      std::vector<Molecule> tstruc_molecule = struc_molecule(_struc);
      std::set<std::string> result;

      Index i, j;

      //For each molecule type
      for(i = 0; i < tstruc_molecule.size(); i++) {
        // For each atomposition in the molecule
        for(j = 0; j < tstruc_molecule[i].size(); j++)
          result.insert(tstruc_molecule[i].atom(j).name());
      }

      return std::vector<std::string>(result.begin(), result.end());
    }

    //************************************************************

    /// Returns an std::vector of each *possible* Molecule in this Structure
    template<typename CoordType>
    std::vector<Molecule> struc_molecule(BasicStructure<CoordType> const &_struc) {

      std::vector<Molecule> tstruc_molecule;
      Index i, j;

      //loop over all Sites in basis
      for(i = 0; i < _struc.basis().size(); i++) {
        //loop over all Molecules in Site
        for(j = 0; j < _struc.basis()[i].occupant_dof().size(); j++) {
          //Collect unique Molecules
          if(!contains(tstruc_molecule, _struc.basis()[i].occupant_dof()[j])) {
            tstruc_molecule.push_back(_struc.basis()[i].occupant_dof()[j]);
          }
        }
      }//end loop over all Sites

      return tstruc_molecule;
    }

    //************************************************************
    /// Returns an std::vector of each *possible* Molecule in this Structure
    template<typename CoordType>
    std::vector<std::string> struc_molecule_name(BasicStructure<CoordType> const &_struc) {

      // get Molecule allowed in struc
      std::vector<Molecule> struc_mol = struc_molecule(_struc);

      // store Molecule names in vector
      std::vector<std::string> struc_mol_name;
      for(int i = 0; i < struc_mol.size(); i++) {
        struc_mol_name.push_back(struc_mol[i].name());
      }

      return struc_mol_name;
    }

    //************************************************************
    /// Returns an std::vector of each *possible* Molecule in this Structure
    template<typename CoordType>
    std::vector<std::vector<std::string> > allowed_molecule_unique_names(BasicStructure<CoordType> const &_struc) {
      using IPair = std::pair<Index, Index>;
      std::map<std::string, std::vector<Molecule> > name_map;
      std::map<std::string, IPair> imap;

      std::vector<std::vector<std::string> > result(_struc.basis().size());
      for(Index b = 0; b < _struc.basis().size(); ++b) {
        for(Index j = 0; j < _struc.basis(b).occupant_dof().size(); ++j) {
          Molecule const &mol(_struc.basis(b).occupant_dof()[j]);
          result[b].push_back(mol.name());
          auto it = name_map.find(mol.name());
          if(it == name_map.end()) {
            name_map[mol.name()].push_back(mol);
            imap[mol.name()] = {b, j};
          }
          else {
            Index i = find_index(it->second, mol);
            if(i == it->second.size()) {
              it->second.push_back(mol);
              if(i == 1) {
                auto inds = imap[mol.name()];
                result[inds.first][inds.second] += ".1";
              }
            }
            if(i > 0)
              result[b][j] += ("." + std::to_string(i + 1));
          }
        }
      }
      return result;
    }

    //************************************************************
    /// Returns a vector with a list of allowed molecule names at each site
    template<typename CoordType>
    std::vector<std::vector<std::string> > allowed_molecule_names(BasicStructure<CoordType> const &_struc) {
      std::vector<std::vector<std::string> > result(_struc.basis().size());

      for(Index b = 0; b < _struc.basis().size(); ++b)
        result[b] = _struc.basis(b).allowed_occupants();

      return result;
    }

    //************************************************************

    template<typename CoordType>
    std::vector<DoFKey> continuous_local_dof_types(BasicStructure<CoordType> const &_struc) {
      std::set<std::string> tresult;

      for(CoordType const &site : _struc.basis()) {
        auto sitetypes = site.dof_types();
        tresult.insert(sitetypes.begin(), sitetypes.end());
      }
      return std::vector<std::string>(tresult.begin(), tresult.end());
    }

    //************************************************************

    template<typename CoordType>
    std::vector<DoFKey> all_local_dof_types(BasicStructure<CoordType> const &_struc) {
      std::set<std::string> tresult;

      for(CoordType const &site : _struc.basis()) {
        auto sitetypes = site.dof_types();
        tresult.insert(sitetypes.begin(), sitetypes.end());
        if(site.occupant_dof().size() > 1) {
          tresult.insert(DoFType::occupation().name());
        }
      }
      return std::vector<std::string>(tresult.begin(), tresult.end());
    }


    //************************************************************

    template<typename CoordType>
    std::vector<DoFKey> global_dof_types(BasicStructure<CoordType> const &_struc) {
      std::vector<std::string> result;
      for(auto const &dof :  _struc.global_dofs())
        result.push_back(dof.first);
      return result;
    }


    //************************************************************

    template<typename CoordType>
    std::vector<SymGroupRepID> occ_symrep_IDs(BasicStructure<CoordType> const &_struc) {
      std::vector<SymGroupRepID> result;
      result.resize(_struc.basis().size());
      for(Index b = 0; b < _struc.basis().size(); ++b) {
        result[b] = _struc.basis()[b].occupant_dof().symrep_ID();
      }
      return result;
    }
    //************************************************************

    template<typename CoordType>
    std::map<DoFKey, DoFSetInfo> global_dof_info(BasicStructure<CoordType> const &_struc) {
      std::map<DoFKey, DoFSetInfo> result;
      for(auto const &dof :  _struc.global_dofs())
        result.emplace(dof.first, dof.second.info());

      return result;
    }

    //************************************************************

    template<typename CoordType>
    std::map<DoFKey, std::vector<DoFSetInfo> > local_dof_info(BasicStructure<CoordType> const &_struc) {
      std::map<DoFKey, std::vector<DoFSetInfo> > result;

      for(DoFKey const &type : continuous_local_dof_types(_struc)) {
        std::vector<DoFSetInfo> tresult(_struc.basis().size(), DoFSetInfo(SymGroupRepID(), Eigen::MatrixXd::Zero(DoF::BasicTraits(type).dim(), 0)));

        for(Index b = 0; b < _struc.basis().size(); ++b) {
          if(_struc.basis()[b].has_dof(type)) {
            tresult[b] = _struc.basis()[b].dof(type).info();
          }
        }
        result.emplace(type, std::move(tresult));
      }
      return result;
    }

    //************************************************************

    template<typename CoordType>
    std::map<DoFKey, Index> local_dof_dims(BasicStructure<CoordType> const &_struc) {
      std::map<DoFKey, Index> result;
      for(DoFKey const &type : continuous_local_dof_types(_struc))
        result[type] = local_dof_dim(type, _struc);

      return result;
    }


    //************************************************************
    template<typename CoordType>
    std::map<DoFKey, Index> global_dof_dims(BasicStructure<CoordType> const &_struc) {
      std::map<DoFKey, Index> result;
      for(auto const &type : _struc.global_dofs())
        result[type.first] = type.second.size();
      return result;
    }

    //************************************************************

    template<typename CoordType>
    Index local_dof_dim(DoFKey const &_name, BasicStructure<CoordType> const &_struc) {
      Index result = 0;
      for(CoordType const &site : _struc.basis()) {
        if(site.has_dof(_name))
          result = max(result, site.dof(_name).size());
      }
      return result;
    }

  }
}

#endif
