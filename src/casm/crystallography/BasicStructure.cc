#include "casm/crystallography/BasicStructure.hh"
#include "casm/basis_set/DoFIsEquivalent.hh"
#include "casm/basis_set/OccupationDoFTraits.hh"
#include "casm/crystallography/LatticeIsEquivalent.hh"
#include "casm/crystallography/LatticePointWithin.hh"
#include "casm/crystallography/Niggli.hh"
#include "casm/crystallography/Site.hh"
#include "casm/crystallography/Molecule.hh"
#include "casm/crystallography/UnitCellCoord.hh"
#include "casm/basis_set/DoF.hh"
#include <boost/filesystem.hpp>

namespace CASM {
  namespace xtal {
    BasicStructure::BasicStructure(const fs::path &filepath) : m_lattice() {
      if(!fs::exists(filepath)) {
        std::cerr << "Error in BasicStructure::BasicStructure(const fs::path &filepath)." << std::endl;
        std::cerr << "  File does not exist at: " << filepath << std::endl;
        exit(1);
      }
      fs::ifstream infile(filepath);
      read(infile);
    }

    //***********************************************************

    BasicStructure::BasicStructure(const BasicStructure &RHS) :
      m_lattice(RHS.lattice()),
      m_title(RHS.title()),
      m_basis(RHS.basis()),
      m_global_dof_map(RHS.m_global_dof_map) {
      for(Index i = 0; i < basis().size(); i++) {
        m_basis[i].set_lattice(lattice(), CART);
      }
    }

    //***********************************************************

    BasicStructure &BasicStructure::operator=(const BasicStructure &RHS) {
      m_lattice = RHS.lattice();
      m_title = RHS.title();
      set_basis(RHS.basis());
      m_global_dof_map = RHS.m_global_dof_map;

      for(Index i = 0; i < basis().size(); i++)
        m_basis[i].set_lattice(lattice(), CART);

      return *this;
    }

    //************************************************************

    DoFSet const &BasicStructure::global_dof(std::string const &_dof_type) const {
      auto it = m_global_dof_map.find(_dof_type);
      if(it != m_global_dof_map.end())
        return (it->second);
      else
        throw std::runtime_error(std::string("In BasicStructure::dof(), this structure does not contain any global DoF's of type " + _dof_type));

    }

    void BasicStructure::reset() {
      set_site_internals();

      within();
    }

    //***********************************************************

    void BasicStructure::update() {
      set_site_internals();
    }

    //***********************************************************

    void BasicStructure::within() {
      for(Index i = 0; i < basis().size(); i++) {
        m_basis[i].within();
      }
      return;
    }

    /**
     * It is NOT wise to use this function unless you have already
     * initialized a superstructure with lattice vectors.
     *
     * It is more wise to use the two methods that call this method:
     * Either the overloaded * operator which does:
     *  SCEL_Lattice * Prim_Structrue = New_Superstructure
     *       --- or ---
     *  New_Superstructure=Prim_BasicStructure.create_superstruc(SCEL_Lattice);
     *
     *  Both of these will return NEW superstructures.
     */

    void BasicStructure::fill_supercell(const BasicStructure &prim) {
      Index i, j;

      auto all_lattice_points = make_lattice_points(prim.lattice(), lattice(), lattice().tol());

      m_basis.clear();

      //loop over basis sites of prim
      for(j = 0; j < prim.basis().size(); j++) {

        //loop over prim_grid points
        for(const auto &lattice_point : all_lattice_points) {
          Coordinate lattice_point_coordinate = make_superlattice_coordinate(lattice_point, prim.lattice(), lattice());

          //push back translated basis site of prim onto superstructure basis
          push_back(prim.basis()[j] + lattice_point_coordinate);

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

    BasicStructure BasicStructure::create_superstruc(const Lattice &scel_lat) const {
      BasicStructure tsuper(scel_lat);
      tsuper.fill_supercell(*this);
      return tsuper;
    }

    //***********************************************************

    void BasicStructure::set_site_internals() {
      for(Index nb = 0; nb < basis().size(); nb++) {
        m_basis[nb].set_basis_ind(nb);
      }
    }

    //***********************************************************

    void BasicStructure::set_lattice(const Lattice &new_lat, COORD_TYPE mode) {

      m_lattice = new_lat;

      for(Index nb = 0; nb < basis().size(); nb++) {
        m_basis[nb].set_lattice(lattice(), mode);
      }
    }

    //***********************************************************


    void BasicStructure::set_title(std::string const &_title) {
      m_title = _title;
    }

    //\Liz D 032514
    //***********************************************************
    /**
     * Allows for the basis elements of a basic structure to be
     * manually set, e.g. as in jsonParser.cc.
     */
    //***********************************************************


    void BasicStructure::set_basis(std::vector<Site> const &_basis, COORD_TYPE mode) {
      m_basis.clear();
      m_basis.reserve(_basis.size());
      for(Site const &site : _basis)
        push_back(site, mode);

    }

    void BasicStructure::clear_basis() {
      m_basis.clear();
      this->reset();

    }

    void BasicStructure::push_back(Site const &_site, COORD_TYPE mode) {
      m_basis.push_back(_site);
      m_basis.back().set_basis_ind(basis().size() - 1);
      m_basis.back().set_lattice(lattice(), mode);
    }

    //************************************************************
    /// Counts sites that allow vacancies
    Index BasicStructure::max_possible_vacancies()const {
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

    void BasicStructure::read(std::istream &stream) {
      int i, t_int;
      char ch;
      std::vector<double> num_elem;
      std::vector<std::string> elem_array;
      bool read_elem = false;
      std::string tstr;
      std::stringstream tstrstream;

      Site tsite(lattice());

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

    void BasicStructure::print_xyz(std::ostream &stream, bool frac) const {
      stream << basis().size() << '\n';
      stream << title() << '\n';
      stream.precision(7);
      stream.width(11);
      stream.flags(std::ios::showpoint | std::ios::fixed | std::ios::right);
      stream << "      a       b       c" << '\n';
      stream << lattice().lat_column_mat() << '\n';
      for(Index i = 0; i < basis().size(); i++) {
        std::string site_label = basis()[i].allowed_occupants().size() == 1 ?  basis()[i].allowed_occupants()[0] : "?" ;
        stream << std::setw(2) << site_label << " ";
        if(frac) {
          stream << std::setw(12) << basis()[i].frac().transpose() << '\n';
        }
        else {
          stream << std::setw(12) << basis()[i].cart() << '\n';
        }
      }

    }

    //***********************************************************

    BasicStructure &BasicStructure::operator+=(const Coordinate &shift) {

      for(Index i = 0; i < basis().size(); i++) {
        m_basis[i] += shift;
      }

      return (*this);
    }


    //***********************************************************

    BasicStructure &BasicStructure::operator-=(const Coordinate &shift) {

      for(Index i = 0; i < basis().size(); i++) {
        m_basis[i] -= shift;
      }
      //factor_group -= shift;
      //asym_unit -= shift;
      return (*this);
    }

    //***********************************************************

    BasicStructure operator*(const Lattice &LHS, const BasicStructure &RHS) {
      BasicStructure tsuper(LHS);
      tsuper.fill_supercell(RHS);
      return tsuper;
    }

    //***********************************************************

    /// \brief Returns true if structure has attributes affected by time reversal
    // private for now, expose if necessary
    bool BasicStructure::is_time_reversal_active() const {
      for(auto const &dof : m_global_dof_map)
        if(dof.second.traits().time_reversal_active())
          return true;
      for(Site const &site : basis())
        if(site.time_reversal_active())
          return true;
      return false;
    }

    //***********************************************************

    std::vector<UnitCellCoord> symop_site_map(SymOp const &_op, BasicStructure const &_struc) {
      return symop_site_map(_op, _struc, _struc.lattice().tol());
    }

    //***********************************************************

    std::vector<UnitCellCoord> symop_site_map(SymOp const &_op, BasicStructure const &_struc, double _tol) {
      std::vector<UnitCellCoord> result;
      // Determine how basis sites transform from the origin unit cell
      for(int b = 0; b < _struc.basis().size(); b++) {
        Site transformed_basis_site = _op * _struc.basis()[b];
        result.emplace_back(UnitCellCoord::from_coordinate(_struc, transformed_basis_site, _tol));
      }
      return result;
    }

    //************************************************************

    /// Returns an std::vector of each *possible* Specie in this Structure
    std::vector<std::string> struc_species(BasicStructure const &_struc) {

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
    std::vector<Molecule> struc_molecule(BasicStructure const &_struc) {

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
    std::vector<std::string> struc_molecule_name(BasicStructure const &_struc) {

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
    std::vector<std::vector<std::string> > allowed_molecule_unique_names(BasicStructure const &_struc) {
      using IPair = std::pair<Index, Index>;
      std::map<std::string, std::vector<Molecule> > name_map;
      std::map<std::string, IPair> imap;

      std::vector<std::vector<std::string> > result(_struc.basis().size());
      for(Index b = 0; b < _struc.basis().size(); ++b) {
        for(Index j = 0; j < _struc.basis()[b].occupant_dof().size(); ++j) {
          Molecule const &mol(_struc.basis()[b].occupant_dof()[j]);
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
    std::vector<std::vector<std::string> > allowed_molecule_names(BasicStructure const &_struc) {
      std::vector<std::vector<std::string> > result(_struc.basis().size());

      for(Index b = 0; b < _struc.basis().size(); ++b)
        result[b] = _struc.basis()[b].allowed_occupants();

      return result;
    }

    //************************************************************

    std::vector<DoFKey> continuous_local_dof_types(BasicStructure const &_struc) {
      std::set<std::string> tresult;

      for(Site const &site : _struc.basis()) {
        auto sitetypes = site.dof_types();
        tresult.insert(sitetypes.begin(), sitetypes.end());
      }
      return std::vector<std::string>(tresult.begin(), tresult.end());
    }

    //************************************************************

    std::vector<DoFKey> all_local_dof_types(BasicStructure const &_struc) {
      std::set<std::string> tresult;

      for(Site const &site : _struc.basis()) {
        auto sitetypes = site.dof_types();
        tresult.insert(sitetypes.begin(), sitetypes.end());
        if(site.occupant_dof().size() > 1) {
          tresult.insert(DoFType::occupation().name());
        }
      }
      return std::vector<std::string>(tresult.begin(), tresult.end());
    }


    //************************************************************

    std::vector<DoFKey> global_dof_types(BasicStructure const &_struc) {
      std::vector<std::string> result;
      for(auto const &dof :  _struc.global_dofs())
        result.push_back(dof.first);
      return result;
    }


    //************************************************************

    std::vector<SymGroupRepID> occ_symrep_IDs(BasicStructure const &_struc) {
      std::vector<SymGroupRepID> result;
      result.resize(_struc.basis().size());
      for(Index b = 0; b < _struc.basis().size(); ++b) {
        result[b] = _struc.basis()[b].occupant_dof().symrep_ID();
      }
      return result;
    }
    //************************************************************

    std::map<DoFKey, DoFSetInfo> global_dof_info(BasicStructure const &_struc) {
      std::map<DoFKey, DoFSetInfo> result;
      for(auto const &dof :  _struc.global_dofs())
        result.emplace(dof.first, dof.second.info());

      return result;
    }

    //************************************************************

    std::map<DoFKey, std::vector<DoFSetInfo> > local_dof_info(BasicStructure const &_struc) {
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

    std::map<DoFKey, Index> local_dof_dims(BasicStructure const &_struc) {
      std::map<DoFKey, Index> result;
      for(DoFKey const &type : continuous_local_dof_types(_struc))
        result[type] = local_dof_dim(type, _struc);

      return result;
    }


    //************************************************************
    std::map<DoFKey, Index> global_dof_dims(BasicStructure const &_struc) {
      std::map<DoFKey, Index> result;
      for(auto const &type : _struc.global_dofs())
        result[type.first] = type.second.size();
      return result;
    }

    //************************************************************

    Index local_dof_dim(DoFKey const &_name, BasicStructure const &_struc) {
      Index result = 0;
      for(Site const &site : _struc.basis()) {
        if(site.has_dof(_name))
          result = max(result, site.dof(_name).size());
      }
      return result;
    }
  }
}
