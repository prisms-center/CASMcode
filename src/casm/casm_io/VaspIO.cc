#include "casm/casm_io/VaspIO.hh"

#include "casm/crystallography/BasicStructure.hh"
#include "casm/crystallography/SimpleStructure.hh"
#include "casm/crystallography/Site.hh"
#include "casm/clex/Configuration.hh"
#include "casm/basis_set/DoF.hh"

namespace CASM {
  namespace VaspIO {

    // --- Definitions -------------------------------------------------------- //



    /// \brief Print POSCAR, providing a range of std::tuple<AtomName, Coordinate, SelectiveDynamics>
    template<typename TupleIterator>
    void PrintPOSCAR::_print(std::ostream &sout,
                             TupleIterator begin,
                             TupleIterator end) const {
      Log log(sout);
      this->_print(log, begin, end);
    }

    /// \brief Print POSCAR, providing a range of std::tuple<AtomName, Coordinate, SelectiveDynamics>
    template<typename TupleIterator>
    void PrintPOSCAR::_print(Log &sout,
                             TupleIterator begin,
                             TupleIterator end) const {

      int tprec = sout.ostream().precision();
      std::ios::fmtflags tflags = sout.ostream().flags();
      int prec = 8;
      sout.ostream().precision(prec);

      typedef std::tuple<std::string, Coordinate, SelectiveDynamics> tuple_type;

      // first filter out all atoms we are going to ignore, the remaining atoms get put in 'atom'
      std::vector<tuple_type> atom;
      for(auto it = begin; it != end; ++it) {
        // if Atom's name is not found in the ignore list, add it to 'atom'
        if(m_ignore.count(std::get<0>(*it)) == 0)
          atom.push_back(*it);

      }

      // print title, scale, and lattice
      sout << sout.indent_str() << m_title << "\n";
      sout << sout.indent_str() << std::fixed << std::setprecision(prec) << m_scale << "\n";

      sout.ostream().flags(std::ios::showpoint | std::ios::fixed | std::ios::right);

      int width = 12;
      width = print_matrix_width(sout, m_lat.lat_column_mat().transpose(), width);
      for(auto it = atom.cbegin(); it != atom.cend(); ++it) {
        Eigen::Vector3d vec;
        if(m_coord_mode == CART) vec = std::get<1>(*it).cart();
        else if(m_coord_mode == FRAC) vec = std::get<1>(*it).frac();
        width = print_matrix_width(sout, vec.transpose(), width);
      }
      Eigen::IOFormat format(prec, width + 1);

      sout << sout.indent_str() << m_lat[0].transpose().format(format) << '\n';
      sout << sout.indent_str() << m_lat[1].transpose().format(format) << '\n';
      sout << sout.indent_str() << m_lat[2].transpose().format(format) << '\n';

      // if after filtering out ignored atoms none are left, return
      if(atom.size() == 0) {
        return;
      }

      // count number each atom, and optionally print atom names line
      std::vector<int> atom_count = {1};
      auto it = atom.cbegin();
      std::string curr_atom = std::get<0>(*it);
      if(m_atom_names) {
        sout << sout.indent_str() << curr_atom << " ";
      }
      ++it;
      for(; it != atom.cend(); ++it) {
        if(std::get<0>(*it) != curr_atom) {
          atom_count.push_back(1);
          curr_atom = std::get<0>(*it);
          if(m_atom_names) {
            sout << curr_atom << " ";
          }
        }
        else {
          atom_count.back()++;
        }
      }
      if(m_atom_names) {
        sout << "\n";
      }

      sout << sout.indent_str();
      for(int i = 0; i < atom_count.size(); i++) {
        sout << atom_count[i] << " ";
      }
      sout << "\n";

      // print 'Selective Dynamics' if using selective dynamics
      if(m_sel_dynamics) {
        sout << sout.indent_str() << "Selective Dynamics\n";
      }

      // print coord mode
      sout << sout.indent_str() << COORD_MODE::NAME(m_coord_mode) << "\n";

      // print all coordinates, and seletive dynamics settings, and atom names if applicable
      for(auto it = atom.cbegin(); it != atom.cend(); ++it) {
        sout << sout.indent_str();
        std::get<1>(*it).print(sout, m_coord_mode, 0, format);
        if(m_sel_dynamics) {
          sout << " " << std::get<2>(*it);
        }
        if(m_append_atom_names) {
          sout << " " << std::get<0>(*it);
        }
        sout << "\n";
      }
      sout << "\n";

      sout.ostream().precision(tprec);
      sout.ostream().flags(tflags);

    }

    template void PrintPOSCAR::_print(
      std::ostream &sout,
      PrintPOSCAR::const_iterator begin,
      PrintPOSCAR::const_iterator end) const;

    template void PrintPOSCAR::_print(
      Log &sout,
      PrintPOSCAR::const_iterator begin,
      PrintPOSCAR::const_iterator end) const;


    /// \brief Construct PrintPOSCAR object
    ///
    /// By default:
    /// - title = struc.title
    /// - scale = 1.0
    /// - coordinate mode = frac (Direct)
    /// - atom names line is printed
    /// - Atoms printed in order appearing in the Structure. (No sorting by atom types)
    /// - No selective dynamics
    /// - {"Va", "va", "VA"} atoms not printed
    /// - not sorted
    ///
    /// Currently:
    /// - no displacement
    /// - assumes all species are atomic
    ///
    PrintPOSCAR::PrintPOSCAR(const BasicStructure<Site> &struc) :
      PrintPOSCAR() {
      m_lat = struc.lattice();

      set_title(struc.title());

      // create tuples collecting (Atom name, Coordinate, SelectiveDynamics) for each site
      for(int i = 0; i < struc.basis().size(); ++i) {
        m_atom_order.push_back(
          tuple_type(
            struc.basis()[i].occ_name(),
            struc.basis()[i],
            SelectiveDynamics()
          )
        );
      }
    }


    /// \brief Construct PrintPOSCAR object
    ///
    /// By default:
    /// - title = config.name()
    /// - scale = 1.0
    /// - coordinate mode = frac (Direct)
    /// - atom names line is printed
    /// - Atoms printed in order appearing in the Configuration. (No sorting by atom types)
    /// - No selective dynamics
    /// - {"Va", "va", "VA"} atoms not printed
    /// - also adds displacements and deformation if they are present as DoF
    ///
    /// Currently:
    /// - assumes molecules should be printed as their individual atoms.
    ///
    PrintPOSCAR::PrintPOSCAR(const SimpleStructure &sstruc, std::string const &_title) :
      PrintPOSCAR() {
      m_lat = Lattice(sstruc.lat_column_mat);
      set_title(_title);

      // create tuples collecting (Atom name, Coordinate, SelectiveDynamics) for each site
      for(Index i = 0; i < sstruc.n_atom(); ++i) {
        m_atom_order.push_back(
          tuple_type(
            sstruc.atom_info.names[i],
            Coordinate(sstruc.atom_info.coords.col(i), lattice(), CART),
            sstruc.atom_info.SD.col(i)
          )
        );
      }

    }

    PrintPOSCAR::PrintPOSCAR(const Configuration &config) :
      PrintPOSCAR(SimpleStructure(config), config.name()) {

    }
    /// \brief Construct PrintPOSCAR object
    ///
    /// By default:
    /// - title = ""
    /// - scale = 1.0
    /// - coordinate mode = frac (Direct)
    /// - atom names line is printed
    /// - Atoms printed in order appearing in the ConfigDoF. (No sorting by atom types)
    /// - No selective dynamics
    /// - {"Va", "va", "VA"} atoms not printed
    ///
    /// Currently:
    /// - no displacement
    /// - assumes all species are atomic
    ///
    /*
    PrintPOSCAR::PrintPOSCAR(const Supercell &scel, const ConfigDoF &configdof) :
      PrintPOSCAR(scel.lattice()) {

      // get occupant name for site i in configdof
      auto occ_name = [&](int i) {
                        return scel.prim().basis()[scel.sublat(i)].site_occupant()[configdof.occ(i)].name();
                      };

      // create tuples collecting (Atom name, Coordinate, SelectiveDynamics) for each site
      for(int i = 0; i < configdof.size(); ++i) {
        m_atom_order.push_back(
                               tuple_type(
                                          occ_name(i),
                                          scel.coord(i), // no displacement yet
                                          SelectiveDynamics()
                                          )
                               );
      }

    }
    */

    /// \brief Default sort is by species name
    void PrintPOSCAR::sort() {
      std::sort(m_atom_order.begin(),
                m_atom_order.end(),
      [ = ](const tuple_type & A, const tuple_type & B) {
        return std::get<0>(A) < std::get<0>(B);
      });
    }

    /// \brief Print POSCAR to stream
    void PrintPOSCAR::print(std::ostream &sout) const {
      _print(sout, m_atom_order.begin(), m_atom_order.end());
    }

    /// \brief Print POSCAR to log (enables indentation)
    void PrintPOSCAR::print(Log &sout) const {
      _print(sout, m_atom_order.begin(), m_atom_order.end());
    }

  }
}
