#include "casm/casm_io/VaspIO.hh"

#include "casm/crystallography/Structure.hh"
#include "casm/clex/Supercell.hh"
#include "casm/clex/Configuration.hh"

namespace CASM {
  namespace VaspIO {

    // --- Definitions -------------------------------------------------------- //

    namespace vaspio_impl {

      /// \brief Print POSCAR, providing a range of std::tuple<AtomName, Coordinate, SelectiveDynamics>
      template<typename TupleIterator>
      void PrintPOSCARBase::_print(std::ostream &sout,
                                   TupleIterator begin,
                                   TupleIterator end) const {

        int tprec = sout.precision();
        std::ios::fmtflags tflags = sout.flags();
        sout.precision(8);

        typedef std::tuple<std::string, Coordinate, SelectiveDynamics> tuple_type;

        // first filter out all atoms we are going to ignore, the remaining atoms get put in 'atom'
        std::vector<tuple_type> atom;
        for(auto it = begin; it != end; ++it) {

          // if Atom's name is not found in the ignore list, add it to 'atom'
          if(m_ignore.cend() == std::find_if(m_ignore.cbegin(),
                                             m_ignore.cend(),
          [&](const std::string & name) {
          if(std::get<0>(*it) == name) {
              return true;
            }
            return false;
          })) {
            atom.push_back(*it);
          }
        }

        // print title, scale, and lattice
        sout << m_title << "\n";
        sout << std::fixed << std::setprecision(8) << m_scale << "\n";

        sout.flags(std::ios::showpoint | std::ios::fixed | std::ios::right);

        sout << ' ' << std::setw(16) << m_lat[0].transpose() << '\n';
        sout << ' ' << std::setw(16) << m_lat[1].transpose() << '\n';
        sout << ' ' << std::setw(16) << m_lat[2].transpose() << '\n';

        // if after filtering out ignored atoms none are left, return
        if(atom.size() == 0) {
          return;
        }

        // count number each atom, and optionally print atom names line
        std::vector<int> atom_count = {1};
        auto it = atom.cbegin();
        std::string curr_atom = std::get<0>(*it);
        if(m_atom_names) {
          sout << curr_atom << " ";
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

        for(int i = 0; i < atom_count.size(); i++) {
          sout << atom_count[i] << " ";
        }
        sout << "\n";

        // print 'Selective Dynamics' if using selective dynamics
        if(m_sel_dynamics) {
          sout << "Selective Dynamics\n";
        }

        // print coord mode
        sout << COORD_MODE::NAME(m_coord_mode) << "\n";

        // print all coordinates, and seletive dynamics settings, and atom names if applicable
        for(auto it = atom.cbegin(); it != atom.cend(); ++it) {
          std::get<1>(*it).print(sout, m_coord_mode);
          if(m_sel_dynamics) {
            sout << " " << std::get<2>(*it);
          }
          if(m_append_atom_names) {
            sout << " " << std::get<0>(*it);
          }
          sout << "\n";
        }
        sout << "\n";

        sout.precision(tprec);
        sout.flags(tflags);

      }

      template void PrintPOSCARBase::_print(
        std::ostream &sout,
        PrintPOSCAR::const_iterator begin,
        PrintPOSCAR::const_iterator end) const;
    }

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
      vaspio_impl::PrintPOSCARBase(struc.lattice()) {

      set_title(struc.title);

      // create tuples collecting (Atom name, Coordinate, SelectiveDynamics) for each site
      for(int i = 0; i < struc.basis.size(); ++i) {
        m_atom_order.push_back(
          tuple_type(
            struc.basis[i].occ_name(),
            struc.basis[i],
            SelectiveDynamics()
          )
        );
      }
    }

    namespace {

      /// \brief returns a lattice with deformation if present
      /// deformed_lattice = F * reference_lattice
      /// F is the deformation matrix
      ///
      Lattice get_lattice(const Configuration &config) {
        const Lattice &ref_lat = config.supercell().lattice();
        if(config.has_deformation()) {
          return Lattice(config.deformation() * ref_lat.lat_column_mat());
        }
        else {
          return ref_lat;
        }
      }

      /// \brief returns a coordinate including displacements and  deformation if present
      /// deformed_coordinate = F * (reference_coordinate + displacement)
      /// F is the deformation matrix
      ///
      Coordinate get_coord(const Configuration &config, Index i, const Lattice &deformed_lat) {
        Coordinate ref_coord = config.supercell().coord(i);
        if(config.has_displacement()) {
          ref_coord.cart() += config.disp(i);
        }
        if(config.has_deformation()) {
          ref_coord.cart() = config.deformation() * ref_coord.const_cart();
        }
        Coordinate deformed_coord(ref_coord.const_cart(), deformed_lat, CART);
        deformed_coord.within();
        return deformed_coord;
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
    /// - assumes all species are atomic
    ///
    PrintPOSCAR::PrintPOSCAR(const Configuration &config) :
      vaspio_impl::PrintPOSCARBase(get_lattice(config)) {

      set_title(config.name());

      const Supercell &scel = config.supercell();

      // create tuples collecting (Atom name, Coordinate, SelectiveDynamics) for each site
      for(int i = 0; i < config.size(); ++i) {
        // Coordinate coord(get_coord(config, i).const_cart(), lattice(), CART);
        // coord.within();
        m_atom_order.push_back(
          tuple_type(
            config.mol(i).name,
            get_coord(config, i, lattice()),
            SelectiveDynamics()
          )
        );
      }

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
    PrintPOSCAR::PrintPOSCAR(const Supercell &scel, const ConfigDoF &configdof) :
      vaspio_impl::PrintPOSCARBase(scel.lattice()) {

      // get occupant name for site i in configdof
      auto occ_name = [&](int i) {
        return scel.prim().basis[scel.sublat(i)].site_occupant()[configdof.occ(i)].name;
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

  }
}
