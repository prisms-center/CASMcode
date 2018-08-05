#ifndef CASM_VaspIO
#define CASM_VaspIO

#include "casm/CASM_global_definitions.hh"
#include "casm/crystallography/BasicStructure.hh"
#include "casm/crystallography/Site.hh"
#include "casm/clex/Supercell.hh"
#include "casm/clex/Configuration.hh"

namespace CASM {

  /// \brief Functions and classes related to VASP input/output
  namespace VaspIO {

    /** \addtogroup VaspIO

      \brief Functions and classes related to VASP input/output

      \ingroup casmIO

    */

    /// \brief Store SelectiveDynamics options
    ///
    /// \ingroup VaspIO
    ///
    class SelectiveDynamics {

    public:

      /// \brief Default Constructor sets to false for all directions
      SelectiveDynamics() {
        m_option[0] = false;
        m_option[1] = false;
        m_option[2] = false;
      }

      SelectiveDynamics(bool x, bool y, bool z) {
        m_option[0] = x;
        m_option[1] = y;
        m_option[2] = z;
      }

      bool &operator[](int i) {
        return m_option[i];
      }

      const bool &operator[](int i) const {
        return m_option[i];
      }


    private:

      bool m_option[3];

    };

    /// \brief Write SelectiveDynamics options
    ///
    /// \relatesalso SelectiveDynamics
    ///
    inline std::ostream &operator<<(std::ostream &sout, const SelectiveDynamics &sel) {

      auto print = [&](bool val) {
        if(val) {
          sout << "T";
        }
        else {
          sout << "F";
        }
      };

      print(sel[0]);
      sout << " ";
      print(sel[1]);
      sout << " ";
      print(sel[2]);
      return sout;
    };


    namespace vaspio_impl {

      class PrintPOSCARBase {

      public:

        /// \brief Construct PrintPOSCAR object
        ///
        /// By default:
        /// - title = ""
        /// - scale = 1.0
        /// - coordinate mode = frac (Direct)
        /// - atom names line is printed
        /// - No selective dynamics
        /// - atom names appended to each coordinate line
        /// - {"Va", "va", "VA"} atoms ignored
        ///
        PrintPOSCARBase() :
          m_title(""),
          m_scale(1.0),
          m_coord_mode(FRAC),
          m_atom_names(true),
          m_sel_dynamics(false),
          m_append_atom_names(true),
          m_ignore {"VA", "Va", "va"} {}

        /// \brief Construct PrintPOSCAR object
        ///
        /// By default:
        /// - title = ""
        /// - scale = 1.0
        /// - coordinate mode = frac (Direct)
        /// - atom names line is printed
        /// - No selective dynamics
        /// - atom names appended to each coordinate line
        PrintPOSCARBase(const Lattice &lat) :
          PrintPOSCARBase() {
          m_lat = lat;
        }

        /// \brief Set title
        void set_title(std::string title) {
          m_title = title;
        }

        /// \brief Set scaling factor
        void set_scale(double s) {
          m_scale = s;
        }

        /// \brief Set coordinate mode to Direct (fractional)
        void set_direct() {
          m_coord_mode = FRAC;
        }
        /// \brief Set coordinate mode to fractional (Direct)
        void set_frac() {
          m_coord_mode = FRAC;
        }

        /// \brief Set coordinate mode to Cartesian
        void set_cart() {
          m_coord_mode = CART;
        }

        /// \brief Set coordinate mode
        void set_coord_mode(COORD_TYPE mode) {
          m_coord_mode = mode;
        }

        /// \brief Set selective dynamics off
        void set_selective_dynamics_off() {
          m_sel_dynamics = false;
        }

        /// \brief Set selective dynamics on
        void set_selective_dynamics_on() {
          m_sel_dynamics = true;
        }

        /// \brief Do not print atom names line
        void set_atom_names_off() {
          m_atom_names = false;
        }

        /// \brief Print atom names line
        void set_atom_names_on() {
          m_atom_names = true;
        }

        /// \brief Do not append atom name to end of each coordinate line
        void set_append_atom_names_off() {
          m_append_atom_names = false;
        }

        /// \brief Append atom name to end of each coordinate line
        void set_append_atom_names_on() {
          m_append_atom_names = false;
        }

        /// \brief Access vector of atom names which should not be printed, such as for vacancies
        std::vector<std::string> &ignore() {
          return m_ignore;
        }

        /// \brief const Access vector of atom names which should not be printed, such as for vacancies
        const std::vector<std::string> &ignore() const {
          return m_ignore;
        }


      protected:

        /// \brief Print POSCAR, provide a range of std::tuple<Atom name, Coordinate, SelectiveDynamics>
        template<typename TupleIterator>
        void _print(std::ostream &sout, TupleIterator begin, TupleIterator end);

      private:

        std::string m_title;
        double m_scale;
        COORD_TYPE m_coord_mode;
        bool m_atom_names;
        bool m_sel_dynamics;
        bool m_append_atom_names;
        Lattice m_lat;

        /// \brief List of atom names which should not be printed (primarily for vacancies)
        std::vector<std::string> m_ignore;

      };

    }

    /// \brief Print POSCAR with formating options
    ///
    /// Example:
    /// \code
    /// std::ostream file("POSCAR");
    /// Configuration config;
    /// PrintPOSCAR printer(config);
    /// printer.title("My system");
    /// printer.set_cart();
    /// printer.sort();
    /// printer.print(file);
    /// file.close();
    /// \endcode
    ///
    /// \ingroup VaspIO
    ///
    class PrintPOSCAR : public vaspio_impl::PrintPOSCARBase {

    public:

      typedef std::string AtomName;

      /// \brief Atom name, Coordinate, SelectiveDynamics
      typedef std::tuple<AtomName, Coordinate, SelectiveDynamics> tuple_type;

      typedef std::vector<tuple_type>::iterator iterator;

      typedef std::vector<tuple_type>::const_iterator const_iterator;


      /// \brief Construct PrintPOSCAR object
      PrintPOSCAR(const BasicStructure<Site> &struc);

      /// \brief Construct PrintPOSCAR object
      PrintPOSCAR(const Configuration &config);

      /// \brief Construct PrintPOSCAR object
      PrintPOSCAR(const Supercell &scel, const ConfigDoF &configdof);

      /// \brief Iterate over tuples of (AtomName, Coordinate, SelectiveDynamics)
      iterator begin() {
        return m_atom_order.begin();
      }

      /// \brief Iterate over tuples of (AtomName, Coordinate, SelectiveDynamics)
      iterator end() {
        return m_atom_order.end();
      }

      /// \brief Iterate over tuples of (AtomName, Coordinate, SelectiveDynamics)
      const_iterator cbegin() const {
        return m_atom_order.cbegin();
      }

      /// \brief Iterate over tuples of (AtomName, Coordinate, SelectiveDynamics)
      const_iterator cend() const {
        return m_atom_order.cend();
      }

      /// \brief Default sort is by atom name
      void sort();

      /// \brief Print POSCAR to stream
      void print(std::ostream &sout);


    private:

      /// \brief (AtomName, Coordinate, SelectiveDynamics)
      std::vector<tuple_type> m_atom_order;

    };


    // --- Definitions -------------------------------------------------------- //

    namespace vaspio_impl {

      /// \brief Print POSCAR, providing a range of std::tuple<AtomName, Coordinate, SelectiveDynamics>
      template<typename TupleIterator>
      void PrintPOSCARBase::_print(std::ostream &sout,
                                   TupleIterator begin,
                                   TupleIterator end) {

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
    inline PrintPOSCAR::PrintPOSCAR(const BasicStructure<Site> &struc) :
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
    ///
    /// Currently:
    /// - no displacement
    /// - assumes all species are atomic
    ///
    inline PrintPOSCAR::PrintPOSCAR(const Configuration &config) :
      vaspio_impl::PrintPOSCARBase(config.get_supercell().get_real_super_lattice()) {

      set_title(config.name());

      const Supercell &scel = config.get_supercell();

      // create tuples collecting (Atom name, Coordinate, SelectiveDynamics) for each site
      for(int i = 0; i < config.size(); ++i) {
        m_atom_order.push_back(
          tuple_type(
            config.get_mol(i).name,
            scel.coord(i), // no displacement yet
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
    inline PrintPOSCAR::PrintPOSCAR(const Supercell &scel, const ConfigDoF &configdof) :
      vaspio_impl::PrintPOSCARBase(scel.get_real_super_lattice()) {

      // get occupant name for site i in configdof
      auto occ_name = [&](int i) {
        return scel.get_prim().basis[scel.get_b(i)].site_occupant()[configdof.occ(i)].name;
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
    inline void PrintPOSCAR::sort() {
      std::sort(m_atom_order.begin(),
                m_atom_order.end(),
      [ = ](const tuple_type & A, const tuple_type & B) {
        return std::get<0>(A) < std::get<0>(B);
      });
    }

    /// \brief Print POSCAR to stream
    inline void PrintPOSCAR::print(std::ostream &sout) {
      _print(sout, m_atom_order.begin(), m_atom_order.end());
    }

  }


}

#endif
