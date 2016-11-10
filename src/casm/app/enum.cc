#include <cstring>

#include "casm/app/casm_functions.hh"
#include "casm/clex/PrimClex.hh"
#include "casm/clex/ConfigEnumIterator.hh"
#include "casm/clex/FilteredConfigIterator.hh"
#include "casm/clex/ConfigIO.hh"
#include "casm/clex/ConfigEnumAllOccupations.hh"
#include "casm/clex/ConfigEnumIterator.hh"
#include "casm/completer/Handlers.hh"

namespace CASM {

  namespace Completer {
    EnumOption::EnumOption(): OptionHandlerBase("enum") {}

    int EnumOption::min_vol() const {
      return m_min_vol;
    }

    int EnumOption::max_vol() const {
      return m_max_vol;
    }

    const std::vector<std::string> &EnumOption::filter_strs() const {
      return m_filter_strs;
    }

    const fs::path &EnumOption::matrix_path() const {
      return m_matrix_path;
    }

    const std::string &EnumOption::lattice_directions() const {
      return m_lattice_directions_str;
    }

    void EnumOption::initialize() {
      add_help_suboption();
      add_scelnames_suboption();

      m_desc.add_options()
      ("min", po::value<int>(&m_min_vol)->default_value(1), "Min volume")
      ("max", po::value<int>(&m_max_vol), "Max volume")
      ("filter", po::value<std::vector<std::string> >(&m_filter_strs)->multitoken()->value_name(ArgHandler::query()), "Filter configuration enumeration so that only configurations matching a 'casm query'-type expression are recorded")
      ("all,a", "Enumerate configurations for all supercells")
      ("supercells,s", "Enumerate supercells")
      ("configs,c", "Enumerate configurations")
      ("matrix,m", po::value<fs::path>(&m_matrix_path)->value_name(ArgHandler::path()), "Specify a matrix to apply to the primitive cell before beginning enumeration")
      ("lattice-directions,z", po::value<std::string>(&m_lattice_directions_str)->default_value("abc"), "Restrict enumeration along a, b or c lattice vectors");

      return;
    }
  }


  // ///////////////////////////////////////
  // 'enum' function for casm
  //    (add an 'if-else' statement in casm.cpp to call this)

  int enum_command(const CommandArgs &args) {

    //casm enum [—supercell min max] [—config supercell ] [—hopconfigs hop.background]
    //- enumerate supercells and configs and hop local configurations

    int min_vol, max_vol;
    std::vector<std::string> scellname_list, filter_expr;
    //double tol;
    COORD_TYPE coordtype = CASM::CART;

    int dims = 3;
    Eigen::Matrix3i G = Eigen::Matrix3i::Identity();    //The matrix that defines the lattice vectors to enumerate over relative to the primitive vectors
    fs::path matrix_path;
    std::string ezmode = "abc"; //Assign default in ::initialize()?
    po::variables_map vm;

    Completer::EnumOption enum_opt;

    try {
      po::store(po::parse_command_line(args.argc, args.argv, enum_opt.desc()), vm); // can throw

      /** --help option
       */
      if(vm.count("help")) {
        std::cout << "\n";
        std::cout << enum_opt.desc() << std::endl;

        return 0;
      }

      if(vm.count("desc")) {
        std::cout << "\n";
        std::cout << enum_opt.desc() << std::endl;

        std::cout << "DESCRIPTION\n" << std::endl;

        std::cout << "  casm enum --supercells --max V                        \n"
                  "  - To enumerate supercells up to volume V (units: number\n"
                  "    of primitive cells)                                 \n"
                  "  - Use --matrix and --lattice-directions for restricted\n"
                  "    enumeration: \n\n";

        std::cout << "  --matrix" << std::endl;
        std::cout << "    - When using --supercells, you may use the this option             " << std::endl;
        std::cout << "      to specify a transformation matrix to apply to your primitive    " << std::endl;
        std::cout << "      lattice vectors before you begin your enumeration.               " << std::endl;
        std::cout << "      For example, if the PRIM used for your project was a primitive   " << std::endl;
        std::cout << "      FCC structure, you could point this option to a file M.txt with  " << std::endl;
        std::cout << "          -1  1  1\n           1 -1  1          \n           1  1 -1" << std::endl;
        std::cout << "      to enumerate only over conventional FCC supercells, since the    " << std::endl;
        std::cout << "      matrix shown above transforms a primitive FCC to a conventional  " << std::endl;
        std::cout << "      FCC." << std::endl;
        std::cout << "    - If this option isn't specified, the lattice vectors to           " << std::endl;
        std::cout << "      enumerate over will simply be the vectors of your PRIM structure." << std::endl;
        std::cout << "    - Note that the --max and --min values you specify will be         " << std::endl;
        std::cout << "      relative to the determinant of your matrix. For the provided     " << std::endl;
        std::cout << "      example, specifying --max 6 would result in supercells up to size" << std::endl;
        std::cout << "      24, since det(M)*6=24 (4 primitive lattice sites per conventional" << std::endl;
        std::cout << "      cell)." << std::endl;
        std::cout << std::endl;
        std::cout << "  --lattice-directions" << std::endl;
        std::cout << "    - When using --supercells, you may use the this option             " << std::endl;
        std::cout << "      restrict the supercell enumeration to 1, 2 or 3 of the lattice   " << std::endl;
        std::cout << "      vectors, to get 1,2 or 3-dimensional supercells. By directly     " << std::endl;
        std::cout << "      specifying combinations of 'a', 'b' and 'c', you determine which " << std::endl;
        std::cout << "      of the lattice vectors you want to enumerate over.               " << std::endl;
        std::cout << "      For example, to enumerate 1-dimensional supercells along the 'c' " << std::endl;
        std::cout << "      direction, simply specify '--lattice-directions c'. If you want       " << std::endl;
        std::cout << "      2-dimensional supercells along the a and c lattice vectors,      " << std::endl;
        std::cout << "      specify '--lattice-directions ac'.                                    " << std::endl;
        std::cout << "    - If this option isn't specified, 3-dimensional supercells will be " << std::endl;
        std::cout << "      enumerated, equivalent to specifying 'lattice-directions abc'.        " << std::endl;
        std::cout << "    - This option can be used in conjunction with the --matrix option. " << std::endl;
        std::cout << "      If this is the case, then the meaning of 'a', 'b' and 'c' changes" << std::endl;
        std::cout << "      from the lattice vectors of your PRIM, to the  vectors of the    " << std::endl;
        std::cout << "      lattice resulting from multiplying your PRIM by the specified    " << std::endl;
        std::cout << "      matrix.\n\n";

        std::cout << "  casm enum --configs --all'                             \n"
                  "  - To enumerate configurations for all supercells.      \n\n"

                  "  casm enum --configs --min MINV --max MAXV              \n"
                  "  - To enumerate configurations for supercells ranging in\n"
                  "    volume from MINV to MAXV (units: number of primitive \n"
                  "    cells).                                              \n\n"

                  "  casm enum --configs --scellname NAME                   \n"
                  "  - To enumerate configurations for a particular         \n"
                  "    supercell.                                           \n\n"

                  "  casm enum --configs [...] --filter '... casm query commands...' \n"
                  "  - To perform restricted enumeration of configurations  \n"
                  "    such that only configurations the return '1' (true)  \n"
                  "    for a particular 'casm query' command are retained.  \n\n";
        return 0;
      }

      po::notify(vm); // throws on error, so do after help in case
      // there are any problems

      min_vol = enum_opt.min_vol();
      max_vol = enum_opt.max_vol();
      scellname_list = enum_opt.supercell_strs();
      filter_expr = enum_opt.filter_strs();
      matrix_path = enum_opt.matrix_path();
      ezmode = enum_opt.lattice_directions();

      if(min_vol != 1 && !vm.count("max")) {
        std::cerr << "\n" << enum_opt.desc() << "\n" << std::endl;
        std::cerr << "Error in 'casm enum'. If --min is given, --max must also be given." << std::endl;
        return ERR_INVALID_ARG;
      }
      if(vm.count("lattice-directions")) {
        bool bad_args = false;
        if(ezmode.size() > 3 || ezmode.size() < 1) {
          bad_args = true;
        }
        for(int i = 0; i < ezmode.size(); i++) {
          if(ezmode[i] != 'a' && ezmode[i] != 'b' && ezmode[i] != 'c') {
            bad_args = true;
            break;
          }
        }
        if(bad_args) {
          std::cerr << std::endl << enum_opt.desc() << std::endl << std::endl;
          std::cerr << "When using --lattice-directions, specify the primitive lattice vectors you want to enumerate over with a string." << std::endl;
          std::cerr << "For example, to enumerate over only a and b, pass 'ab'. To enumerate over only c pass 'c'." << std::endl;

          return ERR_INVALID_ARG;
        }
      }
      if(vm.count("supercells") + vm.count("configs") != 1) {
        std::cerr << "\n" << enum_opt.desc() << "\n" << std::endl;
        std::cerr << "Error in 'casm enum'. Exactly one of either --supercells or --configs must be given." << std::endl;
        return ERR_INVALID_ARG;
      }
      if(vm.count("supercells") && !vm.count("max")) {
        std::cerr << "\n" << enum_opt.desc() << "\n" << std::endl;
        std::cerr << "Error in 'casm enum'. If --supercells is given, --max must be given." << std::endl;
        return ERR_INVALID_ARG;
      }
      if(vm.count("configs") && (vm.count("max") + vm.count("all") != 1)) {
        std::cerr << "\n" << enum_opt.desc() << "\n" << std::endl;
        std::cerr << "Error in 'casm enum'. If --configs is given, exactly one of either --max or --all must be given." << std::endl;
        return ERR_INVALID_ARG;
      }
    }
    catch(po::error &e) {
      std::cerr << enum_opt.desc() << std::endl;
      std::cerr << "ERROR: " << e.what() << std::endl << std::endl;
      return ERR_INVALID_ARG;
    }
    catch(std::exception &e) {
      std::cerr << enum_opt.desc() << std::endl;
      std::cerr << "ERROR: " << e.what() << std::endl << std::endl;
      return ERR_UNKNOWN;

    }

    COORD_MODE C(coordtype);

    const fs::path &root = args.root;
    if(root.empty()) {
      args.err_log.error("No casm project found");
      args.err_log << std::endl;
      return ERR_NO_PROJ;
    }

    // If 'args.primclex', use that, else construct PrimClex in 'uniq_primclex'
    // Then whichever exists, store reference in 'primclex'
    std::unique_ptr<PrimClex> uniq_primclex;
    PrimClex &primclex = make_primclex_if_not(args, uniq_primclex);
    const DirectoryStructure &dir = primclex.dir();
    const ProjectSettings &set = primclex.settings();

    if(vm.count("lattice-directions")) {
      dims = ezmode.size();
      while(ezmode.size() != 3) {
        if(std::find(ezmode.begin(), ezmode.end(), 'a') == ezmode.end()) {
          ezmode.push_back('a');
        }
        if(std::find(ezmode.begin(), ezmode.end(), 'b') == ezmode.end()) {
          ezmode.push_back('b');
        }
        if(std::find(ezmode.begin(), ezmode.end(), 'c') == ezmode.end()) {
          ezmode.push_back('c');
        }
      }
    }
    if(vm.count("matrix")) {
      fs::ifstream matstream(matrix_path);
      for(int i = 0; i < 9; i++) {
        matstream >> G;
      }
    }


    if(vm.count("supercells")) {
      std::cout << "\n***************************\n" << std::endl;

      std::cout << "Generating supercells from " << min_vol << " to " << max_vol << std::endl;
      std::cout << "Lattice vectors ";
      for(int i = 0; i < dims; i++) {
        if(i == dims - 1) {
          std::cout << "and ";
        }
        std::cout << ezmode[i] << " ";
      }
      std::cout << "will be used. (supercells will be " << dims << "-dimensional)" << std::endl;

      if(vm.count("matrix")) {
        std::cout << "The transformation matrix" << std::endl << G << std::endl << "will be applied to the primitive cell before enumeration begins" << std::endl;
      }
      std::cout << std::endl;

      primclex.generate_supercells(min_vol, max_vol, ezmode, G, true);
      std::cout << "\n  DONE." << std::endl << std::endl;

      std::cout << "Write SCEL." << std::endl << std::endl;
      primclex.print_supercells();

    }
    else if(vm.count("configs")) {
      if(vm.count("dimensions") || vm.count("matrix") || (dims < 3)) {
        std::cerr << "Option --configs in conjunction with limited supercell enumeration is currently unsupported" << std::endl;
        return ERR_INVALID_ARG;
      }

      std::cout << "\n***************************\n" << std::endl;
      std::vector<Supercell *> scel_selection;

      //Build the supercell selection based on user input
      if(vm.count("all")) {
        std::cout << "Enumerate all configurations" << std::endl << std::endl;
        for(int j = 0; j < primclex.get_supercell_list().size(); j++) {
          scel_selection.push_back(&(primclex.get_supercell(j)));
        }
      }
      else {
        if(vm.count("max")) {
          std::cout << "Enumerate configurations from volume " << min_vol << " to " << max_vol << std::endl << std::endl;
          for(int j = 0; j < primclex.get_supercell_list().size(); j++) {
            if(primclex.get_supercell(j).volume() >= min_vol && primclex.get_supercell(j).volume() <= max_vol) {
              scel_selection.push_back(&(primclex.get_supercell(j)));
            }
          }
        }
        if(vm.count("scellname")) {
          Index j;
          std::cout << "Enumerate configurations for named supercells" << std::endl << std::endl;
          for(int i = 0; i < scellname_list.size(); i++) {
            if(!primclex.contains_supercell(scellname_list[i], j)) {
              std::cout << "Error in 'casm enum'. Did not find supercell: " << scellname_list[i] << std::endl;
              return 1;
            }
            scel_selection.push_back(&(primclex.get_supercell(j)));
          }
        }
      }
      //\Finished building supercell selection

      if(scel_selection.empty()) {
        std::cout << "Did not find any supercells. Make sure to 'casm enum --supercells' first!" << std::endl << std::endl;
        return ERR_MISSING_DEPENDS;
      }

      //We have the selection. Now do enumeration
      for(auto it = scel_selection.begin(); it != scel_selection.end(); ++it) {
        std::cout << "  Enumerate configurations for " << (**it).get_name() << " ...  " << std::flush;

        ConfigEnumAllOccupations<Configuration> enumerator(**it);
        Index num_before = (**it).get_config_list().size();
        if(vm.count("filter")) {
          try {
            (**it).add_unique_canon_configs(filter_begin(enumerator.begin(), enumerator.end(), filter_expr, set.config_io()), filter_end(enumerator.end()));
          }
          catch(std::exception &e) {
            std::cerr << "Cannot filter configurations using the expression provided: \n" << e.what() << "\nExiting...\n";
            return ERR_INVALID_ARG;
          }
        }
        else {
          (**it).add_unique_canon_configs(enumerator.begin(), enumerator.end());
        }

        std::cout << ((**it).get_config_list().size() - num_before) << " configs." << std::endl;
      }
      std::cout << "  DONE." << std::endl << std::endl;

    }

    std::cout << "Writing config_list..." << std::endl;
    primclex.write_config_list();
    std::cout << "  DONE" << std::endl;

    std::cout << std::endl;

    return 0;
  };

}

