#include <cstring>
#include "casm/app/casm_functions.hh"
#include "casm/app/AppIO.hh"
#include "casm/clusterography/ClusterOrbits.hh"
#include "casm/clex/PrimClex.hh"
#include "casm/clex/ClexBasis.hh"

namespace CASM {


  // ///////////////////////////////////////
  // 'clusters' function for casm
  //    (add an 'if-else' statement in casm.cpp to call this)

  int bset_command(const CommandArgs &args) {

    po::variables_map vm;

    /// Set command line options using boost program_options
    po::options_description desc("'casm bset' usage");
    desc.add_options()
    ("help,h", "Write help documentation")
    ("update,u", "Update basis set")
    ("orbits", "Pretty-print orbit prototypes")
    ("functions", "Pretty-print prototype cluster functions for each orbit")
    ("clusters", "Pretty-print all clusters")
    ("force,f", "Force overwrite");

    try {
      po::store(po::parse_command_line(args.argc, args.argv, desc), vm); // can throw
      bool call_help = false;

      /** --help option
       */
      if(vm.count("help") || call_help) {
        std::cout << "\n";
        std::cout << desc << std::endl;

        std::cout << "DESCRIPTION" << std::endl;
        std::cout << "    Generate and inspect cluster basis functions. A bspecs.json file should be available at\n"
                  << "        $ROOT/basis_set/$current_bset/bspecs.json\n"
                  << "    Run 'casm format --bspecs' for an example file.\n\n" ;

        return 0;
      }

      po::notify(vm); // throws on error, so do after help in case
      // there are any problems
    }
    catch(po::error &e) {
      std::cerr << desc << std::endl;
      std::cerr << "\nERROR: " << e.what() << std::endl << std::endl;
      return ERR_INVALID_ARG;
    }
    catch(std::exception &e) {
      std::cerr << desc << std::endl;
      std::cerr << "\nERROR: "  << e.what() << std::endl;
      return ERR_UNKNOWN;

    }

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

    // not sure how this will work yet...
    std::vector<std::string> dof_keys = {"occupation"};

    if(vm.count("update")) {

      // initialize project info
      DirectoryStructure dir(root);
      ProjectSettings set(root);
      Structure prim(read_prim(dir.prim()));

      std::cout << "\n***************************\n" << std::endl;


      if(!fs::is_regular_file(dir.bspecs(set.bset()))) {
        std::cout << "Error in 'casm bset': No basis set specifications file found at: " << dir.bspecs(set.bset()) << std::endl;
        return ERR_MISSING_INPUT_FILE;
      }


      std::vector<fs::path> filepaths({dir.clust(set.bset()),
                                       dir.basis(set.bset()),
                                       dir.clexulator_src(set.name(), set.bset()),
                                       dir.clexulator_o(set.name(), set.bset()),
                                       dir.clexulator_so(set.name(), set.bset())
                                      });

      bool any_existing_files = false;
      std::for_each(filepaths.cbegin(),
                    filepaths.cend(),
      [&](const fs::path & p) {
        if(fs::exists(p)) {
          if(!any_existing_files) {
            std::cout << "Existing files:\n";
            any_existing_files = true;
          }
          std::cout << "  " << p << "\n";
        }
      });

      std::cout << "\n";

      if(any_existing_files) {
        if(vm.count("force")) {
          std::cout << "Using --force. Will overwrite existing files.\n\n";
          fs::remove(dir.clexulator_src(set.name(), set.bset()));
          fs::remove(dir.clexulator_o(set.name(), set.bset()));
          fs::remove(dir.clexulator_so(set.name(), set.bset()));

          std::cout << "\n***************************\n" << std::endl;

        }
        else {
          std::cout << "Exiting due to existing files.  Use --force to force overwrite.\n\n";
          return ERR_EXISTING_FILE;
        }
      }

      jsonParser bspecs_json;
      std::vector<Orbit<IntegralCluster> > orbits;
      std::unique_ptr<ClexBasis> clex_basis;

      try {
        bspecs_json.read(dir.bspecs(set.bset()));

        args.log.construct("Orbitree");
        args.log << std::endl;

        make_orbits(prim,
                    prim.factor_group(),
                    bspecs_json,
                    set.crystallography_tol(),
                    alloy_sites_filter,
                    PrimPeriodicIntegralClusterSymCompare(set.crystallography_tol()),
                    std::back_inserter(orbits),
                    args.log);

        clex_basis.reset(new ClexBasis(prim));
        clex_basis->generate(orbits.begin(), orbits.end(), bspecs_json, dof_keys);

      }
      catch(std::exception &e) {
        args.err_log << e.what() << std::endl;
        return ERR_INVALID_INPUT_FILE;
      }

      // -- write clust.json ----------------
      {
        jsonParser clust_json;
        write_clust(orbits.begin(), orbits.end(), bspecs_json, clust_json);
        clust_json.write(dir.clust(set.bset()));

        args.log.write(dir.clust(set.bset()).string());
        args.log << std::endl;
      }


      // -- write basis.json ----------------
      {
        jsonParser basis_json;
        write_basis(orbits.begin(), orbits.end(), *clex_basis, basis_json, set.crystallography_tol());
        basis_json.write(dir.basis(set.bset()));

        args.log.write(dir.basis(set.bset()).string());
        args.log << std::endl;
      }


      // -- write global Clexulator
      {
        // get the neighbor list
        PrimNeighborList nlist(
          set.nlist_weight_matrix(),
          set.nlist_sublat_indices().begin(),
          set.nlist_sublat_indices().end()
        );

        // expand the nlist to contain sites in all orbits
        std::set<UnitCellCoord> nbors;
        prim_periodic_neighborhood(orbits.begin(), orbits.end(), std::inserter(nbors, nbors.begin()));
        nlist.expand(nbors.begin(), nbors.end());

        // write source code
        fs::ofstream outfile;
        outfile.open(dir.clexulator_src(set.name(), set.bset()));
        throw std::runtime_error("Error: print_clexulator is being re-implemented");
        //print_clexulator(*clex_basis, nlist, set.global_clexulator(), outfile, set.crystallography_tol());
        outfile.close();

        args.log.write(dir.clexulator_src(set.name(), set.bset()).string());
        args.log << std::endl;
      }

    }
    else if(vm.count("orbits") || vm.count("clusters") || vm.count("functions")) {

      DirectoryStructure dir(root);
      ProjectSettings set(root);

      if(!fs::exists(dir.clust(set.bset()))) {
        std::cerr << "ERROR: No 'clust.json' file found. Make sure to update your basis set with 'casm bset -u'.\n";
        return ERR_MISSING_DEPENDS;
      }

      PrimClex primclex(root, args.log);
      jsonParser clust_json(dir.clust(set.bset()));

      std::vector<Orbit<IntegralCluster> > orbits;
      read_clust(
        std::back_inserter(orbits),
        clust_json,
        primclex.prim(),
        primclex.prim().factor_group(),
        PrimPeriodicIntegralClusterSymCompare(set.crystallography_tol())
      );

      if(vm.count("orbits")) {
        print_clust(orbits.begin(), orbits.end(), args.log, ProtoSitesPrinter());
      }
      if(vm.count("clusters")) {
        print_clust(orbits.begin(), orbits.end(), args.log, FullSitesPrinter());
      }
      if(vm.count("functions")) {
        jsonParser bspecs_json;
        bspecs_json.read(dir.bspecs(set.bset()));

        ClexBasis clex_basis(primclex.prim());
        clex_basis.generate(orbits.begin(), orbits.end(), bspecs_json, dof_keys);

        print_clust(orbits.begin(), orbits.end(), args.log, ProtoFuncsPrinter(clex_basis));
      }
    }
    else {
      std::cerr << "\n" << desc << "\n";
    }

    return 0;
  };

}

