#include <cstring>
#include "casm/app/casm_functions.hh"

#include "casm/app/AppIO.hh"
#include "casm/clusterography/ClusterOrbits.hh"
#include "casm/clusterography/IntegralCluster.hh"
#include "casm/clex/PrimClex.hh"
#include "casm/clex/ClexBasis.hh"

//#include "casm/app/DirectoryStructure.hh"
//#include "casm/app/ProjectSettings.hh"
//#include "casm/crystallography/Structure.hh"
//#include "casm/app/AppIO.hh"
//#include "casm/clex/PrimClex.hh"
//#include "casm/clusterography/jsonClust.hh"


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
        args.log << "\n";
        args.log << desc << std::endl;

        args.log << "DESCRIPTION" << std::endl;
        args.log << "    Generate and inspect cluster basis functions. A bspecs.json file should be available at\n"
                 << "        $ROOT/basis_set/$current_bset/bspecs.json\n"
                 << "    Run 'casm format --bspecs' for an example file.\n\n" ;

        return 0;
      }

      po::notify(vm); // throws on error, so do after help in case
      // there are any problems
    }
    catch(po::error &e) {
      args.err_log << desc << std::endl;
      args.err_log << "\nERROR: " << e.what() << std::endl << std::endl;
      return ERR_INVALID_ARG;
    }
    catch(std::exception &e) {
      args.err_log << desc << std::endl;
      args.err_log << "\nERROR: "  << e.what() << std::endl;
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
    const DirectoryStructure &dir = primclex.dir();
    const ProjectSettings &set = primclex.settings();

    if(vm.count("update")) {

      // initialize project info
      Structure prim = primclex.prim();

      if(!fs::is_regular_file(dir.bspecs(set.bset()))) {
        args.err_log.error("'bspecs.json' file not found");
        args.err_log << "expected basis set specifications file at: " << dir.bspecs(set.bset()) << "\n" << std::endl;
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
            args.log.custom("Found existing files");
            any_existing_files = true;
          }
          args.log << "found: " << p << "\n";
        }
      });

      if(any_existing_files) {
        if(vm.count("force")) {
          args.log << "Using --force. Will overwrite existing files.\n" << std::endl;
          fs::remove(dir.clexulator_src(set.name(), set.bset()));
          fs::remove(dir.clexulator_o(set.name(), set.bset()));
          fs::remove(dir.clexulator_so(set.name(), set.bset()));
        }
        else {
          args.log << "Exiting due to existing files.  Use --force to force overwrite.\n" << std::endl;
          return ERR_EXISTING_FILE;
        }
      }

      jsonParser bspecs_json;
      std::vector<PrimPeriodicIntegralClusterOrbit> orbits;
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
        clex_basis->generate(orbits.begin(), orbits.end(), bspecs_json);

      }
      catch(std::exception &e) {
        args.err_log << e.what() << std::endl;
        return ERR_INVALID_INPUT_FILE;
      }

      // -- write clust.json ----------------
      {
        jsonParser clust_json;
        write_clust(orbits.begin(), orbits.end(), clust_json, ProtoSitesPrinter(), bspecs_json);
        clust_json.write(dir.clust(set.bset()));

        args.log.write(dir.clust(set.bset()).string());
        args.log << std::endl;
      }

      // -- write basis.json ----------------
      {
        jsonParser basis_json;
        write_clust(orbits.begin(), orbits.end(), basis_json, ProtoFuncsPrinter(*clex_basis), bspecs_json);
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

      if(!fs::exists(dir.clust(set.bset()))) {
        args.err_log.error("No 'clust.json' file found");
        args.err_log << "Make sure to update your basis set with 'casm bset -u'.\n" << std::endl;
        return ERR_MISSING_DEPENDS;
      }

      PrimClex primclex(root, args.log);
      jsonParser clust_json(dir.clust(set.bset()));

      std::vector<PrimPeriodicIntegralClusterOrbit> orbits;
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
        clex_basis.generate(orbits.begin(), orbits.end(), bspecs_json);

        print_clust(orbits.begin(), orbits.end(), args.log, ProtoFuncsPrinter(clex_basis));
      }
    }
    else {
      args.err_log.error("Unknown error");
      args.err_log << desc << "\n" << std::endl;
    }

    return 0;
  };

}

