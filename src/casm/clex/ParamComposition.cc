#include "casm/clex/ParamComposition.hh"

#include <cmath>
#include <unistd.h>

#include <set>
#include <map>
#include "casm/misc/algorithm.hh"
#include "casm/misc/CASM_Eigen_math.hh"

namespace CASM {


  //*************************************************************
  //GENERATE Routines

  namespace Local {
    // Takes set of subsystem compomponents, given as [set of indices], and set of sublattices, given as the pairs {allowed components [set of indeces], multiplicity [index]},
    // and total number of components in the larger system. Returns the set of extreme integer compositions for the subsystem
    static std::vector<Eigen::VectorXi> _sub_extremes(std::set<Index> const &subcomponents,  std::map<std::set<Index>, Index>  const &sublats, Index dim) {
      std::vector<Eigen::VectorXi> result;

      std::vector<Index> compon(subcomponents.begin(), subcomponents.end());

      // Count over k-combinations of subsystem components, where k is number of sublattices
      Index k = sublats.size();
      Index n = compon.size();
      Index ncomb = nchoosek(n, k);

      // Combination is stored in 'combo' vector
      for(Index ic = 0; ic < ncomb; ++ic) {
        std::vector<Index> tcombo = index_to_kcombination(ic, k);

        std::vector<Index> combo(tcombo.rbegin(), tcombo.rend());

        // Consider each permutation of the k-combination elements, which specifies the 'direction' in which to maximize the composition
        std::vector<Index> priority;
        for(Index c : combo)
          priority.push_back(compon[c]);

        //std::cout << "combo: " << combo << "\n";
        do {
          //std::cout << "   priority: " << priority << "\n";
          // Maximize composition in direction specified by the current 'priority'
          Eigen::VectorXi tend(Eigen::VectorXi::Zero(dim));
          for(auto const &sublat : sublats) {
            for(Index i : priority) {
              if(sublat.first.count(i)) {
                tend[i] += sublat.second; //increment by multiplicity of the sublattice
                break;
              }
            }
          }
          // Keep unique extrema
          if(!contains(result, tend))
            result.push_back(tend);

          //std::cout << "tend_members.size(): " << tend_members.size() << "\n";

        }
        while(next_permutation(priority.begin(), priority.end()));   //repeat the above for all permutations of priority
      }

      //std::cout << "sub_max:\n";
      //for(auto const & el : result)
      //std::cout << el.transpose() << "\n";
      //std::cout << "\n";
      return result;
    }

    //*************************************************************

    static std::map<std::set<Index>, std::map<std::set<Index>, Index> > _chemical_subsystems(ParamComposition::AllowedOccupants const &_allowed_occs) {
      std::map<std::set<Index>, std::map<std::set<Index>, Index> > result;

      std::vector<std::string> compon = ParamComposition::string_components(_allowed_occs);

      // Convert _allowed_occs to occ_map, which has a pair for each unique sublattice, consisting of the indices of its allowed components and its multiplicity
      std::map<std::set<Index>, Index> occ_map;
      for(auto const &list : _allowed_occs) {
        std::set<Index> tocc;
        for(auto const &occ : list)
          tocc.insert(find_index(compon, occ));
        auto it = occ_map.find(tocc);
        if(it != occ_map.end())
          ++(it->second);
        else
          occ_map[tocc] = 1;
      }

      //std::cout << "occ_map: \n";
      //for(auto const& occ : occ_map){
      //for(Index i : occ.first){
      //  std::cout << " " << i;
      //}
      //std::cout << ": " << occ.second << "\n";
      //}

      // Use flooding algorithm to partition unique sublattices into chemically independent subsystems

      // sublattices that have been partitioned
      std::set<Index> visited;

      // A chemical subsystem is a graph with species as nodes. Each sublattice is a fully connected subgraph.
      // Traverse the edges implied by this assumption to find the nodes belong to the connected subgraph associated with each species.
      // Loop over all species and record the sublattices that generate each subgraph thus obtained
      for(Index i = 0; i < compon.size() && visited.size() < compon.size(); ++i) {

        // Queue of nodes from which to continue search
        std::set<Index> q({i});

        // Set of nodes in the current subgraph
        std::set<Index> s;

        if(visited.count(i))
          continue;

        // tmap will hold the sublattices for the independent subsystem
        std::map<std::set<Index>, Index> tmap;

        while(!q.empty()) {
          Index j = *q.begin();
          q.erase(q.begin());
          if(!visited.count(j)) {
            visited.insert(j);
            // Loop over sublattices
            for(auto const &list : occ_map) {
              // If searched species 'j' is allowed at this sublattice, add its allowed species to the queue and to the set of connected nodes
              // Add the sublattice to the independent subsystem
              if(list.first.count(j)) {
                q.insert(list.first.begin(), list.first.end());
                s.insert(list.first.begin(), list.first.end());
                tmap.emplace(list);
              }
            }
          }
        }
        result[s] = tmap;
      }

      //std::cout << "subsystems result: \n";
      //for(auto const& sub : result){
      //for(Index i : sub.first){
      //  std::cout << " " << i;
      //}
      //std::cout << " ::\n";
      //for(auto const& occ : sub.second){
      //  for(Index i : occ.first){
      //    std::cout << " " << i;
      //  }
      //  std::cout << ": " << occ.second << "\n";
      //}
      //std::cout << "\n";
      //}


      return result;
    }

  }//\end namespace Local

  //*************************************************************

  std::vector<std::string> ParamComposition::string_components(ParamComposition::AllowedOccupants const &_allowed_occs) {
    std::vector<std::string> result;
    for(auto const &site : _allowed_occs) {
      for(auto const &occ : site) {
        if(!contains(result, occ)) {
          result.push_back(occ);
        }
      }
    }
    return result;
  }

  //---------------------------------------------------------------------------

  ParamComposition::ParamComposition(ParamComposition::AllowedOccupants _allowed_occs)
    : m_allowed_occs(std::move(_allowed_occs)),
      m_prim_end_members(0, 0) {
    m_components = string_components(allowed_occs());
    m_comp.resize(2);
    m_comp[0].resize(0, 0);
    m_comp[1].resize(0, 0);
    m_origin.resize(0);
    m_rank_of_space = -1;
  }

  ParamComposition::ParamComposition(ParamComposition::AllowedOccupants _allowed_occs,
                                     const Eigen::MatrixXd &transf_mat,
                                     const Eigen::VectorXd &_origin,
                                     const int &_rank_of_space,
                                     const int &COMP_TYPE) :
    m_allowed_occs(std::move(_allowed_occs)) {
    m_components = string_components(allowed_occs());
    m_rank_of_space = _rank_of_space;
    m_origin = _origin;
    m_comp.resize(2);
    if(COMP_TYPE == PARAM_COMP) {
      m_comp[PARAM_COMP] = transf_mat;
      m_comp[NUMBER_ATOMS] = m_comp[PARAM_COMP].inverse();
    }
    else if(COMP_TYPE == NUMBER_ATOMS) {
      m_comp[NUMBER_ATOMS] = transf_mat;
      m_comp[PARAM_COMP] = m_comp[NUMBER_ATOMS].inverse();
    }
    calc_spanning_end_members();

  };

  //*************************************************************
  /*   GENERATE_END_MEMBERS

       End members are generated by assigning priority values to each
       component. Based on the priority value, the number of atoms for
       every component is maximized. The routine then iterates thorough
       all possible permutations of the priority values to generate all
       possible end members
  */
  //*************************************************************

  void ParamComposition::generate_prim_end_members() {
    std::map<std::set<Index>, std::map<std::set<Index>, Index> > subsystems = Local::_chemical_subsystems(allowed_occs());

    std::vector<Eigen::VectorXi> tresult(1, Eigen::VectorXi::Zero(components().size()));

    for(auto const &subsystem : subsystems) {
      std::vector<Eigen::VectorXi> tsubs = Local::_sub_extremes(subsystem.first, subsystem.second, components().size());
      std::vector<Eigen::VectorXi> tresult2;
      tresult2.reserve(tresult.size()*tsubs.size());
      for(auto const &v1 : tresult) {
        for(auto const &v2 : tsubs) {
          tresult2.push_back(v1 + v2);
        }
      }
      std::swap(tresult, tresult2);
    }

    std::map<Index, std::vector<Eigen::VectorXi > > tsort;
    for(auto const &v : tresult) {
      tsort[v.squaredNorm()].push_back(v);
    }


    //Store tsort as an Eigen::MatrixXd
    //makes it easier to find the rank of the space
    m_prim_end_members.resize(tresult.size(), components().size());
    Index l = 0;
    for(auto it = tsort.rbegin(); it != tsort.rend(); ++it) {
      //std::cout << "sqnorm " << it->first << " has " << it->second.size() << "\n";
      for(auto const &el : it->second) {
        m_prim_end_members.row(l++) = el.cast<double>().transpose();
      }
    }


  }

  //---------------------------------------------------------------------------

  //*********************************************************************
  /* GENERATE_COMPOSITION_SPACE

     generates all possible composition axes that result in positive
     values of the parametric composition.
     ALGORITHM:
       - start by finding the rank of the space that user has defined
         in the PRIM
       - pick one of the end_members as the origin. To enumerate all
         possible axes, we loop through all possible end_members
       - (rank-1) number of end_members are picked as spanning end
         members from the remaining list of end_members that we get
         from the PRIM
       - A composition object is calculated that is then used to
         calculated the parametric Composition given the current
         choice of end members and origin. If it results in positive
         numbers for all the end members that are listed for the PRIM,
         this set of (origin,spanning end members) is pushed back onto
         the allowed list of composition axes
       - The process is repeated for all such unique combinations of
         (origin, spanning end members)
   */
  //*********************************************************************


  void ParamComposition::generate_composition_space(bool verbose) {
    //Eigen object to do the QR decomposition of the list of prim_end_members
    Eigen::FullPivHouseholderQR<Eigen::MatrixXd> qr(m_prim_end_members);
    Eigen::VectorXd torigin; //temp origin
    Eigen::VectorXd test_comp;

    // If there is already a set of enumerated spaces for this
    // Composition object
    if(m_allowed_list.size() != 0) {
      std::cerr << "WARNING in ParamComposition::generate_composition_space, your allowed_list is non-empty. If you are not careful,you may end up with repeats of allowed composition axes" << std::endl;
    }

    // calculate the rank of the space.
    // NOTE to the wise: # of spanning members = rank-1
    m_rank_of_space = qr.rank();
    if(verbose)
      std::cout << "Rank of space : " << m_rank_of_space << std::endl;

    // Count over K-combinations of end members, for each K-combination, select one as origin and the rest as axes
    Index K = m_rank_of_space;
    Index N = m_prim_end_members.rows();
    Index ncomb = nchoosek(N, K);
    Eigen::MatrixXd tmembers(m_prim_end_members.cols(), K);

    // Combination is stored in 'combo' vector
    std::vector<Index> combo;

    //set of spanning end members
    std::vector< Eigen::VectorXd > tspanning;
    for(Index c = 0; c < ncomb; ++c) {
      combo = index_to_kcombination(c, K);
      for(Index i = 0; i < K; ++i) {
        tmembers.col(i) = m_prim_end_members.row(combo[i]).transpose();
      }
      if(c > 100000 && m_allowed_list.size() > 0)
        break;
      if(qr.compute(tmembers).rank() < K)
        continue;

      if(verbose) {
        double p = round(1e4 * double(c) / double(ncomb)) / 100.;
        std::cout << "combo #" << c << ": " << combo << std::endl;
        std::cout << "Found " << m_allowed_list.size() << std::endl;
        std::cout << p << "% complete. Considering possible origins for set \n" << tmembers << std::endl;
      }

      //loop through all the end members and start them off as the origin
      for(Index i_o = 0; i_o < K; ++i_o) {
        torigin = tmembers.col(i_o);

        if(verbose)
          std::cout << "The origin is: " << torigin.transpose() << "\n---" << std::endl;

        //after picking origin, see if remaining end members span the space
        //end members to span the space. Only keep those combinations
        //that result in positive compositions
        tspanning.clear();
        for(Index j = 0; j < K; ++j) {

          if(j != i_o) {
            tspanning.push_back(tmembers.col(j) - torigin);
          }
        }

        //initialize a ParamComposition object with these spanning vectors and origin
        ParamComposition tcomp = calc_composition_object(torigin, tspanning);

        //loop through end members and see what the values of the compositions works out to
        if(verbose)
          std::cout << "Calculated compositions:" << std::endl;

        //flag to test for positive composition axes
        bool is_positive = true;

        for(Index j = 0; is_positive && j < m_prim_end_members.rows(); j++) {
          // Calculates the composition value given the previously
          // initialized Composition object
          test_comp = tcomp.calc(m_prim_end_members.row(j), NUMBER_ATOMS);
          if(verbose)
            std::cout << m_prim_end_members.row(j) << "  :  " << test_comp.transpose() << std::endl;

          for(Index k = 0; k < test_comp.size(); k++) {
            //if the calculated parametric composition value is either
            //less than 0 or is nan(this will occur if the current set
            //of origin and spanning end members form a subspace in
            //the composition space for this PRIM)
            if(test_comp(k) < -TOL) {
              is_positive = false;
              break;
            }
          }
        }
        if(is_positive) {
          //push back this composition object, its good!
          m_allowed_list.push_back(tcomp);
        }

      }
    }
    if(verbose)
      std::cout << "                                                                                                                          " << std::endl;

  }

  //---------------------------------------------------------------------------

  //*************************************************************
  //PRINT Routines
  void ParamComposition::print_composition_formula(std::ostream &stream, const int &stream_width) const {
    int composition_var = (int)'a';
    std::stringstream tstr;
    for(Index i = 0; i < components().size(); i++) {
      bool first_char = true;
      tstr << components()[i] << "(";
      if(!almost_zero(m_origin(i))) {
        first_char = false;
        tstr << m_origin(i);
      }
      for(int j = 0; j + 1 < (m_rank_of_space); j++) {
        if(almost_zero(m_comp[PARAM_COMP](i, j))) {
          continue;
        }
        if(almost_zero(m_comp[PARAM_COMP](i, j) - 1)) {
          if(first_char) {
            tstr << (char)(composition_var + j);
            first_char = false;
          }
          else {
            tstr << '+' << (char)(composition_var + j);
          }
        }
        else if(almost_zero(m_comp[PARAM_COMP](i, j) + 1)) {
          tstr << '-' << (char)(composition_var + j);
          first_char = false;
        }
        else {
          stream << m_comp[PARAM_COMP](i, j) << (char)(composition_var + j);
          first_char = false;
        }
      }
      tstr << ")";
    }

    stream << tstr.str().c_str();

    return;
  }

  //---------------------------------------------------------------------------

  void ParamComposition::print_member_formula(const Eigen::VectorXd &member, std::ostream &stream, const int &stream_width) const {
    std::stringstream tstr;
    for(EigenIndex i = 0; i < member.size(); i++) {
      if(almost_zero(member(i))) {
        continue;
      }
      if(almost_zero(member(i) - 1)) {
        tstr << components()[i];
      }
      else {
        tstr << components()[i] << int(member(i));
      }
    }
    stream << std::setw(stream_width) << tstr.str().c_str();
  }

  //---------------------------------------------------------------------------

  void ParamComposition::print_origin_formula(std::ostream &stream, const int &stream_width) const {
    print_member_formula(m_origin, stream, stream_width);
  }

  //---------------------------------------------------------------------------

  void ParamComposition::print_composition_axes(std::ostream &stream) const {
    stream << "Number of choices of composition axes: " << m_allowed_list.size() << std::endl;
    //Print Header: ORIGIN, <COMPOUNDS AT THE ENDS OF DIFFERENT AXES> , GENERAL FORMULA

    stream << std::setw(10) << "INDEX";
    stream << std::setw(10) << "ORIGIN";
    for(int i = 0; i + 1 < m_rank_of_space; i++) {
      stream << std::setw(10) << (char)((int)'a' + i);
    }
    stream << "    ";
    stream << "GENERAL FORMULA";
    stream << std::endl;

    stream << std::setw(10) << "  ---";
    stream << std::setw(10) << "  ---";
    for(int i = 0; i < (m_rank_of_space - 1); i++) {
      stream << std::setw(10) << "  ---";
    }
    stream << "    ";
    stream << "---" << std::endl;

    for(Index i = 0; i < m_allowed_list.size(); i++) {
      stream << std::setw(10) << i;
      m_allowed_list[i].print_origin_formula(stream, 10);
      std::vector< Eigen::VectorXd > allowed_spanning_end_members;
      allowed_spanning_end_members = m_allowed_list[i].spanning_end_members();
      for(Index j = 0; j < allowed_spanning_end_members.size(); j++) {
        print_member_formula(allowed_spanning_end_members[j], stream, 10);
      }
      stream << "    ";
      m_allowed_list[i].print_composition_formula(stream, 20);
      stream << std::endl;
    }
    //        print_end_member_formula(1,std::cout,10);
    // for(Index i=0;i<allowed_list.size();i++){
    //     allowed_list[i].print_composition_formula(std::cout);
    //     std::cout<<std::endl;
    // }
  }

  //---------------------------------------------------------------------------

  void ParamComposition::print_curr_composition_axes(std::ostream &stream) const {
    //Print Header: ORIGIN, <COMPOUNDS AT THE ENDS OF DIFFERENT AXES> , GENERAL FORMULA

    stream << std::setw(20) << "ORIGIN";
    for(int i = 0; i < (m_rank_of_space - 1); i++) {
      stream << std::setw(10) << (char)((int)'a' + i);
    }
    stream << "    ";
    stream << "GENERAL FORMULA";
    stream << std::endl;

    stream << std::setw(20) << "  ---";
    for(int i = 0; i < (m_rank_of_space - 1); i++) {
      stream << std::setw(10) << "  ---";
    }
    stream << "    ";
    stream << "---" << std::endl;

    print_origin_formula(stream, 20);
    std::vector< Eigen::VectorXd > allowed_spanning_end_members;
    allowed_spanning_end_members = spanning_end_members();
    for(int j = 0; j < allowed_spanning_end_members.size(); j++) {
      print_member_formula(allowed_spanning_end_members[j], stream, 10);
    }
    stream << "    ";
    print_composition_formula(stream, 20);
    stream << std::endl;

  }


  //*************************************************************
  //CALC ROUTINES

  //*************************************************************
  /*
     CALC

     To calculate the composition AFTER having set the origin and
     spanning end members for the class, you need to pass it the
     "given" values i.e. either the parametric composition or the
     number of atoms per primClex unit.
     If you want the PARAM_COMP given NUMBER_ATOMS set the mode to
     PARAM_COMP.
     If you eant the NUMBER_ATOMS given PARAM_COMP set the mode to
     NUMBER_ATOMS.
     i.e. set the mode to whatever is the quantity that you are giving
     the class
   */
  //*************************************************************

  Eigen::VectorXd ParamComposition::calc(const Eigen::VectorXd &tcomp, const int &MODE) {
    if(MODE == PARAM_COMP) {
      return m_origin + m_comp[PARAM_COMP] * tcomp;
    }
    else {
      return (m_comp[NUMBER_ATOMS] * (tcomp - m_origin)).head(m_rank_of_space - 1);
    }
  }

  //*************************************************************

  Eigen::VectorXd ParamComposition::calc_param_composition(const Eigen::VectorXd &num_atoms_per_prim) const {
    return (m_comp[NUMBER_ATOMS] * (num_atoms_per_prim - m_origin)).head(m_rank_of_space - 1);
  }

  //*************************************************************

  Eigen::VectorXd ParamComposition::calc_num_atoms(const Eigen::VectorXd &param_composition) const {
    return m_origin + m_comp[PARAM_COMP] * param_composition;
  }

  //*************************************************************

  //Given an origin and spanning vectors, returns a ParamComposition object that points to the same Prim as (*this)
  ParamComposition ParamComposition::calc_composition_object(const Eigen::VectorXd &torigin, const std::vector< Eigen::VectorXd> tspanning) {
    //holds the temporary transformation matrix that is going to be
    //used to initialize the new composition object
    Eigen::MatrixXd tmat;
    //if(!tspanning.empty() && tspanning[0].size() != components().size()) {
    //std::cerr << "ERROR in ParamComposition::calc_composition_object the spanning vectors are not as long as the number of ";
    //std::cerr << "components in this system. I'm confused and recommend you quit and try again. However, not going to force quit\n";
    //}
    tmat.setIdentity(components().size(), components().size());
    //copy the spanning vectors into tmat
    for(Index i = 0; i < tspanning.size(); i++) {
      tmat.col(i) = tspanning[i];
    }
    //generate an orthogonal set if there aren't as many spanning
    //vectors as the number of components in the system
    if(tspanning.size() < (components().size())) {
      Eigen::FullPivHouseholderQR<Eigen::MatrixXd> qr(tmat.leftCols(tspanning.size()));
      Eigen::MatrixXd torthogonal = qr.matrixQ();
      tmat.rightCols((components().size() - tspanning.size())) = torthogonal.rightCols((components().size() - tspanning.size()));
      //copy the orthogonalized vectors into tmat. We now have a
      //complete spanning set in this component space
    }
    return ParamComposition(allowed_occs(), tmat, torigin, m_rank_of_space, PARAM_COMP);
  }

  //assuming that you have filled in the prim_end_members and the
  //origin. This fills up the transformation matrices
  void ParamComposition::calc_transformation_matrices() {
    Eigen::MatrixXd tmat;
    tmat.resize(components().size(), components().size());
    for(Index i = 0; i < m_spanning_end_members.size(); i++) {
      tmat.col(i) = m_spanning_end_members[i] - m_origin;
    }
    if(m_spanning_end_members.size() < components().size()) {
      Eigen::FullPivHouseholderQR<Eigen::MatrixXd> qr(tmat.leftCols(m_spanning_end_members.size()));
      Eigen::MatrixXd torthogonal = qr.matrixQ();
      tmat.rightCols((components().size() - m_spanning_end_members.size())) = torthogonal.rightCols((components().size() - m_spanning_end_members.size()));
      //copy the orthogonalized vectors into tmat. We now have a
      //complete spanning set in this component space
    }
    if(m_comp.size() != 2) {
      m_comp.resize(2);
    }
    m_comp[PARAM_COMP] = tmat;
    m_comp[NUMBER_ATOMS] = m_comp[PARAM_COMP].inverse();
  }

  //**************************************************************
  /*SPANNING END MEMBERS

    Returns an std::vector< Eigen::VectorXi > that contain the spanning end
    members listed in the same order as they occur in the
    transformation matrix.

  */
  void ParamComposition::calc_spanning_end_members() {
    std::vector< Eigen::VectorXd > tspan_end;
    if(m_rank_of_space <= 0) {
      std::cerr << "WARNING something is wrong in ParamComposition::spanning_end_members. The rank_of_space in the ParamComposition object is <=0. I do not know how to calculate the end_members in such a space" << std::endl;
      m_spanning_end_members = tspan_end;
      return;
    }
    for(int i = 0; i + 1 < m_rank_of_space; i++) {
      tspan_end.push_back((m_comp[PARAM_COMP].col(i) + m_origin));
    }
    m_spanning_end_members = tspan_end;
  }

  //*************************************************************
  //MISCELLANEOUS

  void ParamComposition::select_composition_axes(const Index &choice) {

    //printing the data from the 'choice'
    // std::cout<<"This is the data from choice: "<<std::endl;
    // allowed_list[choice].print_composition_matrices(std::cout);
    if(m_allowed_list.size() < choice + 1) {
      std::cerr << "ERROR in ParamComposition::select_composition_axes. Your value of choice is outside the range of allowed_list" << std::endl;
      exit(666);
    }
    m_comp = m_allowed_list[choice].comp();
    //    components = allowed_list[choice].components();
    m_origin = m_allowed_list[choice].origin();
    m_rank_of_space = m_allowed_list[choice].rank_of_space();
    m_spanning_end_members = m_allowed_list[choice].spanning_end_members();
    //    print_composition_matrices(std::cout);
  }


  //*************************************************************
  //ACCESSORS

  std::string ParamComposition::composition_formula() const {
    std::stringstream ss;
    print_composition_formula(ss, 20);
    return ss.str();
  }


}

