#ifndef BP_Vasp_CC
#define BP_Vasp_CC

#include "casm/BP_C++/BP_Vasp.hh"


namespace BP {

  void BP_POSCAR_class::read_POSCAR(std::string poscar) {
    reset();

    BP_Vec<double> d_list;
    BP_Vec<std::string> s_list;
    BP_Vec<std::string> atom_types_listed_after_positions;
    int i, j;

    BP_Parse infile(poscar);

    header_line = infile.getline();
    scl = BP::stod(infile.getline());


    // get supercell vectors
    for(i = 0; i < 3; i++) {
      d_list = infile.getline_double();
      lat_vectors.add(cart_coord_class(scl * d_list[0], scl * d_list[1], scl * d_list[2]));
    }

    //scl = 1.0;

    //// type/num of atoms
    s_list = infile.getline_string();

    atom_types_exist = 0;
    if(isalpha(s_list[0][0])) {
      // get type of atoms
      atom_types_exist = 1;

      atom_type = s_list;

      s_list = infile.getline_string();

    }

    num_atoms = stoi(s_list);

    set_cum_num_atoms();

    // get Direct / Cartesian
    s_list = infile.getline_string();

    if(s_list[0][0] == 'S' || s_list[0][0] == 's') {
      selective_dynamics = true;
      s_list = infile.getline_string();

    }
    else
      selective_dynamics = false;


    if(s_list[0][0] == 'C' || s_list[0][0] == 'c' || s_list[0][0] == 'K' || s_list[0][0] == 'k') {
      // cartesian

      direct = false;

      for(i = 0; i < sum_num_atoms; i++) {
        s_list = infile.getline_string();
        atom_pos_cart.add(cart_coord_class(BP::stod(s_list[0]), BP::stod(s_list[1]), BP::stod(s_list[2])));

        if(selective_dynamics) {
          sel_dynamics.add();
          sel_dynamics.last().add(next_bool(s_list[3]));
          sel_dynamics.last().add(next_bool(s_list[4]));
          sel_dynamics.last().add(next_bool(s_list[5]));
        }
        else if(s_list.size() > 3 && atom_types_exist == false) {
          // add atom type if listed
          atom_types_listed_after_positions.add(s_list[3]);
          set_write_all_atom_types(true);
        }

        //poscar_index.add(i+1);
      }

    }
    else {
      // direct

      direct = true;

      for(i = 0; i < sum_num_atoms; i++) {
        s_list = infile.getline_string();
        atom_pos_frac.add(frac_coord_class(BP::stod(s_list[0]), BP::stod(s_list[1]), BP::stod(s_list[2])));

        if(selective_dynamics) {
          sel_dynamics.add();
          sel_dynamics.last().add(next_bool(s_list[3]));
          sel_dynamics.last().add(next_bool(s_list[4]));
          sel_dynamics.last().add(next_bool(s_list[5]));
        }
        else if(s_list.size() > 3  && atom_types_exist == false) {
          // add atom type if listed
          atom_types_listed_after_positions.add(s_list[3]);
          set_write_all_atom_types(true);
        }

        //poscar_index.add(i+1);
      }

    }



    // if atom types were listed after the atom coordinates,
    //   need to set cum_num_atoms and atom_type lists
    if(atom_types_listed_after_positions.size() > 0) {
      for(i = 0; i < num_atoms.size(); i++) {
        atom_type.add(atom_types_listed_after_positions[0]);

        for(j = 0; j < num_atoms[i]; j++)
          atom_types_listed_after_positions.ordered_remove(0);

      }

      atom_types_exist = true;

    }

    //std::cout << "finish read_POSCAR" << std::endl;
  };

  void BP_POSCAR_class::write_POSCAR(std::string poscar) {
    BP_Write outfile(poscar);
    outfile.newfile();

    write_POSCAR(outfile.get_ostream());

  };

  void BP_POSCAR_class::write_POSCAR(std::ostream &sout) {

    int w = 20;
    int p = 16;

    sout << header_line << std::endl;
    sout << std::fixed << std::setw(w) << std::setprecision(p) << scl << std::endl;

    for(int i = 0; i < 3; i++) {
      for(int j = 0; j < 3; j++) {
        sout << std::fixed << std::setw(w) << std::setprecision(p) << lat_vectors[i][j] / scl << " ";
      }
      sout << "\n";
    }

    if(atom_types_exist == 1 && format_version == 5) {
      for(int i = 0; i < atom_type.size(); i++)
        sout << atom_type[i] << " ";
      sout << "\n";
    }

    for(int i = 0; i < num_atoms.size(); i++)
      sout << num_atoms[i] << " ";
    sout << "\n";

    if(selective_dynamics)
      sout << "Selective Dynamics\n";

    if(direct) {
      sout << "Direct" << std::endl;

      for(int i = 0; i < atom_pos_frac.size(); i++) {
        for(int j = 0; j < 3; j++)
          sout << std::fixed << std::setw(w) << std::setprecision(p) << atom_pos_frac[i][j] << " ";

        if(selective_dynamics) {
          for(int j = 0; j < sel_dynamics[i].size(); j++) {
            if(sel_dynamics[i][j])
              sout << "T ";
            else
              sout << "F ";
          }
        }

        if(write_all_atom_types && atom_types_exist)
          sout << get_type(i) << " ";

        sout << "\n";
      }
    }
    else {
      sout << "Cartesian" << std::endl;

      for(int i = 0; i < atom_pos_cart.size(); i++) {
        for(int j = 0; j < 3; j++)
          sout << std::fixed << std::setw(w) << std::setprecision(p) << atom_pos_cart[i][j] << " ";

        if(selective_dynamics) {
          for(int j = 0; j < sel_dynamics[i].size(); j++) {
            if(sel_dynamics[i][j])
              sout << "T ";
            else
              sout << "F ";
          }
        }

        if(write_all_atom_types && atom_types_exist)
          sout << get_type(i) << " ";


        sout << "\n";
      }
    }

  };

}

#endif // BP_Vasp_CC

