#ifndef BP_Vasp_HH
#define BP_Vasp_HH

#include <iostream>
#include <math.h>
#include <iomanip>
#include "casm/BP_C++/BP_Parse.hh"
#include "casm/BP_C++/BP_Vec.hh"
#include "casm/BP_C++/BP_coords.hh"
#include "casm/external/Eigen/Dense"

namespace BP {

  class BP_POSCAR_class {
  private:

    double scl;
    BP_Vec< cart_coord_class > lat_vectors;		// periodic lattice vectors: lat_vectors[ vector][ length x,y,z]
    BP_Vec< frac_coord_class > atom_pos_frac;		//position of atoms in direct coordinates: atom_pos[atom][position vec1, vec2, vec3]   atoms are added in type order
    BP_Vec< cart_coord_class > atom_pos_cart;		//position of atoms in direct coordinates: atom_pos[atom][position vec1, vec2, vec3]
    //BP_Vec<int> poscar_index;				//order in poscar (starts with 1:  1,2,3,4,...)
    int sum_num_atoms;						//total number of atoms
    bool atom_types_exist;					// whether or not atom names included in POSCAR; VASP version 4: atom_types_exist == 0, VASP version 5: atom_types_exist == 1
    BP_Vec<std::string> atom_type;				// atom type (list of each type of atom in the POSCAR, not the type of each atom)
    BP_Vec<int> num_atoms;					// # of atoms of each type
    BP_Vec<int> cum_num_atoms;				// cumulative # of atoms of each type up to this one
    BP_Vec< BP_Vec<bool> > sel_dynamics;

    //int nbands;
    //bool pro_exists;
    //float ***pro;
    //double *Eband;
    //double *band_occ;
    //double *band_tot;

    bool direct;
    std::string header_line;
    bool selective_dynamics;
    int format_version;			// 4: don't write atom_type line; 5: do write atom_type line
    int write_all_atom_types;	// write atom_type at end of line for each atom position

  public:

    void read_POSCAR(std::string poscar);
    void write_POSCAR(std::ostream &sout);
    void write_POSCAR(std::string poscar);

    BP_POSCAR_class() {
      reset();

    };

    BP_POSCAR_class(std::string poscar) {
      reset();
      read_POSCAR(poscar);
    };

    ~BP_POSCAR_class() {


    };

    void reset() {
      direct = true;
      header_line = "";
      selective_dynamics = false;
      sum_num_atoms = 0;
      scl = 1.0;
      atom_types_exist = false;
      write_all_atom_types = false;
      format_version = 5;
      atom_pos_frac.clear();
      atom_pos_cart.clear();
      //poscar_index.clear();
      num_atoms.clear();
      cum_num_atoms.clear();
      atom_type.clear();
    };

    //////////////

    void set_format_version(int i) {
      if(i == 4)
        format_version = 4;
      else if(i == 5)
        format_version = 5;

    }

    void set_write_all_atom_types(bool b) {
      write_all_atom_types = b;
    }

    //////////////

    void set_lat_vectors(double i1, const BP_Vec< cart_coord_class > &i2) {
      scl = i1;
      lat_vectors = i2;
    };

    void get_lat_vectors(double &d, BP_Vec< cart_coord_class > &vec) {
      d = scl;
      vec = lat_vectors;
      for(int i = 0; i < 3; i++)
        for(int j = 0; j < 3; j++)
          vec[i][j] /= scl;

    };

    BP_Vec< cart_coord_class > get_lat_vectors() const {
      return lat_vectors;
    };

    ///		Returns prim_class::pr_v as a matrix with columns as pr_v
    ///
    BP_Vec< BP_Vec<double> > get_lat_vectors_inv_matrix() const {
      int i, j;
      // supercell vectors
      Eigen::Matrix3d M, Minv;
      for(i = 0; i < 3; i++) {
        for(j = 0; j < 3; j++) {
          M(j, i) = lat_vectors[i][j];
        }
      }

      // inverse supercell vectors
      Minv = M.inverse();

      BP_Vec< BP_Vec<double> > m = BP_Vec< BP_Vec<double> >(3, BP_Vec<double>(3, 0.0));

      for(i = 0; i < 3; i++)
        for(j = 0; j < 3; j++)
          m[i][j] = Minv(i, j);

      return m;
    };

    ///		Returns prim_class::pr_v as a matrix with columns as pr_v
    ///
    BP_Vec< BP_Vec<double> > get_lat_vectors_matrix() const {
      BP_Vec< BP_Vec<double> > m = BP_Vec< BP_Vec<double> >(3, BP_Vec<double>(3, 0.0));

      for(int i = 0; i < 3; i++)
        for(int j = 0; j < 3; j++)
          m[i][j] = lat_vectors[j][i];

      return m;
    };

    //////////////

    std::string get_header_line() {
      return header_line;
    };

    void set_header_line(std::string s) {
      header_line = s;
    };

    ////////////// get total number of atoms, and acces num_atoms

    unsigned long int size() const {
      return sum_num_atoms;
    };

    unsigned long int get_num_atoms() {
      return sum_num_atoms;
    };

    unsigned long int get_num_atoms(int i) {
      return num_atoms[i];
    };

    unsigned long int get_num_atoms(std::string type) {
      if(atom_types_exist) {
        unsigned long int type_index;
        if(atom_type.find_first(type, type_index)) {
          return num_atoms[type_index];
        }
        else
          return 0;
      }
      else
        return 0;
    };

    ////////////// access atom_type vector (list of which types of atoms are in the list)

    BP_Vec<std::string> get_atom_type() const {
      return atom_type;
    };

    std::string get_atom_type(int i) {
      return atom_type[i];
    };

    unsigned long int atom_type_size() {
      return num_atoms.size();
    };

    bool get_atom_types_exist() {
      return atom_types_exist;
    };

    void set_atom_type(BP_Vec<std::string> s_list) {
      if(num_atoms.size() == s_list.size()) {
        atom_type = s_list;
        atom_types_exist = true;
      }
    };

    void set_atom_type(unsigned long int i, std::string s) {
      if(atom_types_exist) {
        atom_type[i] = s;
        atom_types_exist = true;
      }
    };

    /////////////// given index of atom, return that atom's type or type_index

    unsigned long int get_type_index(unsigned long int i) {
      unsigned long int j = 0;

      while(cum_num_atoms[j] <= i) {
        j++;
      }

      return j;
    };

    std::string get_type(unsigned long int i) const {
      unsigned long int j = 0;


      while(cum_num_atoms[j] <= i) {
        j++;

      }

      return atom_type[j];
    };

    //////////////

    void set_selective_dynamics(bool i1) {
      if(i1) {
        if(!selective_dynamics) {
          selective_dynamics = true;
          sel_dynamics = BP_Vec< BP_Vec<bool> >(sum_num_atoms, BP_Vec<bool>(3, true));
        }
      }
      else {
        selective_dynamics = false;
        sel_dynamics.clear();
      }
    };

    void set_sel_dynamics(int i, BP_Vec<bool> sel) {
      if(!selective_dynamics) {
        selective_dynamics = true;
        sel_dynamics = BP_Vec< BP_Vec<bool> >(sum_num_atoms, BP_Vec<bool>(3, true));
      }

      sel_dynamics[i] = sel;

    };

    void set_sel_dynamics(BP_Vec< BP_Vec<bool> > sel) {
      if(sel.size() != sum_num_atoms)
        return;

      selective_dynamics = true;
      sel_dynamics = sel;

    };

    bool get_selective_dynamics() {
      return selective_dynamics;
    };

    BP_Vec< BP_Vec<bool> > get_sel_dynamics() {
      return sel_dynamics;
    };

    BP_Vec<bool> get_sel_dynamics(unsigned long int  i) {
      return sel_dynamics[i];
    };

    //////////////

    bool get_coord_mode() {
      return direct;
    }

    void set_coord_mode(bool i1) {
      if(i1)
        set_direct();
      else
        set_cart();
    }

    void set_direct() {
      if(!direct) {
        direct = true;

        atom_pos_frac.clear();

        for(unsigned long int i = 0; i < atom_pos_cart.size(); i++) {
          atom_pos_frac.add(cart_to_frac(lat_vectors, atom_pos_cart[i]));

        }

      }
    };

    void set_cart() {
      if(direct) {
        direct = false;

        atom_pos_cart.clear();

        for(unsigned long int  i = 0; i < atom_pos_frac.size(); i++) {
          atom_pos_cart.add(frac_to_cart(lat_vectors, atom_pos_frac[i]));

        }

      }
    };

    bool is_direct() {
      return direct;
    };

    bool is_cart() {
      return !direct;
    };

    ////////////// add an atom, giving frac_coords

    void add(std::string type, frac_coord_class f) {
      if(!direct) {
        cart_coord_class c = frac_to_cart(lat_vectors, f);
        add(type, c);
        return;
      }

      unsigned long int type_index;

      atom_types_exist = true;

      if(!atom_type.find_first(type, type_index)) {
        type_index = atom_type.size();

        atom_type.add(type);
        num_atoms.add(0);

        set_cum_num_atoms();
      }


      add(type_index, f);

    };

    void add(unsigned long int  type_index, frac_coord_class f) {
      if(!direct) {
        cart_coord_class c = frac_to_cart(lat_vectors, f);
        add(type_index, c);
        return;
      }

      if(type_index >= atom_type.size())
        return;

      unsigned long int atom_index = cum_num_atoms[type_index];
      atom_pos_frac.add_in_place(f, atom_index);

      if(selective_dynamics) {
        BP_Vec<bool> sel = BP_Vec<bool>(3, true);
        sel_dynamics.add_in_place(sel, atom_index);
      }

      num_atoms[type_index]++;
      set_cum_num_atoms();
    };

    void add(std::string type, frac_coord_class f, BP_Vec<bool> sel) {
      if(!direct) {
        cart_coord_class c = frac_to_cart(lat_vectors, f);
        add(type, c, sel);
        return;
      }

      unsigned long int type_index;

      atom_types_exist = true;

      if(!atom_type.find_first(type, type_index)) {
        type_index = atom_type.size();

        atom_type.add(type);
        num_atoms.add(0);

        set_cum_num_atoms();
      }


      add(type_index, f, sel);

    };

    void add(unsigned long int  type_index, frac_coord_class f, BP_Vec<bool> sel) {
      if(type_index >= atom_type.size())
        return;

      if(!selective_dynamics)
        return;

      unsigned long int atom_index = cum_num_atoms[type_index];
      atom_pos_frac.add_in_place(f, atom_index);

      sel_dynamics.add_in_place(sel, atom_index);

      num_atoms[type_index]++;
      set_cum_num_atoms();
    };

    ////////////// add an atom, giving cart_coords

    void add(std::string type, cart_coord_class c) {
      if(direct) {
        frac_coord_class f = cart_to_frac(lat_vectors, c);
        add(type, f);
        return;
      }

      unsigned long int type_index;

      atom_types_exist = true;

      if(!atom_type.find_first(type, type_index)) {

        type_index = atom_type.size();

        atom_type.add(type);
        num_atoms.add(0);

        set_cum_num_atoms();

      }


      add(type_index, c);

    };

    void add(unsigned long int  type_index, cart_coord_class c) {
      if(direct) {
        frac_coord_class f = cart_to_frac(lat_vectors, c);
        add(type_index, f);
        return;
      }

      if(type_index >= atom_type.size())
        return;

      unsigned long int atom_index = cum_num_atoms[type_index];
      atom_pos_cart.add_in_place(c, atom_index);

      if(selective_dynamics) {
        BP_Vec<bool> sel = BP_Vec<bool>(3, true);
        sel_dynamics.add_in_place(sel, atom_index);
      }

      num_atoms[type_index]++;
      set_cum_num_atoms();
    };

    void add(std::string type, cart_coord_class c, BP_Vec<bool> sel) {
      if(direct) {
        frac_coord_class f = cart_to_frac(lat_vectors, c);
        add(type, f, sel);
        return;
      }

      unsigned long int type_index;

      atom_types_exist = true;

      if(!atom_type.find_first(type, type_index)) {
        type_index = atom_type.size();

        atom_type.add(type);
        num_atoms.add(0);

        set_cum_num_atoms();
      }


      add(type_index, c, sel);

    };

    void add(unsigned long int  type_index, cart_coord_class c, BP_Vec<bool> sel) {
      if(type_index >= atom_type.size())
        return;

      if(!selective_dynamics)
        return;

      unsigned long int atom_index = cum_num_atoms[type_index];
      atom_pos_cart.add_in_place(c, atom_index);

      sel_dynamics.add_in_place(sel, atom_index);

      num_atoms[type_index]++;
      set_cum_num_atoms();
    };

    ////////////// get atom coordinate

    frac_coord_class frac(unsigned long int i) const {
      return atom_pos_frac[i];
    };

    cart_coord_class cart(unsigned long int i) const {
      return atom_pos_cart[i];
    };

    /////////////// remove atoms

    void remove(unsigned long int  index) {
      if(direct) {
        atom_pos_frac.remove(index);
      }
      else {
        atom_pos_cart.remove(index);
      }

      num_atoms[ get_type_index(index)]--;
      clean_atom_types();
      set_cum_num_atoms();
    }

    void remove_all() {
      atom_pos_frac.clear();
      atom_pos_cart.clear();
      num_atoms.clear();
      atom_type.clear();
      cum_num_atoms.clear();
      sum_num_atoms = 0;

    };

    /////////////// transform atom positions

    void in_cell() {
      bool orig_coord_mode = get_coord_mode();

      set_direct();

      for(unsigned long int i = 0; i < atom_pos_frac.size(); i++)
        atom_pos_frac[i].in_unit();

      set_coord_mode(orig_coord_mode);
    };

    void translate_frac(frac_coord_class f) {
      bool orig_coord_mode = get_coord_mode();

      set_direct();

      for(unsigned long int i = 0; i < atom_pos_frac.size(); i++)
        atom_pos_frac[i] = atom_pos_frac[i] + f;

      set_coord_mode(orig_coord_mode);
    };

    void translate_cart(cart_coord_class c) {
      bool orig_coord_mode = get_coord_mode();

      set_cart();

      for(unsigned long int i = 0; i < atom_pos_cart.size(); i++)
        atom_pos_cart[i] = atom_pos_cart[i] + c;

      set_coord_mode(orig_coord_mode);
    };

    //make_supercell()

    ///////////////

  private:

    void clean_atom_types() {
      unsigned long int  i = 0;
      while(i < atom_type.size()) {
        if(num_atoms[i] == 0) {
          num_atoms.ordered_remove(i);
          atom_type.ordered_remove(i);
        }
        else {
          i++;
        }

      }
    };

    void set_cum_num_atoms() {
      cum_num_atoms.clear();

      if(num_atoms.size() == 0) {
        sum_num_atoms = 0;
      }
      else {
        cum_num_atoms.add(num_atoms[0]);
        for(unsigned long int i = 1; i < num_atoms.size(); i++)
          cum_num_atoms.add(num_atoms[i] + cum_num_atoms[i - 1]);
        sum_num_atoms = cum_num_atoms.last();
      }
    };

  };

}

#endif // BP_Vasp_HH

