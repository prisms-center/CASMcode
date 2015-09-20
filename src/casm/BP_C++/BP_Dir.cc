#ifndef BP_Dir_CC
#define BP_Dir_CC

#include "casm/BP_C++/BP_Dir.hh"


namespace BP {

  /// Constructor, initializes at 'getcwd()'
  ///
  BP_Dir::BP_Dir() {
    char path[1024];
    _joint = "/";

    _abs_name = std::string(getcwd(path, sizeof(path)));
    _name = without_path(_abs_name);  //tokenize( _abs_name, _joint).last();


    generate_files_dirs();
  };

  std::string BP_Dir::without_path(std::string s) const {
    return tokenize(_abs_name, _joint).last();
  };

  bool BP_Dir::is_absolute_path(std::string s) const {
    //return ( s[0] == '/')
    return (s.substr(0, _joint.size()) == _joint);
  };

  void BP_Dir::make_absolute_path(std::string &s) const {
    if(is_absolute_path(s)) return;
    else s = _abs_name + _joint + s;
  };

  std::string BP_Dir::abs_name(std::string s) const {
    make_absolute_path(s);
    return s;
  };

  /// change current directory
  ///
  bool BP_Dir::cd(std::string s) {
    BP_Dir orig_dir = *this;
    bool all_ok = true;

    // set
    if(is_absolute_path(s)) {
      //absolute path

      char path[1024];
      _abs_name = std::string(s);
      _name = tokenize(_abs_name, _joint).last();

    }
    else {
      // relative path
      BP_Vec<std::string> path_list;
      path_list = tokenize(s, _joint);
      for(int i = 0; i < path_list.size(); i++) {
        all_ok &= change_directory(path_list[i]);
        if(!all_ok) break;
      }
    }

    if(all_ok) {
      all_ok &= generate_files_dirs();
    }

    if(!all_ok) {
      *this = orig_dir;
      return false;
    }

    return true;

  };

  /// Make a new directory 's'
  ///
  bool BP_Dir::mkdir(std::string s) {
    // make absolute
    make_absolute_path(s);

    if(exists(s)) return false;

    //S_IRWXU: read, write, execute/search by owner
    //S_IRUSR: read permission, owner
    //S_IWUSR: write permission, owner
    //S_IXUSR: execute/search permission, owner
    //S_IRWXG: read, write, execute/search by group
    //S_IRGRP: read permission, group
    //S_IWGRP: write permission, group
    //S_IXGRP: execute/search permission, group
    //S_IRWXO: read, write, execute/search by others
    //S_IROTH: read permission, others
    //S_IWOTH: write permission, others
    //S_IXOTH: execute/search permission, others
    //S_ISUID: set-user-ID on execution
    //S_ISGID: set-group-ID on execution
    //S_ISVTX: on directories, restricted deletion flag

    // ::mkdir is in global namespace
    ::mkdir(s.c_str(), S_IRWXU | S_IRGRP | S_IXGRP | S_IROTH | S_IXOTH);

    refresh();

    return true;

  };

  /// Rename a file or directory 's1' -> 's2'
  ///
  bool BP_Dir::rename(std::string s1, std::string s2) {
    // translate to absolute paths
    make_absolute_path(s1);
    make_absolute_path(s2);

    //std::cout << "s1: " << s1 << std::endl;
    //std::cout << "s2: " << s2 << std::endl;

    if(::rename(s1.c_str(), s2.c_str()) != 0) {
      //std::cout << "mv not OK" << std::endl;
      return false;
    }
    //std::cout << "mv OK" << std::endl;
    refresh();
    return true;
  };

  /// Copy a file 's1' to 's2'
  ///
  bool BP_Dir::copy(std::string s1, std::string s2) {
    // translate to absolute paths
    make_absolute_path(s1);
    make_absolute_path(s2);

    std::ifstream source;
    std::ofstream destination;

    source.open(s1.c_str(), std::ios::binary);
    destination.open(s2.c_str(), std::ios::binary);

    if(!source) {
      std::cout << "source error" << std::endl;
      return false;
    }

    if(!destination) {
      std::cout << "destination error" << std::endl;
      return false;
    }

    destination << source.rdbuf();

    source.close();
    destination.close();

    refresh();
    return true;
  };

  /// Copy a directory 's1' & all sub-directories/files to 's2'
  ///
  bool BP_Dir::copydir(std::string s1, std::string s2) {
    // translate to absolute paths
    make_absolute_path(s1);
    make_absolute_path(s2);

    // check if s1 is a directory
    if(!is_dir(s1)) return false;

    // check that s2 isn't a file
    if(is_file(s2)) return false;

    bool all_ok = true;
    unsigned long int i;

    // mkdir s2
    // open dir1 at s1 and dir2 at s2
    // copy all files in dir1 to dir2
    // copydir all directories in dir1 to dir2

    // mkdir s2
    if(!exists(s2)) mkdir(s2);

    // open dir1 at s1 and dir2 at s2
    BP_Dir dir1, dir2;
    dir1.cd(s1);
    dir2.cd(s2);

    // copy all the files in dir1 to dir2
    for(i = 0; i < dir1.files().size(); i++) {
      all_ok &= copy(dir1.abs_name() + _joint + dir1.files()[i], dir2.abs_name() + _joint + dir1.files()[i]);
      if(!all_ok) return false;
    }

    // copy all the directories in dir1 to dir2
    for(i = 0; i < dir1.dirs().size(); i++) {
      all_ok &= copydir(dir1.abs_name() + _joint + dir1.dirs()[i], dir2.abs_name() + _joint + dir1.dirs()[i]);
      if(!all_ok) return false;
    }

    refresh();
    return true;
  };

  /// Catenate a files 's1' + 's2' -> write to file s3
  ///
  bool BP_Dir::cat(std::string s1, std::string s2, std::string s3) {
    // translate to absolute paths
    make_absolute_path(s1);
    make_absolute_path(s2);
    make_absolute_path(s3);


    std::ifstream source1, source2;
    std::ofstream destination;

    source1.open(s1.c_str(), std::ios::binary);
    source2.open(s2.c_str(), std::ios::binary);
    destination.open(s3.c_str(), std::ios::binary);

    if(!source1) {
      std::cout << "source1 error" << std::endl;
      return false;
    }

    if(!source2) {
      std::cout << "source2 error" << std::endl;
      return false;
    }

    if(!destination) {
      std::cout << "destination error" << std::endl;
      return false;
    }

    destination << source1.rdbuf() << source2.rdbuf();

    source1.close();
    source2.close();
    destination.close();

    refresh();
    return true;
  };


  // -----Remove files & directories---------------------------
  // All of these only act on files and directories in the current directory (_abs_name)

  /// Remove an empty directory
  ///
  bool BP_Dir::rmdir(std::string s) {
    if(s == "*")
      return rmdir(dirs_list);
    else
      return remove_dir(s);
  }

  /// Remove a list of empty directories
  ///
  bool BP_Dir::rmdir(BP_Vec<std::string> s_list) {
    // try to delete each directory in s_list
    bool all_ok = true;
    for(unsigned long int i = 0; i < s_list.size(); i++) {
      all_ok &= remove_dir(s_list[i]);
      if(!all_ok) return false;
    }

    return all_ok;
  }


  /// Remove a directory & all sub-directories/files
  ///
  bool BP_Dir::recurs_rmdir(std::string s) {
    if(s == "*")
      return recurs_rmdir(dirs_list);
    else
      return recurs_remove_dir(s);
  }

  /// Remove a list of directories & all sub-directories/files
  ///
  bool BP_Dir::recurs_rmdir(BP_Vec<std::string> s_list) {
    // for all directories list in
    // remove directory

    // try to delete each directory in s_list
    bool all_ok = true;
    for(unsigned long int i = 0; i < s_list.size(); i++) {
      all_ok &= recurs_remove_dir(s_list[i]);
      if(!all_ok) return false;
    }

    return all_ok;
  }


  /// Remove a file
  ///
  bool BP_Dir::rm(std::string s) {
    if(s == "*")
      return rm(files_list);
    else
      return remove_file(s);
  }

  /// Remove a list of files
  ///
  bool BP_Dir::rm(BP_Vec<std::string> s_list) {
    // try to delete each file in s_list
    bool all_ok = true;
    for(unsigned long int i = 0; i < s_list.size(); i++) {
      all_ok &= remove_file(s_list[i]);
      if(!all_ok) return false;
    }

    return all_ok;
  }


  /// Print file 's' to stream
  ///
  bool BP_Dir::print(std::string s, std::ostream &sout) const {
    if(!is_file(s)) return false;

    // make absolute
    make_absolute_path(s);

    BP_Parse file(s);
    do {
      sout << file.getline_exact();
    }
    while(!file.eof());
    return true;
  };

  /// Check if 's' is a file
  ///
  bool BP_Dir::is_file(std::string s) const {
    //unsigned long int index;
    //if( files_list.find_first(s, index))
    //	return true;
    //return false;

    if(s.size() == 0) return false;

    // make absolute
    make_absolute_path(s);

    struct stat st;
    if(lstat(s.c_str(), &st) != 0) {
      return false;
    }
    else {
      if(S_ISDIR(st.st_mode)) {
        //std::cout << " DIRECTORY" << std::endl;
        return false;

      }
      else {
        //std::cout << " FILE" << std::endl;
        return true;
      }
    }

  };

  /// Check if 's' is a directory
  ///
  bool BP_Dir::is_dir(std::string s) const {
    //unsigned long int index;
    //if( dirs_list.find_first(s, index))
    //	return true;
    //return false;

    if(s.size() == 0) return false;

    // make absolute
    make_absolute_path(s);

    struct stat st;
    if(lstat(s.c_str(), &st) != 0) {
      return false;
    }
    else {
      if(S_ISDIR(st.st_mode)) {
        //std::cout << " DIRECTORY" << std::endl;
        return true;

      }
      else {
        //std::cout << " FILE" << std::endl;
        return false;
      }
    }

  };

  /// Check if 's' is a file or a directory
  ///
  bool BP_Dir::exists(std::string s) const {
    //if( is_file(s) || is_dir(s))
    //	return true;
    //return false;

    if(s.size() == 0) return false;

    // make absolute
    make_absolute_path(s);

    struct stat st;
    if(lstat(s.c_str(), &st) != 0) {
      return false;
    }
    else {
      return true;
    }
  };

  /// Refresh the list of files and directories in the current directory
  ///
  void BP_Dir::refresh() {
    generate_files_dirs();
  };


  //-----------Private-------------------------------//

  ///		Opens the directory _abs_name and finds all files and directories
  ///			- if it suceeds, returns true
  ///			- if it fails, returns false, and everything is unchanged
  ///
  bool BP_Dir::generate_files_dirs() {
    /// add new directories
    DIR *curr_dir;
    struct dirent *obj;
    struct stat st;
    char path[1024];
    std::string obj_name;
    BP_Vec< std::string> new_files_list;
    BP_Vec< std::string> new_dirs_list;



    // open the curr directory
    curr_dir = opendir(_abs_name.c_str());
    if(!curr_dir) {
      return false;
      //std::cout << "BP_Dir Error in generate_files_dirs, opendir() failed" << std::endl;
      //exit(1);
    }
    errno = 0;

    while((obj = readdir(curr_dir))) {
      obj_name = std::string(obj->d_name);

      lstat((_abs_name + _joint + std::string(obj->d_name)).c_str(), &st);
      if(S_ISDIR(st.st_mode)) {
        //std::cout << " DIRECTORY" << std::endl;
        if(obj_name != "." && obj_name != "..") {
          //std::cout << "add dir: " << obj_name << std::endl;
          new_dirs_list.add(obj_name);
        }
      }
      else {
        //std::cout << " FILE" << std::endl;
        new_files_list.add(obj_name);
      }
    }

    if(errno) {
      return false;
      //std::cout << "BP_Dir Error in generate_files_dirs, readdir() failed" << std::endl;
      //exit(1);
    }

    closedir(curr_dir);

    // clear files and dirs
    files_list = new_files_list;
    dirs_list = new_dirs_list;


    return true;
  };

  bool BP_Dir::change_directory(std::string s) {
    if(s == "..") {
      BP_Vec<std::string> path_list = tokenize(_abs_name, "/");

      if(path_list.size() == 1) {
        _name = "";
        _abs_name = _joint;

        return generate_files_dirs();
      }
      else if(path_list.size() > 1) {
        _name = path_list[ path_list.size() - 2];
        _abs_name = "";
        for(int i = 0; i < path_list.size() - 1; i++)
          _abs_name += _joint + path_list[i];

        return generate_files_dirs();
      }

    }
    else {
      unsigned long int index;
      if(dirs_list.find_first(s, index)) {
        _name = s;
        _abs_name += _joint + s;

        return generate_files_dirs();
      }
    }
    return false;
  };

  bool BP_Dir::remove_dir(std::string s) {
    // directory to remove must be in current directory
    unsigned long int i, index;

    //std::cout << "about to remove directory: " << (_abs_name + "/" + s).c_str() << std::endl;
    //BP_pause();

    // remove a file
    if(dirs_list.find_first(s, index)) {
      if(::rmdir((_abs_name + _joint + s).c_str()) != 0) {
        // error removing the directory
        return false;
      }

      refresh();
      return true;


    }

    // directory not found
    return false;
  }

  bool BP_Dir::recurs_remove_dir(std::string s) {
    // directory to remove must be in current directory
    unsigned long int i, index;
    bool all_ok = true;
    // first cd to s
    // then remove all files
    // then recurs_remove all directories
    // then remove directory

    if(dirs_list.find_first(s, index)) {

      // first cd to s
      all_ok &= cd(s);
      if(!all_ok) return false;

      // then remove all files in that directory
      all_ok &= rm("*");
      if(!all_ok) return false;

      // then recurs_remove all directories
      while(dirs_list.size() > 0) {
        all_ok &= recurs_remove_dir(dirs_list[0]);
        if(!all_ok) return false;
      }
      // then move back up one
      all_ok &= cd("..");
      if(!all_ok) return false;

      // then remove the directory
      all_ok &= rmdir(s);
      if(!all_ok) return false;

    }
    else {
      all_ok = false;
    }

    return all_ok;




  }

  bool BP_Dir::remove_file(std::string s) {
    // file to remove must be in current directory

    unsigned long int i, index;

    std::cout << "about to remove file: " << (_abs_name + _joint + s).c_str() << std::endl;
    //BP_pause();

    // remove a file
    if(files_list.find_first(s, index)) {
      if(std::remove((_abs_name + _joint + s).c_str()) != 0) {
        // error removing the file
        return false;
      }

      refresh();
      return true;


    }

    std::cout << "file not found" << std::endl;

    // file not found
    return false;

  }


  /*
  class BP_Dir_Edge
  {
  	public:
  	bool val;
  };

  // cd, mkdir, get files and directories
  ///< \ingroup BP BP_Dir
  class BP_Dir_Tree : private BP_Graph< BP_Dir, BP_Dir_Edge>
  {
  	private:
  	BP_GVec<BP_Dir>					dir_list;
  	BP_GVec<BP_Dir_Edge>			edge_list;

  	BP_GVec_Member<BP_Dir>			*curr;
  	BP_GVec_Member<BP_Dir>			*ref;

  	public:

  	BP_Dir_Tree();
  	//~BP_Dir_Tree();

  	void cd( std::string s);
  	void cd( BP_GVec_Member<BP_Dir> *dir);
  	void mkdir(std::string s);

  	std::string cwd();
  	BP_GVec_Member<BP_Dir>* cwdp();

  	BP_Vec<std::string> files();
  	BP_Vec<std::string> dirs();

  	bool exists( std::string s);
  	bool is_file( std::string s);
  	bool is_dir( std::string s);

  	void expand_up();

  	void print( std::ostream &sout);
  	void print( std::ostream &sout, BP_GVec_Member<BP_Dir> *start);
  	void print_files( std::ostream &sout);
  	void print_files( std::ostream &sout, BP_GVec_Member<BP_Dir> *start);

  	private:
  	void generate(BP_GVec_Member<BP_Dir> *start);
  	void print_dirs( std::ostream &sout, BP_GVec_Member<BP_Dir> *dir);
  	void print_dirs_files( std::ostream &sout, BP_GVec_Member<BP_Dir> *dir);

  };

  BP_Dir_Tree::BP_Dir_Tree()
  {
  	char path[1024];

  	dir_list.add();
  	dir_list.last().name = ".";
  	//dir_list.last().rel_name = ".";
  	dir_list.last().abs_name = std::string( getcwd( path, sizeof(path)) );
  	//std::cout << "'" << dir_list.last().abs_name << "'" << std::endl;
  	ref = dir_list.last_member();
  	curr = ref;
  	add_vertex( ref);
  	generate( ref);

  };

  void BP_Dir_Tree::print( std::ostream &sout)
  {
  	print_dirs( sout, ref);
  };

  void BP_Dir_Tree::print( std::ostream &sout,  BP_GVec_Member<BP_Dir> *start)
  {
  	print_dirs( sout, start);
  };

  void BP_Dir_Tree::print_files( std::ostream &sout)
  {
  	print_dirs_files( sout, ref);
  };

  void BP_Dir_Tree::print_files( std::ostream &sout,  BP_GVec_Member<BP_Dir> *start)
  {
  	print_dirs_files( sout, start);
  };

  BP_Vec<std::string> BP_Dir_Tree::files()
  {
  	return curr->get_val().files;
  };

  BP_Vec<std::string> BP_Dir_Tree::dirs()
  {
  	return curr->get_val().dirs;
  };

  void BP_Dir_Tree::expand_up()
  {
  	BP_Dir new_dir;
  	BP_GVec_Member<BP_Dir> *new_ref;
  	unsigned long int i;

  	BP_Vec<std::string> dirs = tokenize( ref->get_val().abs_name, "/");
  	ref->get_val().name = dirs.last();

  	new_dir.name = ".";
  	new_dir.abs_name = "";

  	for( i=0; i<dirs.size()-1; i++)
  	{
  		new_dir.abs_name += "/" + dirs[i];
  	}

  	new_ref = dir_list.add(new_dir);
  	add_vertex( new_ref);
  	connect_vertices( new_ref, ref, edge_list.add(), 1);


  	/// add new directories
  	DIR * new_ref_dir;
  	struct dirent *obj;
  	struct stat st;
  	char path[1024];
  	std::string name;


  	// open the new_ref directory
  	new_ref_dir = opendir(new_ref->get_val().abs_name.c_str());
  	if (!new_ref_dir)
  	{
  		std::cout << "BP_Dir_Tree() Error, opendir() failed" << std::endl;
  		exit(1);
  	}
  	errno=0;


  	while( obj = readdir(new_ref_dir) )
  	{
  		//std::cout << obj->d_name << " ";
  		lstat( (new_ref->get_val().abs_name + "/" + std::string(obj->d_name)).c_str(), &st);
  		if(S_ISDIR(st.st_mode))
  		{
  			//std::cout << " DIRECTORY" << std::endl;
  			name = std::string( obj->d_name);
  			if( name != "." && name != "..")
  			{
  				//std::cout << "add dir: " << name << std::endl;
  				new_ref->get_val().dirs.add( name);

  				//BP_pause();

  				// don't add ref, since it already exists
  				if( name != ref->get_val().name)
  				{
  					dir_list.add();
  					dir_list.last().name = name;
  					//dir_list.last().rel_name = start->get_val().rel_name + "/" + name;
  					dir_list.last().abs_name = new_ref->get_val().abs_name + "/" + name;

  					add_vertex( dir_list.last_member());
  					connect_vertices( new_ref, dir_list.last_member(), edge_list.add(), 1);

  					generate(dir_list.last_member());
  				}

  			}
  		}
  		else
  		{
  			//std::cout << " FILE" << std::endl;
  			name = std::string( obj->d_name);
  			new_ref->get_val().files.add( name);

  		}
  	}

  	if (errno){
  		std::cout << "BP_Dir_Tree() Error, readdir() failed" << std::endl;
  		exit(1);
  	}

  	closedir(new_ref_dir);




  	ref = new_ref;
  	curr = new_ref;

  	return;
  };

  void BP_Dir_Tree::cd(std::string s)
  {
  	if( s == "..")
  	{
  		if( num_incident_edges(curr, 1) == 1)
  			curr = get_neighbor(curr, 0, 1);
  		else
  			expand_up();
  	}
  	else
  	{
  		unsigned long int index;
  		if( curr->get_val().dirs.find_first( s, index))
  		{

  			for( unsigned long int i=0; i<num_incident_edges(curr, -1); i++)
  				if( get_neighbor(curr, i, -1)->get_val().name == s)
  				{
  					curr = get_neighbor(curr, i, -1);
  					break;
  				}
  		}
  	}
  };

  void BP_Dir_Tree::cd(BP_GVec_Member<BP_Dir> *dir)
  {
  	curr = dir;
  };

  void BP_Dir_Tree::mkdir(std::string s)
  {
  	mkdir( (curr->get_val().abs_name + "/" + s).c_str());
  };

  std::string BP_Dir_Tree::cwd()
  {
  	return curr->get_val().abs_name;
  };

  BP_GVec_Member<BP_Dir>* BP_Dir_Tree::cwdp()
  {
  	return curr;
  };

  bool BP_Dir_Tree::is_file( std::string s)
  {
  	unsigned long int index;
  	return curr->get_val().files.find_first(s, index);
  };

  bool BP_Dir_Tree::is_dir( std::string s)
  {
  	unsigned long int index;
  	return curr->get_val().files.find_first(s, index);
  };

  bool BP_Dir_Tree::exists( std::string s)
  {
  	if( is_file(s) || is_dir(s))
  		return true;

  	return false;
  };


  void BP_Dir_Tree::generate(BP_GVec_Member<BP_Dir> *start)
  {
  	//std::cout << "begin generate: " << start->get_val().name << std::endl;

  	DIR * start_dir;
  	struct dirent *obj;
  	struct stat st;
  	char path[1024];
  	std::string name;


  	// open the start directory
  	start_dir = opendir(start->get_val().abs_name.c_str());
  	if (!start_dir)
  	{
  		std::cout << "BP_Dir_Tree() Error, opendir() failed" << std::endl;
  		exit(1);
  	}
  	errno=0;


  	while( obj = readdir(start_dir) )
  	{
  		//std::cout << obj->d_name << " ";
  		lstat( (start->get_val().abs_name + "/" + std::string(obj->d_name)).c_str(), &st);
  		if(S_ISDIR(st.st_mode))
  		{
  			//std::cout << " DIRECTORY" << std::endl;
  			name = std::string( obj->d_name);
  			if( name != "." && name != "..")
  			{
  				//std::cout << "add dir: " << name << std::endl;
  				start->get_val().dirs.add( name);

  				//BP_pause();

  				dir_list.add();
  				dir_list.last().name = name;
  				//dir_list.last().rel_name = start->get_val().rel_name + "/" + name;
  				dir_list.last().abs_name = start->get_val().abs_name + "/" + name;

  				add_vertex( dir_list.last_member());
  				connect_vertices( start, dir_list.last_member(), edge_list.add(), 1);

  				generate(dir_list.last_member());
  			}
  		}
  		else
  		{
  			//std::cout << " FILE" << std::endl;
  			name = std::string( obj->d_name);
  			start->get_val().files.add( name);

  		}
  	}

  	if (errno){
  		std::cout << "BP_Dir_Tree() Error, readdir() failed" << std::endl;
  		exit(1);
  	}

  	closedir(start_dir);

  };

  void BP_Dir_Tree::print_dirs( std::ostream &sout, BP_GVec_Member<BP_Dir> *dir)
  {
  	std::cout << dir->get_val().abs_name << std::endl;
  	for( unsigned long int i=0; i<num_incident_edges(dir,-1); i++)
  	{
  		print_dirs( sout, get_neighbor( dir, i, -1));
  	}

  };

  void BP_Dir_Tree::print_dirs_files( std::ostream &sout, BP_GVec_Member<BP_Dir> *dir)
  {
  	unsigned long int i;

  	for( i=0; i<dir->get_val().files.size(); i++)
  		std::cout << dir->get_val().abs_name + "/" + dir->get_val().files[i] << std::endl;

  	for( i=0; i<num_incident_edges(dir,-1); i++)
  	{
  		print_dirs_files( sout, get_neighbor( dir, i, -1));
  	}

  };
  */
}

#endif // BP_Dir_CC

