#ifndef BP_Parse_CC
#define BP_Parse_CC

#include "casm/BP_C++/BP_Parse.hh"

namespace BP {

  // general use functions

  void cut_start_whtspace(std::string &s1) {
    std::string s2;
    int beginning = 1;
    int i = 0;
    while(i < s1.size()) {
      if(beginning == 1) {
        if(s1[i] == ' ' || s1[i] == '\t') i++;
        else {
          beginning = 0;
        }
      }
      if(beginning == 0) {
        s2 += s1[i];
        i++;
      }
    }
    s1 = s2;

  }

  void cut_end_whtspace(std::string &s1) {
    int i = s1.size();

    while(i-- >= 0)
      if(s1[i] != ' ' && s1[i] != '\t')
        break;
    s1.resize(++i);

  }

  std::string cut_start_whtspace(const std::string &s1) {
    std::string s2;
    int beginning = 1;
    int i = 0;
    while(i < s1.size()) {
      if(beginning == 1) {
        if(s1[i] == ' ' || s1[i] == '\t') i++;
        else {
          beginning = 0;
        }
      }
      if(beginning == 0) {
        s2 += s1[i];
        i++;
      }
    }

    return s2;
  }

  std::string cut_end_whtspace(const std::string &s1) {
    std::string s2 = s1;
    int i = s2.size();

    while(i-- >= 0)
      if(s2[i] != ' ' && s2[i] != '\t')
        break;
    s2.resize(++i);
    return s2;
  }

  std::string itos(int i1) {
    std::stringstream ss;
    ss << i1;
    return ss.str();
  }

  std::string dtos(double i1) {
    std::stringstream ss;
    ss << i1;
    return ss.str();
  }

  std::string dtos(double i1, int prec) {
    std::stringstream ss;
    ss.precision(prec);
    ss << i1;
    return ss.str();
  }

  std::string dtos_fixed(double i1, int prec) {
    std::stringstream ss;
    ss.flags(std::ios::fixed);
    ss.precision(prec);
    ss << i1;
    return ss.str();
  }

  std::string dtos_scientific(double i1, int prec) {
    std::stringstream ss;
    ss.flags(std::ios::scientific);
    ss.precision(prec);
    ss << i1;
    return ss.str();
  }

  int stoi(const std::string &s1) {
    return atoi(s1.c_str());
  }

  long int stoli(const std::string &s1) {
    return strtol(s1.c_str(), NULL, 10);
  }

  unsigned long int stouli(const std::string &s1) {
    return strtoul(s1.c_str(), NULL, 10);
  }

  double stod(const std::string &s1) {
    return atof(s1.c_str());
  }

  std::complex<double> stocd(const std::string &s1) {
    std::istringstream ss(s1);
    std::complex<double> c;
    ss >> c;
    return c;
  }

  bool stob(const std::string &str) {
    if(str == "N" || str == "n" || str == "F" || str == "f" || str == "0")
      return false;
    else if(str == "No" || str == "NO" || str == "no")
      return false;
    else if(str == "False" || str == "FALSE" || str == "false")
      return false;
    else if(str == "T" || str == "t" || str == "Y" || str == "y" || str == "1")
      return true;
    else if(str == "Yes" || str == "YES" || str == "yes")
      return true;
    else if(str == "True" || str == "TRUE" || str == "true")
      return true;
    else {
      std::cout << "Error reading bool: " << str << std::endl;
      std::cout << "True options are: 'T', 't', 'True', 'TRUE', 'true', 'Y', 'y', 'Yes', 'YES', 'yes', '1'" << std::endl;
      std::cout << "False options are: 'F', 'f', 'False', 'FALSE', 'false', 'N', 'n', 'No', 'NO', 'no', '0'" << std::endl;
      exit(1);
    }
  }

  BP_Vec<int> stoi(const BP_Vec<std::string> &s_list) {
    BP_Vec<int> i_list;
    i_list.capacity(s_list.size());
    for(int i = 0; i < s_list.size(); i++)
      i_list.add(BP::stoi(s_list[i]));

    return i_list;
  }

  BP_Vec<long int> stoli(const BP_Vec<std::string> &s_list) {
    BP_Vec<long int> i_list;
    i_list.capacity(s_list.size());
    for(int i = 0; i < s_list.size(); i++)
      i_list.add(stoli(s_list[i]));

    return i_list;
  }

  BP_Vec<unsigned long int> stouli(const BP_Vec<std::string> &s_list) {
    BP_Vec<unsigned long int> i_list;
    i_list.capacity(s_list.size());
    for(int i = 0; i < s_list.size(); i++)
      i_list.add(stouli(s_list[i]));

    return i_list;
  }

  BP_Vec<double> stod(const BP_Vec<std::string> &s_list) {
    BP_Vec<double> d_list;
    d_list.capacity(s_list.size());
    for(int i = 0; i < s_list.size(); i++)
      d_list.add(BP::stod(s_list[i]));

    return d_list;
  }

  BP_Vec<std::complex<double> > stocd(const BP_Vec<std::string> &s_list) {
    BP_Vec< std::complex<double> > cd_list;
    cd_list.capacity(s_list.size());
    for(int i = 0; i < s_list.size(); i++)
      cd_list.add(stocd(s_list[i]));

    return cd_list;
  }

  BP_Vec<bool> stob(const BP_Vec<std::string> &s_list) {
    BP_Vec<bool> b_list;
    b_list.capacity(s_list.size());
    for(int b = 0; b < s_list.size(); b++)
      b_list.add(stob(s_list[b]));

    return b_list;
  }

  bool contains(const std::string &s1, const std::string &pattern) {
    if(s1.find(pattern) < s1.size())
      return true;
    else
      return false;
  }

  bool contains(const std::string &s1, const BP_Vec<std::string> &pattern_list) {
    for(int i = 0; i < pattern_list.size(); i++)
      if(contains(s1, pattern_list[i]))
        return true;
    return false;
  }

  void cut_comment(std::string &s1, const std::string &com) {
    cut_comment(s1, BP_Vec<std::string>(1, com));
  }

  void cut_comment(std::string &s1, const BP_Vec<std::string> &com_list) {
    //std::cout << "begin cut_comment( std::string &s1, const BP_Vec<std::string> &com_list)" << std::endl;
    BP_Vec<std::string> lines = tokenize(s1, BP_Vec<std::string>(1, "\n"));

    //std::cout << "lines1: " << std::endl;
    //for( int i=0; i<lines.size(); i++)
    //	std::cout << "  i: " << i << "  : '" << lines[i] << "'" << std::endl;

    for(int i = 0; i < lines.size(); i++)
      lines[i] = trim(lines[i], com_list);

    //std::cout << "lines2: " << std::endl;
    //for( int i=0; i<lines.size(); i++)
    //	std::cout << "  i: " << i << "  : '" << lines[i] << "'" << std::endl;


    std::string s2;
    for(int i = 0; i < lines.size(); i++)
      s2 += lines[i];

    //std::cout << "s2: '" << s2 << "'" << std::endl;

    s1 = s2;
  }

  std::string cut_comment(const std::string &s1, const std::string &com) {
    return cut_comment(s1, BP_Vec<std::string>(1, com));
  }

  std::string cut_comment(const std::string &s1, const BP_Vec<std::string> &com_list) {
    BP_Vec<std::string> lines = tokenize(s1, BP_Vec<std::string>(1, "\n"));
    for(int i = 0; i < lines.size(); i++)
      lines[i] = trim(lines[i], com_list);
    std::string s2;
    for(int i = 0; i < lines.size(); i++)
      s2 += lines[i];
    return s2;

  }

  /*
  std::string trim(const std::string &s1, const BP_Vec<std::string> &delim_list)
  {	// remove comments (delim_list) from end of std::string (s1)

  	std::cout << "begin trim(const std::string &s1, const BP_Vec<std::string> &delim_list)" << std::endl;
  	int i;
  	std::string s2;
  	std::string match;
  	i=0;
  	while( i<s1.size())
  	{
  		match = any_match(s1.substr(i,s1.size()-i),delim_list);
  		//std::cout << "i: " << i << "  match: " << match << std::endl;
  		if(match.size() != 0)	// if a delimiter
  		{
  			return s2;
  		}
  		else	// if not a delimiter
  		{
  			s2 += s1[i];
  		}

  		i++;

  	}


  	return s2;

  }
  */

  std::string trim(const std::string &s1, const BP_Vec<std::string> &delim_list) {
    // remove comments (delim_list) from end of std::string (s1)

    //std::cout << "begin trim(const std::string &s1, const BP_Vec<std::string> &delim_list)" << std::endl;
    int i;
    std::string s2;
    std::string match;

    std::string::size_type len, max_delim = 0;
    for(i = 0; i < delim_list.size(); i++)
      if(delim_list[i].size() > max_delim)
        max_delim = delim_list[i].size();

    i = 0;
    while(i < s1.size()) {

      len = std::min(s1.size() - i, max_delim);
      match = any_match(s1.substr(i, len), delim_list);
      //match = any_match(s1.substr(i,s1.size()-i),delim_list);
      //std::cout << "i: " << i << "  match: " << match << std::endl;
      if(match.size() != 0) {	// if a delimiter
        return s2;
      }
      else {	// if not a delimiter
        s2 += s1[i];
      }

      i++;

    }


    return s2;

  }

  void clean(std::string &s1) {
    cut_end_whtspace(s1);
    cut_start_whtspace(s1);
  }

  void clean(std::string &s1, const BP_Vec<std::string> &com_list) {
    cut_comment(s1, com_list);
    cut_end_whtspace(s1);
    cut_start_whtspace(s1);
  }

  void clean(BP_Vec<std::string> &s1) {
    for(int i = 0; i < s1.size(); i++)
      clean(s1[i]);
  }

  void clean(BP_Vec<std::string> &s1, const BP_Vec<std::string> &com_list) {
    for(int i = 0; i < s1.size(); i++)
      clean(s1[i], com_list);
  }

  std::string clean(const std::string &s1, const BP_Vec<std::string> &com_list) {
    return cut_start_whtspace(cut_end_whtspace(cut_comment(s1, com_list)));
  }

  BP_Vec<std::string> tokenize_whtspace(const std::string &s1) {
    //std::cout << "begin tokenize_whtspace(const std::string &s1)" << std::endl;
    BP_Vec<std::string> delim_list;
    delim_list.add(" ");
    delim_list.add("\t");
    delim_list.add("\f");
    delim_list.add("\n");
    delim_list.add("\r");
    return tokenize(s1, delim_list);
  }

  BP_Vec<std::string> tokenize_whtspace(const std::string &s1, const std::string &str) {
    //std::cout << "begin tokenize_whtspace(const std::string &s1, const std::string &str)" << std::endl;
    BP_Vec<std::string> delim_list(1, str);
    delim_list.add(" ");
    delim_list.add("\t");
    delim_list.add("\f");
    delim_list.add("\n");
    delim_list.add("\r");
    return tokenize(s1, delim_list);
  }

  BP_Vec<std::string> tokenize_whtspace(const std::string &s1, const BP_Vec<std::string> &s_list) {
    //std::cout << "begin tokenize_whtspace(const std::string &s1, const BP_Vec<std::string> &s_list)" << std::endl;
    BP_Vec<std::string> delim_list(s_list);
    delim_list.add(" ");
    delim_list.add("\t");
    delim_list.add("\f");
    delim_list.add("\n");
    delim_list.add("\r");
    return tokenize(s1, delim_list);
  }

  BP_Vec<std::string> tokenize(const std::string &s1, const std::string &str) {
    return tokenize(s1, BP_Vec<std::string>(1, str));
  }

  /*BP_Vec<std::string> tokenize(const std::string &s1, const BP_Vec<std::string> &delim_list)
  {
  	std::cout << "begin tokenize()" << std::endl;
  	std::cout << "          s1: " << s1 << std::endl;
  	std::cout << "  delim_list: " << delim_list << std::endl;
  	int i,j;
  	BP_Vec<std::string> tok_list;
  	std::string s2;
  	std::string match;
  	bool token_started = 0;
  	i=0;
  	while( i<s1.size())
  	{
  		match = any_match(s1.substr(i,s1.size()-i),delim_list);
  		std::cout << "i: " << i << "  match: " << match << std::endl;
  		if(match.size() != 0)	// if a delimiter
  		{
  			if( token_started == 0)
  			{
  				// combine delimiters
  			}
  			else
  			{
  				// finish token
  				tok_list.add(s2);
  				token_started = 0;
  			}
  		}
  		else	// if not a delimiter
  		{
  			if( token_started == 0)
  			{
  				// start token
  				s2.clear();
  				s2 += s1[i];
  				token_started = 1;
  			}
  			else
  			{
  				// add to token
  				s2 += s1[i];
  			}
  		}
  		if( match.size() == 0)	i++;
  		else i += match.size();

  		std::cout << "s2: '" << s2 << "'" << std::endl;
  	}

  	if( token_started)
  		tok_list.add(s2);
  	//std::cout << "tok_list: " << std::endl;
  	//	for( i=0; i<tok_list.size(); i++)
  	//		std::cout << "  i: " << i << "  '" << tok_list[i] << "'" << std::endl;
  	std::cout << "tok_list: " << tok_list << std::endl;
  	std::cout << "finish tokenize()" << std::endl;
  	return tok_list;

  }*/

  std::string any_match(const std::string &s1, const BP_Vec<std::string> &str_list) {
    for(int i = 0; i < str_list.size(); i++) {
      if(s1.substr(0, str_list[i].size()) == str_list[i])
        return str_list[i];
    }

    return std::string();

  }

  BP_Vec<std::string> tokenize(const std::string &s1, const BP_Vec<std::string> &delim_list) {
    //std::cout << "begin tokenize()" << std::endl;
    //std::cout << "          s1: " << s1 << std::endl;
    //std::cout << "  delim_list: " << delim_list << std::endl;
    int i, j;
    BP_Vec<std::string> tok_list;
    std::string s2;
    std::string match;

    std::string::size_type len, max_delim = 0;
    for(i = 0; i < delim_list.size(); i++)
      if(delim_list[i].size() > max_delim)
        max_delim = delim_list[i].size();

    bool token_started = 0;
    i = 0;


    while(i < s1.size()) {
      len = std::min(s1.size() - i, max_delim);
      match = any_match(s1.substr(i, len), delim_list);
      //std::cout << "i: " << i << "  match: " << match << std::endl;
      if(match.size() != 0) {	// if a delimiter
        if(token_started == 0) {
          // combine delimiters
        }
        else {
          // finish token
          tok_list.add(s2);
          token_started = 0;
        }
      }
      else {	// if not a delimiter
        if(token_started == 0) {
          // start token
          s2.clear();
          s2 += s1[i];
          token_started = 1;
        }
        else {
          // add to token
          s2 += s1[i];
        }
      }
      if(match.size() == 0)	i++;
      else i += match.size();

      //std::cout << "s2: '" << s2 << "'" << std::endl;
    }

    if(token_started)
      tok_list.add(s2);
    //std::cout << "tok_list: " << std::endl;
    //	for( i=0; i<tok_list.size(); i++)
    //		std::cout << "  i: " << i << "  '" << tok_list[i] << "'" << std::endl;
    //std::cout << "tok_list: " << tok_list << std::endl;
    //std::cout << "finish tokenize()" << std::endl;
    return tok_list;

  }

  std::string next_string(std::string &s1) {
    //return first token from s1 and cut it from s1

    //std::cout << "begin next_string" << std::endl;

    BP_Vec<std::string> delim_list;
    delim_list.add(" ");
    delim_list.add("\t");
    delim_list.add("\f");
    delim_list.add("\n");
    delim_list.add("\r");
    int i, j;
    std::string s2;
    std::string match;
    bool token_started = 0;
    i = 0;
    while(i < s1.size()) {
      match = any_match(s1.substr(i, s1.size() - i), delim_list);
      //std::cout << "i: " << i << "  match: " << match << std::endl;
      if(match.size() != 0) {	// if a delimiter
        if(token_started == 0) {
          // combine delimiters
        }
        else {
          // finish token
          i += match.size();
          s1 = s1.substr(i, s1.size() - i);
          return s2;
        }
      }
      else {	// if not a delimiter
        if(token_started == 0) {
          // start token
          s2.clear();
          s2 += s1[i];
          token_started = 1;
        }
        else {
          // add to token
          s2 += s1[i];
        }
      }
      if(match.size() == 0)	i++;
      else i += match.size();

      //std::cout << "s2: '" << s2 << "'" << std::endl;
    }


    s1 = s1.substr(i, s1.size() - i);
    return s2;
  };

  double next_double(std::string &s1) {
    return BP::stod(next_string(s1));
  };

  int next_int(std::string &s1) {
    return BP::stoi(next_string(s1));
  };

  long int next_lint(std::string &s1) {
    return stoli(next_string(s1));
  };

  unsigned long int next_ulint(std::string &s1) {
    return stouli(next_string(s1));
  };

  bool next_bool(std::string &s1) {
    std::string str = next_string(s1);
    if(str == "N" || str == "n" || str == "F" || str == "f" || str == "0")
      return false;
    else if(str == "No" || str == "NO" || str == "no")
      return false;
    else if(str == "False" || str == "FALSE" || str == "false")
      return false;
    else if(str == "T" || str == "t" || str == "Y" || str == "y" || str == "1")
      return true;
    else if(str == "Yes" || str == "YES" || str == "yes")
      return true;
    else if(str == "True" || str == "TRUE" || str == "true")
      return true;
    else {
      std::cout << "Error reading bool: " << str << std::endl;
      std::cout << "True options are: 'T', 't', 'True', 'TRUE', 'true', 'Y', 'y', 'Yes', 'YES', 'yes', '1'" << std::endl;
      std::cout << "False options are: 'F', 'f', 'False', 'FALSE', 'false', 'N', 'n', 'No', 'NO', 'no', '0'" << std::endl;
      exit(1);
    }
  };

  std::string &operator>>(std::string &s, int &i) {	///< \ingroup BP_Parse
    i = next_int(s);
    return s;
  };

  std::string &operator>>(std::string &s, long int &i) {	///< \ingroup BP_Parse
    i = next_lint(s);
    return s;
  };

  std::string &operator>>(std::string &s, double &i) {	///< \ingroup BP_Parse
    i = next_double(s);
    return s;
  };

  std::string &operator>>(std::string &s, std::string &i) {	///< \ingroup BP_Parse
    i = next_string(s);
    return s;
  };

  unsigned long int nlines(std::string filename) {
    // count number of lines in a file
    // \f, \n, \r are all counted

    std::ifstream infile(filename.c_str());
    if(!infile) {
      std::cout << "BP_Parse::nlines() error: cannot open '" << filename << "'" << std::endl;
      exit(1);
    }

    unsigned long int N = 0;

    char c;
    infile.get(c);
    while(!infile.eof()) {
      if(c == '\f' || c == '\n' || c == '\r')
        N += 1;
      infile.get(c);

    }

    infile.close();

    return N;
  };

  void read_file_into_array(std::string filename, const BP_Vec<std::string> &com_list, BP_Vec< BP_Vec< int > > &array) {
    // read a file containing an array of whitespace separated values into a double array
    //
    // rows by have varying numbers of fields
    // ignore comments indicated by com_list

    array.erase();
    array.capacity(nlines(filename));

    BP_Parse infile(filename);
    infile.set_com(com_list);

    BP_Vec<int> i_list;

    while(!infile.eof()) {
      i_list = infile.getline_int();
      if(i_list.size() != 0)
        array.add(i_list);
    }

    //array.capacity( array.size());
  };

  void read_file_into_array(std::string filename, const BP_Vec<std::string> &com_list, BP_Vec< BP_Vec< double > > &array) {
    // read a file containing an array of whitespace separated values into a double array
    //
    // rows by have varying numbers of fields
    // ignore comments indicated by com_list

    array.erase();
    array.capacity(nlines(filename));

    BP_Parse infile(filename);
    infile.set_com(com_list);

    BP_Vec<double> d_list;

    while(!infile.eof()) {
      d_list = infile.getline_double();
      if(d_list.size() != 0)
        array.add(d_list);
    }

    //array.capacity( array.size());
  };

  void read_file_into_array(std::string filename, const BP_Vec<std::string> &com_list, BP_Vec< BP_Vec< std::string > > &array) {
    // read a file containing an array of whitespace separated values into a double array
    //
    // rows by have varying numbers of fields
    // ignore comments indicated by com_list

    array.erase();
    array.capacity(nlines(filename));

    BP_Parse infile(filename);
    infile.set_com(com_list);

    BP_Vec<std::string> s_list;

    while(!infile.eof()) {
      s_list = infile.getline_string();
      if(s_list.size() != 0)
        array.add(s_list);
    }

    //array.capacity( array.size());
  };


  /// ///////////////////////////////////////////////
  /// ///////////////////////////////////////////////
  /// class BP_Parse functions (read text files)

  BP_Parse::BP_Parse() {

  };

  BP_Parse::BP_Parse(std::string s1) {
    //std::cout << "s1: " << s1 << std::endl;
    infile.open(s1.c_str());
    if(!infile) {
      std::cout << "BP_Parse error: cannot open '" << s1 << "'" << std::endl;
      exit(1);
    }

  };

  BP_Parse::~BP_Parse() {
    if(infile.is_open()) infile.close();
  };

  bool BP_Parse::try_open(std::string s1) {
    infile.close();
    infile.open(s1.c_str());
    if(!infile) {
      return false;
    }

    return true;
  }

  void BP_Parse::close() {
    infile.close();
  };

  void BP_Parse::reset(std::string s1) {
    infile.close();
    infile.open(s1.c_str());
    if(!infile) {
      std::cout << "BP_Parse error: cannot open '" << s1 << "'" << std::endl;
      exit(1);
    }
  };

  void BP_Parse::add_com(std::string s1) {
    com_list.add(s1);
  };

  void BP_Parse::set_com(const BP_Vec<std::string> &s_list) {
    com_list = s_list;
  };

  void BP_Parse::clear_com() {
    com_list.clear();
  };

  void BP_Parse::remove_com(int i) {
    com_list.remove(i);
  };

  BP_Vec<std::string> BP_Parse::get_com() {
    return com_list;
  };

  std::string BP_Parse::com(int i) {
    return com_list[i];
  };

  std::string BP_Parse::getline() {
    // returns std::string that includes entire line before any comments,
    //   does not include '\f', '\n', '\r'
    //std::cout << "begin getline()" << std::endl;

    std::string str = "";

    char c;
    infile.get(c);
    while(c != '\f' && c != '\n' && c != '\r' && !infile.eof()) {
      str += c;
      infile.get(c);
      //std::cout << "str: " << str << " c: " << c << std::endl;

    }

    //std::getline(infile,str);
    //std::cout << "gl: '" << str << "'" << std::endl;
    cut_comment(str, com_list);
    //std::cout << "gl,cc: '" << str << "'" << std::endl;
    return str;

  };

  std::string BP_Parse::getline_all() {
    // returns std::string that includes entire line including comments
    //  does not include '\f', '\n', '\r'

    std::string str = "";

    char c;
    infile.get(c);
    while(c != '\f' && c != '\n' && c != '\r' && !infile.eof()) {
      str += c;
      infile.get(c);
      //std::cout << "str: " << str << " c: " << c << std::endl;

    }
    //std::cout << "gl: '" << str << "'" << std::endl;
    //std::cout << "gl,cc: '" << str << "'" << std::endl;
    return str;

  };

  std::string BP_Parse::getline_exact() {
    // returns std::string that includes entire line including comments
    //   also includes '\f', '\n', '\r'

    std::string str = "";

    char c;
    do {
      infile.get(c);
      if(!infile.eof())
        str += c;
      //std::cout << "str: " << str << " c: " << c << std::endl;

    }
    while(c != '\f' && c != '\n' && c != '\r' && !infile.eof());

    //std::cout << "gl: '" << str << "'" << std::endl;
    //std::cout << "gl,cc: '" << str << "'" << std::endl;
    return str;

  };

  BP_Vec<std::string> BP_Parse::getline_string() {
    //std::cout << "begin getline_string()" << std::endl;
    return tokenize_whtspace(getline());

  };

  BP_Vec<int> BP_Parse::getline_int() {
    return BP::stoi(tokenize_whtspace(getline()));

  };

  BP_Vec<double> BP_Parse::getline_double() {
    //std::string str;
    //BP_Vec<std::string> s_list;
    //BP_Vec<double> d_list;

    //str = getline();
    //std::cout << "str1: " << str << std::endl;
    //s_list = tokenize_whtspace(str);
    //std::cout << "s_list: " << str << std::endl;
    //for( int i=0; i<s_list.size(); i++)
    //	std::cout << "  i: " << i << "  " << s_list[i] << std::endl;

    //d_list = stod(s_list);
    //std::cout << "d_list: " << str << std::endl;
    //for( int i=0; i<d_list.size(); i++)
    //	std::cout << "  i: " << i << "  " << d_list[i] << std::endl;

    //return d_list;

    return stod(tokenize_whtspace(getline()));

  };

  BP_Vec< std::complex<double> > BP_Parse::getline_cdouble() {
    //std::string str;
    //BP_Vec<std::string> s_list;
    //BP_Vec<double> d_list;

    //str = getline();
    //std::cout << "str1: " << str << std::endl;
    //s_list = tokenize_whtspace(str);
    //std::cout << "s_list: " << str << std::endl;
    //for( int i=0; i<s_list.size(); i++)
    //	std::cout << "  i: " << i << "  " << s_list[i] << std::endl;

    //d_list = stod(s_list);
    //std::cout << "d_list: " << str << std::endl;
    //for( int i=0; i<d_list.size(); i++)
    //	std::cout << "  i: " << i << "  " << d_list[i] << std::endl;

    //return d_list;
    std::cout << "begin getline_cdouble()" << std::endl;
    return stocd(tokenize_whtspace(getline()));

  };

  std::string BP_Parse::next_string() {
    BP_Vec<std::string> delim_list;
    delim_list.add(" ");
    delim_list.add("\t");
    delim_list.add("\f");
    delim_list.add("\n");
    delim_list.add("\r");

    std::string str, match, s;
    char c;

    bool token_started = 0;
    do {
      infile.get(c);
      s = c;
      match = any_match(s, delim_list);
      if(match.size() == 0) {
        token_started = 1;
        str += s;
      }

    }
    while(match.size() == 0 || !token_started);

    return str;

  };

  int BP_Parse::next_int() {
    /*BP_Vec<std::string> delim_list;
    delim_list.add(" ");
    delim_list.add("\t");
    delim_list.add("\f");
    delim_list.add("\n");
    delim_list.add("\r");

    std::string str, match,s;
    char c;

    bool token_started = 0;
    do
    {
    	infile.get(c);
    	s = c;
    	match = any_match(s,delim_list);
    	if( match.size() == 0)
    	{
    		token_started = 1;
    		str += s;
    	}

    }while( match.size() == 0 || !token_started);

    return BP::stoi(str);
    */
    return BP::stoi(next_string());
  };

  long int BP_Parse::next_lint() {
    /*BP_Vec<std::string> delim_list;
    delim_list.add(" ");
    delim_list.add("\t");
    delim_list.add("\f");
    delim_list.add("\n");
    delim_list.add("\r");

    std::string str, match,s;
    char c;

    bool token_started = 0;
    do
    {
    	infile.get(c);
    	s = c;
    	match = any_match(s,delim_list);
    	if( match.size() == 0)
    	{
    		token_started = 1;
    		str += s;
    	}

    }while( match.size() == 0 || !token_started);

    return stoli(str);
    */
    return stoli(next_string());

  };

  unsigned long int BP_Parse::next_ulint() {
    /*BP_Vec<std::string> delim_list;
    delim_list.add(" ");
    delim_list.add("\t");
    delim_list.add("\f");
    delim_list.add("\n");
    delim_list.add("\r");

    std::string str, match,s;
    char c;

    bool token_started = 0;
    do
    {
    	infile.get(c);
    	s = c;
    	match = any_match(s,delim_list);
    	if( match.size() == 0)
    	{
    		token_started = 1;
    		str += s;
    	}


    }while( match.size() == 0 || !token_started);

    return stouli(str);
    */
    return stouli(next_string());
  };

  double BP_Parse::next_double() {
    /*
    BP_Vec<std::string> delim_list;
    delim_list.add(" ");
    delim_list.add("\t");
    delim_list.add("\f");
    delim_list.add("\n");
    delim_list.add("\r");

    std::string str, match,s;
    char c;

    bool token_started = 0;
    do
    {
    	infile.get(c);
    	s = c;
    	match = any_match(s,delim_list);
    	if( match.size() == 0)
    	{
    		token_started = 1;
    		str += s;
    	}

    }while( match.size() == 0 || !token_started);

    return stod(str);
    */
    return BP::stod(next_string());
  };

  std::complex<double> BP_Parse::next_cdouble() {
    return stocd(next_string());
  };

  bool BP_Parse::eof() {
    return infile.eof();
  };

  std::ifstream &BP_Parse::get_istream() {
    return infile;
  };


  /// ///////////////////////////////////////////////
  /// ///////////////////////////////////////////////
  /// class BP_bParse functions (read binary files)

  BP_bParse::BP_bParse(std::string s1) {
    //std::cout << "s1: " << s1 << std::endl;
    infile.open(s1.c_str());
    if(!infile) {
      std::cout << "BP_Parse error: cannot open '" << s1 << "'" << std::endl;
      exit(1);
    }

  };

  BP_bParse::~BP_bParse() {
    if(infile.is_open()) infile.close();
  };

  void BP_bParse::close() {
    infile.close();
  };

  void BP_bParse::reset(std::string s1) {
    infile.close();
    infile.open(s1.c_str());
    if(!infile) {
      std::cout << "BP_Parse error: cannot open '" << s1 << "'" << std::endl;
      exit(1);
    }
  };

  char BP_bParse::next_char() {
    char c;
    infile.read((char *) &c, sizeof(char));
    return c;

  };

  int BP_bParse::next_int() {
    int i;
    infile.read((char *) &i, sizeof(int));
    return i;

  };

  unsigned long int BP_bParse::next_ulint() {
    unsigned long int i;
    infile.read((char *) &i, sizeof(unsigned long int));
    return i;
  };

  long int BP_bParse::next_lint() {
    long int i;
    infile.read((char *) &i, sizeof(long int));
    return i;
  };

  double BP_bParse::next_double() {
    double d;
    infile.read((char *) &d, sizeof(double));
    return d;
  };

  void BP_bParse::peek() {
    infile.peek();
  };

  bool BP_bParse::eof() {
    return infile.eof();
  };


  /// ///////////////////////////////////////////////
  /// ///////////////////////////////////////////////
  /// class BP_Write (write text files)

  //public:
  BP_Write::BP_Write() {

  };

  BP_Write::BP_Write(std::string s) {
    if(outfile.is_open()) outfile.close();

    filename = s;
    append();

  };


  BP_Write::~BP_Write() {
    if(outfile.is_open()) outfile.close();
  };

  void BP_Write::changefile(std::string s) {
    if(outfile.is_open()) outfile.close();

    filename = s;
    append();

  };

  void BP_Write::changefile(std::string s, std::string mode) {
    if(outfile.is_open()) outfile.close();
    filename = s;
    append();

  };

  void BP_Write::newfile() {
    if(outfile.is_open()) outfile.close();

    outfile.open(filename.c_str());
  };

  void BP_Write::newfile(std::string s) {
    if(outfile.is_open()) outfile.close();

    filename = s;
    outfile.open(filename.c_str());
  };

  std::string BP_Write::name() {
    return filename;
  };

  void BP_Write::close() {
    if(outfile.is_open()) outfile.close();
  };

  std::ofstream &BP_Write::get_ostream() {
    return outfile;
  };

  // private:
  void BP_Write::append() {
    if(!outfile.is_open()) {
      outfile.open(filename.c_str(), std::ios::app);
    }
  };


  /// ///////////////////////////////////////////////
  /// ///////////////////////////////////////////////
  /// class BP_bWrite (write binary files)

  //public:
  BP_bWrite::BP_bWrite() {

  };

  BP_bWrite::BP_bWrite(std::string s) {
    if(outfile.is_open()) outfile.close();

    filename = s;
    append();

  };

  BP_bWrite::~BP_bWrite() {
    if(outfile.is_open()) outfile.close();
  };

  void BP_bWrite::changefile(std::string s) {
    if(outfile.is_open()) outfile.close();

    filename = s;
    append();

  };

  void BP_bWrite::changefile(std::string s, std::string mode) {
    if(outfile.is_open()) outfile.close();

    filename = s;
    append();

  };

  void BP_bWrite::newfile() {
    if(outfile.is_open()) outfile.close();

    outfile.open(filename.c_str());
  };

  std::string BP_bWrite::name() {
    return filename;
  };

  void BP_bWrite::close() {
    if(outfile.is_open()) outfile.close();
  };

  void BP_bWrite::write(const char *s, unsigned long int n) {
    append();
    outfile.write(s, n);
  };

  void BP_bWrite::write_char(char c) {
    append();
    outfile.write(&c, sizeof(char));
  };

  void BP_bWrite::write_int(int i) {
    append();
    outfile.write((char *) &i, sizeof(int));
  };

  void BP_bWrite::write_lint(long int i) {
    append();
    outfile.write((char *) &i, sizeof(long int));
  };

  void BP_bWrite::write_ulint(unsigned long int i) {
    append();
    outfile.write((char *) &i, sizeof(unsigned long int));
  };

  void BP_bWrite::write_double(double d) {
    append();
    outfile.write((char *) &d, sizeof(double));
  };

  //private:
  void BP_bWrite::append() {
    if(!outfile.is_open()) {
      outfile.open(filename.c_str(), std::ios::app | std::ios::binary);
    }
  };



}

#endif // BP_Parse_CC

