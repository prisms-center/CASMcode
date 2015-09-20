#ifndef BP_zParse_HH
#define BP_zParse_HH

#include <iostream>
#include <istream>
#include <fstream>
#include <string>
#include "casm/BP_C++/BP_Vec.hh"
#include "casm/BP_C++/BP_Parse.hh"
#include "casm/external/gzstream/gzstream.h"

namespace BP {


  class BP_zParse {
  private:
    igzstream infile;
    BP_Vec<std::string> com_list;

  public:

    BP_zParse(std::string s1) {
      //std::cout << "s1: " << s1 << std::endl;
      infile.open(s1.c_str());
      if(!infile) {
        std::cout << "BP_zParse error: cannot open " << s1 << std::endl;
        exit(1);
      }

    };

    ~BP_zParse() {
      //if( infile.is_open())
      infile.close();
    };

    void close() {
      infile.close();
    };

    void reset(std::string s1) {
      infile.close();
      infile.open(s1.c_str());
      if(!infile) {
        std::cout << "BP_Parse error: cannot open " << s1 << std::endl;
        exit(1);
      }
    };

    void add_com(std::string s1) {
      com_list.add(s1);
    };

    void clear_com() {
      com_list.clear();
    };

    void remove_com(int i) {
      com_list.remove(i);
    };

    BP_Vec<std::string> get_com() {
      return com_list;
    };

    std::string com(int i) {
      return com_list[i];
    };

    std::string getline() {
      std::string str;
      std::getline(infile, str);
      //std::cout << "gl: '" << str << "'" << std::endl;
      cut_comment(str, com_list);
      //std::cout << "gl,cc: '" << str << "'" << std::endl;
      return str;

    };

    std::string getline_all() {
      std::string str;
      std::getline(infile, str);
      //std::cout << "gl: '" << str << "'" << std::endl;
      //std::cout << "gl,cc: '" << str << "'" << std::endl;
      return str;

    };

    BP_Vec<std::string> getline_string() {
      return tokenize_whtspace(getline());

    };

    BP_Vec<int> getline_int() {
      return stoi(tokenize_whtspace(getline()));

    };

    BP_Vec<double> getline_double() {
      std::string str;
      BP_Vec<std::string> s_list;
      BP_Vec<double> d_list;

      str = getline();
      //std::cout << "str1: " << str << std::endl;
      s_list = tokenize_whtspace(str);
      //std::cout << "s_list: " << str << std::endl;
      //for( int i=0; i<s_list.size(); i++)
      //	std::cout << "  i: " << i << "  " << s_list[i] << std::endl;

      d_list = stod(s_list);
      //std::cout << "d_list: " << str << std::endl;
      //for( int i=0; i<d_list.size(); i++)
      //	std::cout << "  i: " << i << "  " << d_list[i] << std::endl;

      return d_list;

      //return stod( tokenize_whtspace( getline()));

    };

    std::string next_string() {
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

    int next_int() {
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

      return BP::stoi(str);

    };

    long int next_lint() {
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

      return stoli(str);

    };

    unsigned long int next_ulint() {
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

      return stouli(str);

    };

    double next_double() {
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

      return BP::stod(str);

    };

    bool eof() {
      return infile.eof();
    };

    template <class T> igzstream &operator >>(T &t) {
      infile >> t;
      return infile;
    };

  };


  class BP_bzParse {
  private:
    igzstream infile;

  public:

    BP_bzParse(std::string s1) {
      //std::cout << "s1: " << s1 << std::endl;
      infile.open(s1.c_str());
      if(!infile) {
        std::cout << "BP_Parse error: cannot open " << s1 << std::endl;
        exit(1);
      }

    };

    ~BP_bzParse() {
      //if( infile.is_open())
      infile.close();
    };

    void close() {
      infile.close();
    };

    void reset(std::string s1) {
      infile.close();
      infile.open(s1.c_str());
      if(!infile) {
        std::cout << "BP_Parse error: cannot open " << s1 << std::endl;
        exit(1);
      }
    };

    char next_char() {
      char c;
      infile.read((char *) &c, sizeof(char));
      return c;

    };

    int next_int() {
      int i;
      infile.read((char *) &i, sizeof(int));
      return i;

    };

    unsigned long int next_ulint() {
      unsigned long int i;
      infile.read((char *) &i, sizeof(unsigned long int));
      return i;
    };

    long int next_lint() {
      long int i;
      infile.read((char *) &i, sizeof(long int));
      return i;
    };

    double next_double() {
      double d;
      infile.read((char *) &d, sizeof(double));
      return d;
    };

    void peek() {
      infile.peek();
    };

    bool eof() {
      return infile.eof();
    };

  };



  class BP_zWrite {
  private:
    ogzstream outfile;
    std::string filename;

  public:

    BP_zWrite() {

    };

    BP_zWrite(std::string s) {
      //if( outfile.is_open())
      outfile.close();

      filename = s;

    };


    ~BP_zWrite() {
      //if( outfile.is_open())
      outfile.close();
    };

    void changefile(std::string s) {
      //if( outfile.is_open())
      outfile.close();

      filename = s;

    };

    void changefile(std::string s, std::string mode) {
      //if( outfile.is_open())
      outfile.close();

      filename = s;

    };

    void newfile() {
      //if( outfile.is_open())
      outfile.close();

      outfile.open(filename.c_str());
    };

    std::string name() {
      return filename;
    };

    void close() {
      //if( outfile.is_open())
      outfile.close();
    };

    template <class T> ogzstream &operator <<(const T &t) {
      outfile << t;
      return outfile;
    };



  };

  class BP_bzWrite {
  private:
    ogzstream outfile;
    std::string filename;

  public:

    BP_bzWrite() {

    };

    BP_bzWrite(std::string s) {
      //if( outfile.is_open())
      outfile.close();

      filename = s;

    };

    ~BP_bzWrite() {
      //if( outfile.is_open())
      outfile.close();
    };

    void changefile(std::string s) {
      //if( outfile.is_open())
      outfile.close();

      filename = s;

    };

    void changefile(std::string s, std::string mode) {
      //if( outfile.is_open())
      outfile.close();

      filename = s;

    };

    void newfile() {
      //if( outfile.is_open())
      outfile.close();

      outfile.open(filename.c_str());
    };

    std::string name() {
      return filename;
    };

    void close() {
      //if( outfile.is_open())
      outfile.close();
    };



    void write(const char *s, unsigned long int n) {
      outfile.write(s, n);
    };

    void write_char(char c) {
      outfile.write(&c, sizeof(char));
    };

    void write_int(int i) {
      outfile.write((char *) &i, sizeof(int));
    };

    void write_lint(long int i) {
      outfile.write((char *) &i, sizeof(long int));
    };

    void write_ulint(unsigned long int i) {
      outfile.write((char *) &i, sizeof(unsigned long int));
    };

    void write_double(double d) {
      outfile.write((char *) &d, sizeof(double));
    };



  };

}

#endif // BP_zParse_HH

