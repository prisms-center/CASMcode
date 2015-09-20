#ifndef BP_Parse_HH
#define BP_Parse_HH

#include <iostream>
#include <istream>
#include <fstream>
#include <string>
#include <cstdio>
#include <cstdlib>
#include <sstream>
#include <complex>
#include "casm/BP_C++/BP_useful.hh"
#include "casm/BP_C++/BP_Vec.hh"

namespace BP {

  void			cut_start_whtspace(std::string &s1);											///< \ingroup BP_Parse
  void			cut_end_whtspace(std::string &s1);											///< \ingroup BP_Parse
  std::string			cut_start_whtspace(const std::string &s1);									///< \ingroup BP_Parse
  std::string			cut_end_whtspace(const std::string &s1);										///< \ingroup BP_Parse
  std::string			itos(int i1);															///< \ingroup BP_Parse
  std::string			dtos(double i1);														///< \ingroup BP_Parse
  std::string			dtos(double i1, int prec);												///< \ingroup BP_Parse
  std::string			dtos_fixed(double i1, int prec);										///< \ingroup BP_Parse
  std::string			dtos_scientific(double i1, int prec);									///< \ingroup BP_Parse
  int				stoi(const std::string &s1);												///< \ingroup BP_Parse
  long int		stoli(const std::string &s1);												///< \ingroup BP_Parse
  unsigned long int	stouli(const std::string &s1);											///< \ingroup BP_Parse
  double			stod(const std::string &s1);												///< \ingroup BP_Parse
  std::complex<double> stocd(const std::string &s1);												///< \ingroup BP_Parse
  bool			stob(const std::string &str);
  BP_Vec<int>		stoi(const BP_Vec<std::string> &s_list);									///< \ingroup BP_Parse
  BP_Vec<long int>	stoli(const BP_Vec<std::string> &s_list);								///< \ingroup BP_Parse
  BP_Vec<unsigned long int>	stouli(const BP_Vec<std::string> &s_list);						///< \ingroup BP_Parse
  BP_Vec<double>	stod(const BP_Vec<std::string> &s_list);									///< \ingroup BP_Parse
  BP_Vec<std::complex<double> > stocd(const BP_Vec<std::string> &s_list);							///< \ingroup BP_Parse
  BP_Vec<bool>	stob(const BP_Vec<std::string> &s_list);
  std::string			any_match(const std::string &s1, const BP_Vec<std::string> &str_list);			///< \ingroup BP_Parse
  bool			contains(const std::string &s1, const std::string &pattern);						///< \ingroup BP_Parse
  bool			contains(const std::string &s1, const BP_Vec<std::string> &pattern);				///< \ingroup BP_Parse
  void			cut_comment(std::string &s1, const std::string &com);							///< \ingroup BP_Parse
  void			cut_comment(std::string &s1, const BP_Vec<std::string> &com_list);				///< \ingroup BP_Parse
  std::string			cut_comment(const std::string &s1, const std::string &com);						///< \ingroup BP_Parse
  std::string			cut_comment(const std::string &s1, const BP_Vec<std::string> &com_list);			///< \ingroup BP_Parse
  std::string			trim(const std::string &s1, const BP_Vec<std::string> &delim_list);				///< \ingroup BP_Parse
  void			clean(std::string &s1);														///< \ingroup BP_Parse
  void			clean(std::string &s1, const BP_Vec<std::string> &com_list);						///< \ingroup BP_Parse
  void			clean(BP_Vec<std::string> &s1);												///< \ingroup BP_Parse
  void			clean(BP_Vec<std::string> &s1, const BP_Vec<std::string> &com_list);				///< \ingroup BP_Parse
  std::string			clean(const std::string &s1, const BP_Vec<std::string> &com_list);				///< \ingroup BP_Parse
  BP_Vec<std::string>	tokenize_whtspace(const std::string &s1);									///< \ingroup BP_Parse
  BP_Vec<std::string>	tokenize_whtspace(const std::string &s1, const std::string &delim);				///< \ingroup BP_Parse
  BP_Vec<std::string>	tokenize_whtspace(const std::string &s1, const BP_Vec<std::string> &delim_list);	///< \ingroup BP_Parse
  BP_Vec<std::string>	tokenize(const std::string &s1, const std::string &delim);						///< \ingroup BP_Parse
  BP_Vec<std::string>	tokenize(const std::string &s1, const BP_Vec<std::string> &delim_list);			///< \ingroup BP_Parse
  std::string			next_string(std::string &s1);												///< \ingroup BP_Parse
  double			next_double(std::string &s1);												///< \ingroup BP_Parse
  int				next_int(std::string &s1);													///< \ingroup BP_Parse
  long int		next_lint(std::string &s1);													///< \ingroup BP_Parse
  unsigned long int	next_ulint(std::string &s1);											///< \ingroup BP_Parse
  bool			next_bool(std::string &s1);													///< \ingroup BP_Parse
  std::string			&operator>>(std::string &s, int &i);											///< \ingroup BP_Parse
  std::string			&operator>>(std::string &s, long int &i);										///< \ingroup BP_Parse
  std::string			&operator>>(std::string &s, double &i);										///< \ingroup BP_Parse
  std::string			&operator>>(std::string &s, std::string &i);										///< \ingroup BP_Parse

  unsigned long int	nlines(std::string filename);
  void			read_file_into_array(std::string filename, const BP_Vec<std::string> &com_list, BP_Vec< BP_Vec< int > > &array);	///< \ingroup BP_Parse
  void			read_file_into_array(std::string filename, const BP_Vec<std::string> &com_list, BP_Vec< BP_Vec< double > > &array);	///< \ingroup BP_Parse
  void			read_file_into_array(std::string filename, const BP_Vec<std::string> &com_list, BP_Vec< BP_Vec< std::string > > &array);	///< \ingroup BP_Parse

  // write text files
  ///< \ingroup BP BP_Parse
  class BP_Parse {
  private:
    std::ifstream infile;
    BP_Vec<std::string> com_list;

  public:

    // // Initialize a BP_Parse object to read file 'filename':
    // BP_Parse infile( filename);
    //
    // // or:
    // BP_Parse infile;
    // infile.try_open(filename);
    //
    // // Add comments to ignore while reading the file
    // infile.add_com("//");
    // infile.add_com("#");
    //
    // // After doing this, lines will be read up to the comment
    // //   For instance, consider the line:
    // //     "test, this is a test // really, it's a test"
    // std::string str1 = infile.getline();
    // cout << str1 << endl; // would give this: "test, this is a test "
    //
    // // To include the comments, use:
    // std::string str2 = infile.getline_all();
    // cout << str2 << endl; // would give this: "test, this is a test // really, it's a test"
    //
    // // To also include trailing characters '\f', '\n', or '\r', use:
    // std::string str3 = infile.getline_exact();
    //
    // // To parse a line by whitespace, and return a vector of values use:
    // BP_Vec<std::string> = infile.getline_string();	// return a vector of strings
    // BP_Vec<std::string> = infile.getline_int();		// return a vector of ints
    // BP_Vec<std::string> = infile.getline_double();	// return a vector of doubles
    // BP_Vec<std::string> = infile.getline_cdouble();	// return a vector of complex doubles
    //
    // // To parse a line by whitespace and some other value:
    // std::string str = infile.getline();
    // BP_Vec<std::string> = tokenize_whtspace(str, ",");
    // BP_Vec<int> = stoi(tokenize_whtspace(str, ","));
    //
    // // To parse a line by some other value (but not whitespace):
    // std::string str = infile.getline();
    // BP_Vec<std::string> = tokenize(str, ",");
    // BP_Vec<int> = stoi(tokenize(str, ","));
    //
    // // To read the next value from the string (as separated by whitespace):
    // std::string str = infile.next_string();
    // int i = infile.next_int();
    // // etc.
    //
    // // To pass the stream to some function that reads from the file:
    // //   for example: void SomeClass::read(std::istream &instream);
    // something.read( infile.get_istream());
    //
    // // You can also use the stream extraction operator ( >> ),
    // //   but this does check for comments
    // std::string str;
    // infile >> str;


    BP_Parse();
    BP_Parse(std::string s1);
    ~BP_Parse();

    bool try_open(std::string s1);
    void close();
    void reset(std::string s1);
    void add_com(std::string s1);
    void set_com(const BP_Vec<std::string> &s_list);
    void clear_com();
    void remove_com(int i);
    BP_Vec<std::string> get_com();
    std::string com(int i);
    std::string getline();
    std::string getline_all();
    std::string getline_exact();
    BP_Vec<std::string> getline_string();
    BP_Vec<int> getline_int();
    BP_Vec<double> getline_double();
    BP_Vec< std::complex<double> > getline_cdouble();
    std::string next_string();
    int next_int();
    long int next_lint();
    unsigned long int next_ulint();
    double next_double();
    std::complex<double> next_cdouble();
    bool eof();
    std::ifstream &get_istream();
    template <class T> std::ifstream &operator >>(T &t);

  };

  template <class T> std::ifstream &BP_Parse::operator >>(T &t) {
    infile >> t;
    return infile;
  };



  // read binary files
  ///< \ingroup BP BP_Parse
  class BP_bParse {
  private:
    std::ifstream infile;

  public:

    BP_bParse(std::string s1);
    ~BP_bParse();

    void close();
    void reset(std::string s1);
    char next_char();
    int next_int();
    unsigned long int next_ulint();
    long int next_lint();
    double next_double();
    void peek();
    bool eof();

  };


  // write text files
  ///< \ingroup BP BP_Parse
  class BP_Write {
  private:
    std::ofstream outfile;
    std::string filename;

  public:

    // // Initialize a BP_Write object:
    // //   By default, files are opened in "append" mode
    // BP_Write file("outfile.txt");
    //
    // // or:
    // BP_Write file;
    // file.changefile("outfile.txt");
    //
    // // To overwrite a file, use:
    // file.newfile();
    //
    // // or:
    // file.newfile("newfilename.txt");
    //
    // // You can write using the stream insertion operator ( << ):
    // file << "Write this to the file." << endl;
    //
    // // You can't do an 'endl' without other output,
    // // So this is not allowed: file << endl;
    // // But this is allowed: file << " " << endl;
    // // So instead of 'endl' do this:
    // file << '\n' << flush;
    //
    // // To pass the stream to some function that writes to the file:
    // //   for example: void SomeClass::write(std::ostream &outstream);
    // something.write( file.get_ostream());



    BP_Write();
    BP_Write(std::string s);
    ~BP_Write();

    void changefile(std::string s);
    void changefile(std::string s, std::string mode);
    void newfile();
    void newfile(std::string s);
    std::string name();
    void close();
    std::ofstream &get_ostream();
    template <class T> std::ofstream &operator <<(const T &t);

  private:
    void append();

  };

  template <class T> std::ofstream &BP_Write::operator <<(const T &t) {
    append();
    outfile << t;
    return outfile;
  };

  // write binary files
  ///< \ingroup BP BP_Parse
  class BP_bWrite {
  private:
    std::ofstream outfile;
    std::string filename;

  public:

    BP_bWrite();
    BP_bWrite(std::string s);
    ~BP_bWrite();

    void changefile(std::string s);
    void changefile(std::string s, std::string mode);
    void newfile();
    std::string name();
    void close();
    void write(const char *s, unsigned long int n);
    void write_char(char c);
    void write_int(int i);
    void write_lint(long int i);
    void write_ulint(unsigned long int i);
    void write_double(double d);

  private:
    void append();

  };



}

#endif // BP_Parse_HH

