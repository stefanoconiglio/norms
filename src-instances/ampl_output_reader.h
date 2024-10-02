#include "inclusions.h"

using namespace std;

#ifndef __AMPL_OUTPUT_READER__
#define __AMPL_OUTPUT_READER__

class ampl_output_reader
{

protected:
	fstream kfile;
	string path;


public:
	ampl_output_reader(string file_path) : kfile(file_path.c_str(), ios_base::out | ios_base::in), path(file_path) {}

public:

	// GENERAL READING METHODS
	inline
	string read_token() {string str_token; kfile >> str_token; return str_token;}
	inline
	double read_token_dbl() {double token; kfile >> token; return token;}

	// SPECIFIC AMPL OUTPUT READING METHODS
	
	bool read_value(const string & name, double & value);
	bool read_vector(const string & name, vector<double>& vect, int size = -1);
	bool read_matrix(const string & name, vector<LaVectorDouble> & matrix, int rows, int columns);

// MOVERS
public:
	inline
	void seek_begin() {kfile.seekg(0);}
	inline
	void seek_end() {kfile.seekg(0, ios::end);}
	inline
	void seek_next() {string str_token; kfile >> str_token;}

	//! Atm I don't know what it does...
	bool seek_by_string(const string & str_to_find, int i_th = 1);
};

#endif
