#include "ampl_output_reader.h"


bool
ampl_output_reader::read_value(const string & name, double & value)
{
	if (!seek_by_string(name))
		return false;

	seek_next();	// skip '='
	value = read_token_dbl();
	return true;
}


bool
ampl_output_reader::read_vector(const string & name, vector<double> & v, int size)
{
	if (!seek_by_string(name))
		return false;
	
	seek_next();	// skip '[*]'
	seek_next();	// skip ':='

	int v_index = 0;
	while (true)
	{
		string index(read_token());
		if (index == string(";")) break;	// end of vector
		v_index = atoi(index.c_str());
		v[v_index-1] = read_token_dbl();
	}
	return true;
}



bool
ampl_output_reader::read_matrix(const string & name, vector<LaVectorDouble> & m, int rows, int columns)
{
	if (!seek_by_string(name))
		return false;
	
	string head_1(read_token());
	if (head_1 == string("[*,*]"))
	{
		seek_next();	// skip ':'
		for (int i = 0; i < columns; i++)
			seek_next();	// skip column index
		
		seek_next();	// skip ':='
		for (int i = 0; i < rows; i++)
		{
			seek_next();	// skip row index
			for (int ii = 0; ii < columns; ii++)
				m[i](ii) = read_token_dbl();
		}
	}
	else
	{
		for (int i = 0; i < rows; i++)
		{
			for (int ii = 0; ii < columns; ii++)
			{
				seek_next();	// skip row index
				seek_next();	// skip column index
				m[i](ii) = read_token_dbl();
			}
		}
	}

	return true;
}

bool
ampl_output_reader::seek_by_string(const string & str_to_find, int i_th)
{
	seek_begin();
	string str_token;
	
	while (!kfile.eof())
	{
		kfile >> str_token;
		if (!str_token.compare(str_to_find))
		{
			if (i_th == 1)
				return true;
			else
				i_th--;
		}
	}
	return false;
}
