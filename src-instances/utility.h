#include"inclusions.h"

#include <time.h>

using namespace std;

#ifndef __UTILITY__
#define __UTILITY__

namespace utility {

//! increments a counter in base-k by 1
void increment(vector <double> & x, int base, int position);

double sum(vector<double> & x);

double average(vector<double> & x);
double variance(vector<double> & x);
double dev_std(vector<double> & x);
double min(vector<double>& vect);
double max(vector<double>& vect);
 
//! checks the eigenvalues-eigenvectors relation
void eigenval_check(LaGenMatDouble & A, double lambda, LaVectorDouble & x, double tol);

//! return the machine epsilon
double machine_epsilon();

//! 1-norm normalization, i.e. the normalized vector x sums up to 1:
void normalize_1norm(vector<double> & x);
void normalize_1norm(LaVectorDouble & x);

//! 2-norm normalization
void normalize_2norm(LaVectorDouble & x, double & gamma);
void normalize_2norm(LaVectorDouble & x);

//! simple greather than ordering class, for double vectors
class greater_than_ordering {
public:
	greater_than_ordering() {}
	bool operator()(const int& a, const int& b)
	{
		return (values[a] < values[b]);
	}

	void set_vectors(const vector<double>& values)
	{
		this->values = values;
	}

private:
	vector<double> values;
};

//! sorts a vector according to a greather than ordering; returns the ordered index vector
vector<int> sort_gt(vector<double>& vect);


} // end of namespace utility

class Hystorical
{
public:
	vector<double> current;
	vector<double> best;
	
	void print(string file_name);
};


//! prints a vector of class T elements: note: it works even for matrices! it is automatically recursively called!!
template <class T>
ostream& operator<< (ostream& os, vector<T>& vect)
{
// 	typename std::vector<T>::const_iterator it;

	os << "[";
// 	for (it = vect.begin(); it != vect.end(); it++)
	for (int i = 0; i < int(vect.size()); i++) {
		os << vect[i];
		if (i < int(vect.size())-1) cout << " ";
	}
	os << "]";
	return os;
}

//! prints a LaVectorDouble
ostream& operator<< (ostream& os, const LaVectorDouble & vec);

//! prints a LaGenMatDouble
ostream& operator<< (ostream& os, const LaGenMatDouble & M);

//empty ostream
struct nullstream: 
  std::ostream { 
 nullstream(): std::ios(0), std::ostream(0) {} 
}; 

#include <stdio.h>
#include <time.h>


#endif
