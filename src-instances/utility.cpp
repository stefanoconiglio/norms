#include "utility.h"

using namespace std;

double
utility::sum( vector<double> & x )
{
  double sum = 0;
  for (int i = 0; i < int(x.size()); i++)
    sum += x[i];
  return sum;
}

void
utility::normalize_1norm(vector<double> & x)
{
	double sum = 0;
	for (int i = 0; i < int(x.size()); i++)
		sum += fabs(x[i]);
	
	for (int i = 0; i < int(x.size()); i++)
		x[i] /= sum;
}

void
utility::normalize_1norm(LaVectorDouble & x)
{
	double sum = 0;
	for (int i = 0; i < int(x.size()); i++)
	  sum += fabs(x(i));
	
	for (int i = 0; i < int(x.size()); i++)
	  x(i) /= sum;
}

void
utility::normalize_2norm(LaVectorDouble & x, double & gamma)
{
	double norm = Blas_Norm2(x);
	for (int i = 0; i < int(x.size()); i++)
	{
		x(i) /= norm;
	}
	gamma /= norm;
}

void
utility::normalize_2norm(LaVectorDouble & x)
{
	double norm = Blas_Norm2(x);
	for (int i = 0; i < int(x.size()); i++)
	{
		x(i) /= norm;
	}
}

double
utility::average(vector <double> & x)
{
	double sum = 0;
	vector<double>::iterator it;
	for (it = x.begin(); it != x.end(); it++)
		sum += *it;
	sum /= x.size();
	return sum;
}

double
utility::variance(vector<double> & x)
{
	if (x.size() == 1)
		return 0;

	double sum = 0;
	double avg = average(x);
	vector<double>::iterator it;
	for (it = x.begin(); it != x.end(); it++)
		sum += pow(*it - avg,2);
	sum /= (x.size()-1);
	return sum;
}

double
utility::dev_std(vector<double> & x)
{
  return sqrt(variance(x));
}

double
utility::min(vector<double>& vect) {
  double ret = 1e300;
  for (int i = 0; i < int(vect.size()); i++)
    if ( vect[i] < ret )
      ret = vect[i];
 
  return ret;
}


double
utility::max(vector<double>& vect) {
  double ret = -1e300;
  for (int i = 0; i < int(vect.size()); i++)
    if ( vect[i] > ret )
      ret = vect[i];
 
  return ret;
}


void
utility::increment(vector <double> & x, int base, int position)
{
	if (position == -1)
		return;

	if ( x[position] < base-1)
		x[position]++;
	else
		x[position] = 0;
		increment(x, base, position-1);
}


void
utility::eigenval_check(LaGenMatDouble & A, double lambda, LaVectorDouble & x, double tol)
{
	LaVectorDouble Ax(x.size());
	Ax = A * x;
	LaVectorDouble lambdaX(x.size());
	lambdaX = x;
	lambdaX *= lambda;
	assert(Blas_Norm2(Ax-lambdaX) <= tol);
}

double
utility::machine_epsilon()
{
	double guess = 1e-10;
	double sum;
	double diff;
	while (1)
	{
		sum = 1 + guess;
		diff = sum - 1;
		if (diff == 0)
			break;
		else
			guess /= (double)10;
	}
	return guess;
}

vector<int>
utility::sort_gt(vector<double>& vect)
{
	int n = vect.size();
	vector<int> ordered_vector;
	for (int i = 0; i < n; i++)
	{
		ordered_vector.push_back(i);
	}

	greater_than_ordering gt;
	gt.set_vectors(vect);
	sort(ordered_vector.begin(), ordered_vector.end(), gt);

	return ordered_vector; 
} 

void
Hystorical::print(string file_name) {
	ofstream hystorical(file_name.c_str());
	hystorical << "h_c = [";
	for (int i = 0; i < int(current.size()); i++)
		hystorical << current[i] << " ";
	hystorical << "];" << endl;
	hystorical << "h_b = [";
	for (int i = 0; i < int(current.size()); i++)
		hystorical << current[i] << " ";
	hystorical << "];" << endl;
}

ostream&
operator<< (ostream& os, const LaVectorDouble & vect)
{
	cout << "[";
	for (int i = 0; i < int(vect.size()); i++)
	{
		cout << vect(i);
		if (i < int(vect.size())-1) cout << " ";
	}
	cout << "]";
	return os;
}

ostream&
operator<< (ostream& os, const LaGenMatDouble & M)
{
	cout << "[ ";
	for (int p = 0; p < M.rows(); p++)
	{	
		for (int q = 0; q < M.cols(); q++)
		{
			cout << M(p,q) << " ";
			if (q < M.cols()-1) cout << " ";
		}
		if ( p< M.rows()-1)
			cout << ";" << endl;
	}
	cout << "]";
	return os;
}

