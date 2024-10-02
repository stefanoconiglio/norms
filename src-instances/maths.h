#include"inclusions.h"

#ifndef __MATHS__
#define __MATHS__

class Maths
{
	public:
	//! computes the conditioning number via SVD decomposition (In Place SVD: matrix A is modified during the computation)
	static double cond(LaGenMatDouble A);
	
	//! solves a (generally overetermined) linear system via SVD decomposition (IP SVD, which modifies matrix A)
	//! the name mimics that of a Lapack++ function, though no such function is implementend in Lapack++
	static double LaSVDLinearSolveIP (const LaGenMatDouble & A, LaVectorDouble & x, const LaVectorDouble & b);

	//! same as the other one, but using weights in p
	static double LaSVDLinearSolveIP (const LaGenMatDouble & A, LaVectorDouble & x, const LaVectorDouble & b, const LaVectorDouble & p);

	//! computes the smallest eigenvalue and its respective eigenvalue by means of the SVD
	inline
	static void SVDSmallestEigenVector(const LaGenMatDouble & A, LaVectorDouble & v, double & lambda)
	{return SVDSmallestOrLargestEigenVector(A, v, lambda, true);}

	//! computes the largest eigenvalue and its respective eigenvalue by means of the SVD
	inline
	static void SVDLargestEigenVector(const LaGenMatDouble & A, LaVectorDouble & v, double & lambda)
	{return SVDSmallestOrLargestEigenVector(A, v, lambda, false);}

	//! returns the smallerst/largest eivenvalue/vector w.r.t. bool smallest == true or false
	static void
	SVDSmallestOrLargestEigenVector(const LaGenMatDouble & A, LaVectorDouble & v, double & lambda, bool smallest);
	
	//! updates the largest eigenvalue (and relative eigvenvector) with the inverse power method
	static void PowerMethod(LaGenMatDouble & a, LaVectorDouble & v, double & lambda, double toll);
	
	//! updates the smallest eigenvalue (and relative eigvenvector) with the inverse power method
	static void InversePowerMethod(LaGenMatDouble & a, LaVectorDouble & v, double & lambda, double toll);
	
	//! computes Drineas' SVD sampling probabilities
	static vector<double> SVDSamplingProbabilities(const LaGenMatDouble & A, const LaVectorDouble & b);

};

#endif


