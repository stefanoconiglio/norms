#include "maths.h"

extern int extended_output;

bool restricted_drineas_prob = false;	//! if true, only the first inequality of the original paper (2006) is used to define the sampling probabilities


// used to keep track of the number of times the methods are called; it was needed to make some comparison with the SVD method for computing the eigenvalue, or something like that
int nLeastEigPowersCalled;
int nLeastEigPowersIterationsCalled;

void
Maths::SVDSmallestOrLargestEigenVector(const LaGenMatDouble & A, LaVectorDouble & v, double & lambda, bool smallest)
{
  assert (A.rows() == A.cols());
  int m = A.rows();
	
  //! new construction of B() and eigenvalues computation
  LaVectorDouble real_eigenvals(m);
  LaVectorDouble img_eigenvals(m);
  LaGenMatDouble eigenvectors(m,m);
  LaEigSolve(A, real_eigenvals, img_eigenvals, eigenvectors);

  //! fetching the min eigenval from the eigenval 'matrix'
  int min_eigenval_idx = 0;
  int max_eigenval_idx = 0;
  int chosen_eigenval_idx = 0;
  double min_eigenval_value = 1e300;
  double max_eigenval_value = -1e300;

  for (int i = 0; i < m; i++)
    {			
      if (real_eigenvals(i) < min_eigenval_value)
	{	
	  min_eigenval_idx = i;
	  min_eigenval_value = real_eigenvals(i);
	}

      if (real_eigenvals(i) > max_eigenval_value)
	{	
	  max_eigenval_idx = i;
	  max_eigenval_value = real_eigenvals(i);
	}
    }

  chosen_eigenval_idx = smallest ? min_eigenval_idx : max_eigenval_idx;

  //! writing the optimal w values
  for (int i = 0; i < m; i++)
    v(i) = eigenvectors(i, chosen_eigenval_idx);
	
#ifdef __DEBUG__
  //! asserting the minimality of the chosen eigenvalue
  for (int d = 0; d < D; d++)
    {
      assert (min_eigenval_value <= real_eigenvals(d));
      assert (max_eigenval_value >= real_eigenvals(d));
    }
#endif
}

double 
Maths::cond(LaGenMatDouble A)
{
  double cond;
  int min_dim = A.rows() < A.cols() ? A.rows() : A.cols();
  LaVectorDouble SIGMA(min_dim);
  LaGenMatDouble U = LaGenMatDouble(A.rows(), A.rows());
  LaGenMatDouble VT = LaGenMatDouble(A.cols(), A.cols());
  LaSVD_IP(A, SIGMA, U, VT);
  cond = SIGMA(0) / SIGMA(min_dim-1);
  return cond;
}

double
Maths::LaSVDLinearSolveIP(const LaGenMatDouble & A, LaVectorDouble & x, const LaVectorDouble & b, const LaVectorDouble & p)
{
  LaGenMatDouble A_weighted(A);
  LaVectorDouble b_weighted(b);
  for (int i = 0; i < A.rows(); i++)
    {
      for (int j = 0; j < A.cols(); j++)
	{
	  A_weighted(i,j) = A_weighted(i,j) * p(i);
	}
      b_weighted(i) = b_weighted(i) * p(i);
    }
  return (LaSVDLinearSolveIP(A_weighted, x, b_weighted));
}

double
Maths::LaSVDLinearSolveIP (const LaGenMatDouble & A, LaVectorDouble & x, const LaVectorDouble & b)
{
  //! initializations: dimensions
  int m = A.rows();
  int n = A.cols();

  //! initializing the needed matrixes
  LaGenMatDouble A_backup(A);		//! reminder: LaSVD_IP destroys the A matrix
  LaGenMatDouble I_m = LaGenMatDouble::eye(m,m);
  LaGenMatDouble I_n = LaGenMatDouble::eye(n,n);
  LaGenMatDouble U(m,m);
  LaGenMatDouble UT(m,m);
  LaGenMatDouble VT(n,n);
  LaGenMatDouble V(n,n);
  int dim_diag_SIGMA = m < n ? m : n; //! length of the diagonal of the SIGMA matrix
  LaVectorDouble diag_SIGMA = LaVectorDouble(dim_diag_SIGMA) ;
  LaGenMatDouble SIGMA = LaGenMatDouble(m,n);

  LaVectorDouble overlined_b = LaVectorDouble(m);
  LaVectorDouble overlined_x = LaVectorDouble(n);

  //! initializing x... one never knows what ppl passes to functions...
  x = LaVectorDouble(n);

  //! computing the SVD
  LaSVD_IP (A_backup,diag_SIGMA,U,VT);
  for (int i = 0; i < diag_SIGMA.size(); i++)
    if (diag_SIGMA(i) < 1e-10)
      diag_SIGMA(i) = 0;

  //! construct the full SIGMA matrix
  int idx = 0;
  for (int i = 0; i<m; i++)
    for (int j = 0; j < n; j++)
      //if (idx <= diag_SIGMA.size())
        SIGMA (i,j) = (i != j) ? 0 : diag_SIGMA (idx++);
  //else
  //SIGMA (i,j) = 0;

  //! compute UT
  Blas_Mat_Trans_Mat_Mult ( U, I_m, UT );

  //! compute V
  Blas_Mat_Trans_Mat_Mult ( VT, I_n, V );

  //! compute overlined_b
  overlined_b = UT*b;

  //! solve the linear system SIGMA * overlined_x = overlined_b
  for (int i = 0; i < diag_SIGMA.size(); i++)
    {
      if (diag_SIGMA(i) == 0)
	overlined_x(i) = 0;
      else
	overlined_x(i) = overlined_b(i) / diag_SIGMA(i);
    }

  //BOFFO: diag_SIGMA is a min(m,n) vector; since the system is underdetermined, m < n and hence diag_SIGMA is a m-component vector;
  //since x is a n-component vector, it's first m components will equal overlined_b(i)/sigma(i); the remaining n-m will take zero value
  //so as to produce a minimal norm vector!
  
  for (int i = diag_SIGMA.size(); i<n; i++)
    overlined_x(i) = 0;
	
  //! solve the linear system VT * x = overlined_x
  LaLinearSolve (VT, x, overlined_x);

  //! compute the 2-norm of the residual vector
  return (Blas_Norm2 (A*x - b));
}

void
Maths::InversePowerMethod(LaGenMatDouble & A, LaVectorDouble & q, double & lambda, double toll)
{
  double		error = 1e300;
  LaVectorDouble	z(q.size());
  LaVectorDouble	old_q(q.size());

  //! If A (A_k) is square, since the pseudohouseholder matrix is SINGULAR, it follows that
  //! A^T * pseudohouseholder * A is SINGULAR as well. (It doens't happen, contrarily to what Mr Frontini
  //! said, in case A is of maximum rank and rectangular).

  //! As B is singular, its smallest eigenvalue is clearly NULL and 
  //! of course B has a corresponding eigenvector which is NOT null, as it is part of
  //! its null space. So asking the Power Method to compute it is not foolish.
  //! The problem is that (clearly I don't know why) the Power Methods shows not a single
  //! convergent sequence of vector q, but TWO different convergent subsequences of such vectors
  //! Our task is then in choosing only ONE of the two convergent subsequences, and use the final
  //! vector as the eigenvector we are looking for.
  //! It looks like that such sequences differ only in a sign; then, since we do not know INSIDE
  //! this function wheter A (A_k) is square or not, we can simply take the first element of q and
  //! change its sign in case it is negative; ALWAYS.

  nLeastEigPowersCalled++;

  //! as is should be, normalize q
  q *= 1/Blas_Norm2(q);

  while (error > toll)
    {
      old_q = q;
      LaLinearSolve (A, z, q);

      q = z * (1/Blas_Norm2(z));

      //! sign change
      if (q(0) < 0)
	q = q * (-1);

      error = Blas_Norm2(q - old_q);

      nLeastEigPowersIterationsCalled++;
    }

  lambda = q * z;
  lambda = 1 / lambda;
}

void
Maths::PowerMethod(LaGenMatDouble & A, LaVectorDouble & q, double & lambda, double toll)
{
  double		error = 1e300;
  LaVectorDouble	z(q.size());
  LaVectorDouble	old_q(q.size());

  //! If A (A_k) is square, since the pseudohouseholder matrix is SINGULAR, it follows that
  //! A^T * pseudohouseholder * A is SINGULAR as well. (It doens't happen, contrarily to what Mr Frontini
  //! said, in case A is of maximum rank and rectangular).

  //! As B is singular, its smallest eigenvalue is clearly NULL and 
  //! of course B has a corresponding eigenvector which is NOT null, as it is part of
  //! its null space. So asking the Power Method to compute it is not foolish.
  //! The problem is that (clearly I don't know why) the Power Methods shows not a single
  //! convergent sequence of vector q, but TWO different convergent subsequences of such vectors
  //! Our task is then in choosing only ONE of the two convergent subsequences, and use the final
  //! vector as the eigenvector we are looking for.
  //! It looks like that such sequences differ only in a sign; then, since we do not know INSIDE
  //! this function wheter A (A_k) is square or not, we can simply take the first element of q and
  //! change its sign in case it is negative; ALWAYS.

  nLeastEigPowersCalled++;
  //! as is should be, normalize q
  q *= 1/Blas_Norm2(q);

  while (error > toll)
    {
      old_q = q;
      z = A * q;

      q = z * (1/Blas_Norm2(z));

      //! sign change
      if (q(0) < 0)
	q = q * (-1);

      error = Blas_Norm2(q - old_q);
		
      nLeastEigPowersIterationsCalled++;
    }

  lambda = q * z;
  lambda = 1 / lambda;
}

vector<double>
Maths::SVDSamplingProbabilities(const LaGenMatDouble & A, const LaVectorDouble & b)
{
#ifdef __DEBUG__
  cout << "SVD_l2_norm_regression_with_sampling started" << endl;
#endif

  //! final probabilities values
  double			add_1, add_2, add_3;
  double			beta_1, beta_2, beta_3;

  //! sampling helper structures
  LaGenMatDouble		S;				//! sampling matrix
	
  //! initializing matrix dimensions
  int m = A.rows();
  int n = A.cols();

  //! m >> n is needed in order to sample
  assert (m > n);

  //! sampling report variables
  vector<double>	obj_function(m, 0);
  vector<int> ordered_vector;

  //! dummy variable x
  LaVectorDouble x(n);
	
  //! SVD matrices
  LaVectorDouble diag_SIGMA(n);
  LaGenMatDouble U(m, m);
  LaGenMatDouble VT(n, n);
  LaGenMatDouble I_m = LaGenMatDouble::eye(m, m);
  LaGenMatDouble UA(m, n);
  LaVectorDouble UAT(n, m);

  //! sampling variables
  vector<double> row_probabilities(m, 0);
  LaGenMatDouble UAc(m, m - n);
  LaGenMatDouble UAcT(m - n, m);
  LaVectorDouble UAc_UAcT_b(m);
  double U_sum_row_norm_UAc_UAcT_b = 0;
  double U_sum_row_norm_2 = 0;
  double U_sum_row_UAc_UAcT_b_2 = 0;
	
  //! backup matrix creation
  LaGenMatDouble A_backup(A);
	
  //! computing the SVD of A
  LaSVD_IP(A_backup, diag_SIGMA, U, VT);
	
  //! Splitting the square matrix U into the thin-SVD m x (m-n) matrix UA and its complementary UAc
  for (int i = 0; i < m; i++)
    {
      //! first half goes into UA
      for (int j = 0; j < n; j++)
	UA(i, j) = U(i, j);

      //! second half goes into UAc
      for (int j = n; j < m; j++)
	UAc(i, j-n) = U(i, j);			
    }

  //! computing UAT
  Blas_Mat_Trans_Mat_Mult(UA, I_m, UAT);

  //! computing UAcT
  Blas_Mat_Trans_Mat_Mult(UAc, I_m, UAcT);

  //! computing the UAc * UAcT matrix, used to evaluate the projection of b onto the space complementary to the column space
  UAc_UAcT_b = UAc * UAcT * b;

  //! compute some parameters based on the sum of 'row information'
  for (int i = 0; i < m; i++)
    {
      U_sum_row_norm_2 += pow( Blas_Norm2( U.row(i) ), 2 );	
      U_sum_row_norm_UAc_UAcT_b += ( Blas_Norm2( U.row(i) ) * UAc_UAcT_b(i) );
      U_sum_row_UAc_UAcT_b_2 += pow( UAc_UAcT_b(i), 2 );
    }

  //! computing the sampling probabilities all the probabilities
  for (int i = 0; i < m; i++)
    {
      double U_row_norm = Blas_Norm2( U.row(i) );
      double U_row_norm_2 = pow( U_row_norm, 2 );

      add_1 = U_row_norm_2 / U_sum_row_norm_2;
		
      add_2 = U_row_norm  * UAc_UAcT_b(i) / U_sum_row_norm_UAc_UAcT_b;

      //! since add_2 could be < 0 and this would violate the necessary conditions indicated in the paper, we force it to be 0 in case it is found negative
      if (add_2 < 0)
	add_2 = 0;

      add_3 = pow( UAc_UAcT_b(i), 2 ) / U_sum_row_UAc_UAcT_b_2;

      beta_1 = (double)1/3;
      beta_2 = (double)1/3;
      beta_3 = (double)1/3;

      //! computing the final sampling probability
      if (!restricted_drineas_prob)
	{
	  row_probabilities[i] = beta_1 * add_1 + beta_2 * add_2 + beta_3 * add_3;
	}
      else
	{
	  row_probabilities[i] = add_1;
	}
			
    }
  return row_probabilities;
}

