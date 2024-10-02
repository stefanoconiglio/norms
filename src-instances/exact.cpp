#include "exact.h"


LaGenMatDouble
exact::B(const LaGenMatDouble & A)
{
  int m = A.rows();
  int n = A.cols();

  assert (m > 0);

  //! initializing the needed matrixes
  LaGenMatDouble B(n, n);
  LaGenMatDouble A_trans(n, m);
  LaGenMatDouble I_m = LaGenMatDouble::eye(m, m);

  LaGenMatDouble PseudoHouseholder(m,m);
  PseudoHouseholder = double(-1)/m;
  for (int i = 0; i < m; i++)
    {
      PseudoHouseholder(i,i) = PseudoHouseholder(i,i)+1;
    }

  //! computation of the transposed matrix
  Blas_Mat_Trans_Mat_Mult(A, I_m, A_trans);

  //! final B computation
  B = A_trans * PseudoHouseholder * A;

  return B;
}

LaGenMatDouble
exact::B(const LaGenMatDouble & A, const LaVectorDouble & p)
{
  int m = A.rows();
  int n = A.cols();

  assert (m > 0);

  //! initializing the needed matrixes
  LaGenMatDouble B(n, n);
  LaGenMatDouble A_trans(n, m);
  LaGenMatDouble I_m = LaGenMatDouble::eye(m, m);

  //! create pp and pp_t vectors, which are ASSURED to be respectedly, a column and a row vector
  LaGenMatDouble pp;
  LaGenMatDouble pp_t;
  pp.resize(p.size(), 1);
  pp_t.resize(1, p.size());
  LaVectorDouble p_t(1, p.size());
  for (int i = 0; i < p.size(); i++)
    {
      pp(i,0) = p(i);
      pp_t(0,i) = p(i);
    }

  //! computation of the transposed matrix
  Blas_Mat_Trans_Mat_Mult(A, I_m, A_trans);

  LaGenMatDouble PseudoHousolder(m,m);
  PseudoHousolder = pp*pp_t;
  PseudoHousolder *= 1/( p * p );

  //! Multiplying by diag(p) so as to obtain PseudoHouseholder --> diag(p) * PseudoHolder * diag(p)
  for (int i = 0; i < m; i++)
    {
      if (p(i) == 1)
	continue;
      for (int j = 0; j < m; j++)
	{
	  if (p(j) == 1)
	    continue;
	  PseudoHousolder(i,j) = PseudoHousolder(i,j) * p(i);
	  PseudoHousolder(i,j) = PseudoHousolder(i,j) * p(j);
	}
    }

  //! final B computation
  B = A_trans * PseudoHousolder * A;

  return B;
}


/*
void
exact::plain_recalculate_l2_norm_hyperplane_parameters(Problem_instance & pi, int j) {
  assert (pi.problem_type.problem == HCLUSTERING);

  int m_j = pi.M(j);
  int n = pi.N();
  if (m_j == 0) {
    cout << "m_j = 0; returning" << endl;
    return;
  }

  LaGenMatDouble A = pi.A(j);

  //! in case A is squared or less than squared (we mean: the system Aw - gamma = 0 is under determined), exact solution (so that Aw-gamma = 0, w != 0) is computed; in case an infinity of solutions exist, the minimal 2-norm one is returned (via SVD, since we already have code to compute it.. =)
  //if (m_j <= n) {
  //  calculate_exact_fitting_parameters(pi, j);
  //  return;
  //}

  LaGenMatDouble B(n,n);
	
  B = exact::B(A);

  double lambda;

  if (pi.ops.power_tolerance == 0)
    Maths::SVDSmallestEigenVector(B, pi.w[j], lambda);
  else
    Maths::InversePowerMethod(B, pi.w[j], lambda, pi.ops.power_tolerance);

  pi.gamma[j] = exact::gamma(A, pi.w[j]);


  cout << "exact::plain_recalculate_l2_norm_hyperplane_parameters" << endl;
  cout << "pi.w[j] " << pi.w[j] << endl;
  cout << "pi.gamma[j] " << pi.gamma[j] << endl;
}

*/

/*
void
exact::recalculate_l2_norm_centroids(Problem_instance & pi, int j)
{
  pi.w[j] = pi.c(j);
}
*/

/*
void
exact::calculate_exact_fitting_parameters(Problem_instance & pi, int j) {
  
  abort ();

  cout << "exact::calculate_exact_fitting_parameters" << endl;

  LaGenMatDouble A = pi.A(j);
 
  cout << "A " << A << endl;

  LaVectorDouble e = LaGenMatDouble::ones(A.rows(),1);
  Maths::LaSVDLinearSolveIP(A, pi.w[j], e);
  pi.gamma[j] = 1;
  utility::normalize_2norm(pi.w[j], pi.gamma[j]);
  cout << "pi.w[j] " << pi.w[j] << endl;
  cout << "pi.gamma[j] " << pi.gamma[j] << endl;
}

*/


double
exact::gamma(const LaGenMatDouble & A, const LaVectorDouble & w)
{
  int m = A.rows();
  LaVectorDouble vec = LaGenMatDouble::ones(1,m) * A;
  double scalar = vec * w;
  return (scalar / m);
}


  /*
double
exact::gamma(const LaGenMatDouble & A, const LaVectorDouble & w, const LaVectorDouble & p)
{
  int m = A.rows();
  LaVectorDouble vec = LaGenMatDouble::ones(1,m) * A;
  double scalar = vec * w;
  return (scalar / (p*p));
}
  */


void
exact::update_l2_norm_hyperplane_parameters(Problem_instance & pi, int j)
{
  LaGenMatDouble A = pi.A(j);
  LaVectorDouble e = LaVectorDouble::ones(A.rows(),1);
  LaGenMatDouble B = exact::B(pi.A(j), e);
  double lambda;

  Maths::InversePowerMethod(B, pi.w[j], lambda, pi.ops.power_tolerance); // computes and sets w[j]
  pi.gamma[j] = exact::gamma(pi.A(j), pi.w[j]);
}




vector<bool>
exact::sample_by_policy(double sampling_threshold, vector<double> & row_probabilities, base_generator_type & generator )
{
  int m = row_probabilities.size();
  vector<bool> row_sampling(m, true);

  //!	sampling m times using the cumulative
  if (sampling_threshold == 0)
    {	
      vector<int> row_cumulative_sampling_counter = random_generator::sampleByCumulative(row_probabilities, m, generator );
      for (int i = 0; i < m; i++)
	{
	  if ( row_cumulative_sampling_counter[i] > 0 )
	    row_sampling[i] = false;
	}
    }
  //!	sampling using the average (multiplied by the threshold) as threshold
  else if (sampling_threshold > 0)
    {		
      double probability_threshold = 0;
      for (int i = 0; i < m; i++)
	{
	  probability_threshold += row_probabilities[i];
	}
      probability_threshold /= m;
      probability_threshold *= sampling_threshold;

      for (int i = 0; i < m; i++)
	{
	  if ( row_probabilities[i] > probability_threshold)
	    row_sampling[i] = false;
	}
    }
  //!	sampling using the number of points as threshold
  else if (sampling_threshold < 0)
    {
      int points_threshold = (int)((double)(m) * sampling_threshold * (-1));
      if (points_threshold > m)
	points_threshold = m;

      //! sorting probabilities
      vector<int> ordered_vector = utility::sort_gt(row_probabilities);
      for (int i = 0; i < points_threshold; i++)
	{
	  int idx = ordered_vector[ordered_vector.size()-i-1];
	  row_sampling[idx] = false;
	}
    }
  return row_sampling;
}


vector<double>
exact::exact_l2_norm_hyperplane_sampling_probabilities(const LaGenMatDouble & A)
{
#ifdef __DEBUG__
  cout << "exact_l2_norm_orthogonal_sampling_probabilities started" << endl;
#endif

  int 	m = A.rows();
  int	n = A.cols();

  LaGenMatDouble zeroed_A = A;
  LaVectorDouble weights(m);
  weights = 1;
  LaGenMatDouble B(n,n);
  LaVectorDouble eigenvector(n);
  vector<double> eigenvalues(m,-1);

  double true_eigenvalue;
  B =  exact::B(A, weights);
  Maths::InversePowerMethod(B, eigenvector, true_eigenvalue, 1e-6);
	
  for (int i = 1; i < m; i++)
    {
      if (i > 0)
	weights(i-1) = 1;
      if (i > 0)
	{
	  for (int l = 0; l < n; l++)
	    zeroed_A(i-1,l) = A(i-1,l);
	}
      for (int l = 0; l < n; l++)
	zeroed_A(i,l) = 0;
      B = exact::B(zeroed_A, weights);
      Maths::InversePowerMethod(B, eigenvector, eigenvalues[i], 1e-6);
    }

  vector<double> deltas(m,0);
  for (int i = 0; i< m; i++)
    {
      deltas[i] = true_eigenvalue - eigenvalues[i];
    }
  utility::normalize_1norm(deltas);
  return deltas;
}

vector<double>
exact::exact_l2_norm_regression_sampling_probabilities(const LaGenMatDouble & A, const LaVectorDouble & b)
{
#ifdef __DEBUG__
  cout << "exact_l2_norm_regression_sampling_probabilities started" << endl;
#endif

  int 	m = A.rows();
  int	n = A.cols();

  LaGenMatDouble zeroed_A = A;
  LaVectorDouble zeroed_b = b;
  LaVectorDouble x(n);

  double solution_value = Maths::LaSVDLinearSolveIP(A, x, b);

  vector<double> values(m,solution_value);
	
  for (int i = 1; i < m; i++)
    {
      if (i > 0)
	{
	  for (int l = 0; l < n; l++)
	    zeroed_A(i-1,l) = A(i-1,l);
	  zeroed_b(i-1) = b(i-1);
	}
      for (int l = 0; l < n; l++)
	zeroed_A(i,l) = 0;
      zeroed_b(i) = 0;
      values[i] -= Maths::LaSVDLinearSolveIP(zeroed_A, x, zeroed_b);
    }

  utility::normalize_1norm(values);
  return values;
}

/*
void
exact::recalculate_l2_norm_linear_submodel_parameters(Problem_instance & pi, int j, LaVectorDouble * p, bool flag_do_sampling, base_generator_type & generator )
{
  assert (pi.problem_type.problem == LINEAR_REGRESSION);
#ifdef __DEBUG__
  cout << "exact::recalculate_l2_norm_affine_submodel_parameters(Problem_instance & pi, int j, LaVectorDouble * p, bool flag_do_sampling)" << endl;
#endif

  //! variables needed for the weighted case
  bool flag_weighted_update;

  int m_j = pi.M(j);
  int n = pi.N();
  assert (m_j > 0);

  //! in case A is squared or less than squared (we mean: the system Aw - gamma = 0 is under determined), exact solution (so that Aw-gamma = 0, w != 0) is computed; in case an infinity of solutions exist, the minimal 2-norm one is returned (via SVD, since we already have code to compute it.. =)
  if (m_j <= n)
    {
      calculate_exact_fitting_parameters(pi, j);
      return;
    }
  //! nothing to be done if there are assigned points

  //! set the flag_weighted_update
  (p != NULL) ? flag_weighted_update = true : flag_weighted_update = false;

  LaGenMatDouble A = pi.A_lin_reg(j);
  LaVectorDouble b = pi.b_reg(j);

  
  //vector<bool> row_sampling(m_j, true);
  //vector<double> row_probabilities;
  //LaVectorDouble row_sampling_weights(m_j);

  //if (flag_do_sampling)
  //  {
  //    if (pi.ops.AOR.sampling == SVD)
//	row_probabilities = Maths::SVDSamplingProbabilities(A, b);
//
//      else if (pi.ops.AOR.sampling == EXACT)
//	row_probabilities = exact_l2_norm_regression_sampling_probabilities(A, b);
//
//      //! sampled according to the policy specified in double sampling_treshold
//      row_sampling = sample_by_policy(pi.ops.AOR.sampling_threshold, row_probabilities, generator );
//    }
  

  LaVectorDouble x(n);		// x_0, ..., x_{n-2} contains the first n-1 elements of w[j]; w[n-1] is equal to 1; x_{n-1} contains the value of gamma, which SHOULD BE NULL!

  //if (flag_do_sampling)
  //  Maths::LaSVDLinearSolveIP(A, x, b, row_sampling_weights);
  //else if (flag_weighted_update)
  if (flag_weighted_update)
    Maths::LaSVDLinearSolveIP(A, x, b, *p);
  else
    Maths::LaSVDLinearSolveIP(A, x, b);

#ifdef __DEBUG__
  cout << "A" << endl << A << endl;
  cout << "b" << endl << b << endl;
  cout << "x" << endl << x << endl;
#endif

  for (int l = 0; l < n-1; l++) {
    pi.w[j](l) = x(l);
  }
  pi.w[j](n-1) = -1;
  pi.gamma[j] = x(n-1);
  assert (pi.gamma[j] == 0);
}
*/

 /*

void
exact::recalculate_l2_norm_affine_submodel_parameters(Problem_instance & pi, int j, LaVectorDouble * p, bool flag_do_sampling, base_generator_type & generator )
{
  assert (pi.problem_type.problem == AFFINE_REGRESSION);
#ifdef __DEBUG__
  cout << "exact::recalculate_l2_norm_linear_submodel_parameters(Problem_instance & pi, int j, LaVectorDouble * p, bool flag_do_sampling)" << endl;
#endif

  //! variables needed for the weighted case
  bool flag_weighted_update;

  int m_j = pi.M(j);
  int n = pi.N();
  assert (m_j > 0);

  //! in case A is squared or less than squared (we mean: the system Aw - gamma = 0 is under determined), exact solution (so that Aw-gamma = 0, w != 0) is computed; in case an infinity of solutions exist, the minimal 2-norm one is returned (via SVD, since we already have code to compute it.. =)
  if (m_j <= n)
    {
      calculate_exact_fitting_parameters(pi, j);
      return;
    }
  //! nothing to be done if there are assigned points

  //! set the flag_weighted_update
  (p != NULL) ? flag_weighted_update = true : flag_weighted_update = false;


  LaGenMatDouble A = pi.A_aff_reg(j);
  LaVectorDouble b = pi.b_reg(j);

  
  vector<bool> row_sampling(m_j, true);
  vector<double> row_probabilities;
  LaVectorDouble row_sampling_weights(m_j);

  // if (flag_do_sampling)
  //   {
  //     if (pi.ops.AOR.sampling == SVD)
  // 	row_probabilities = Maths::SVDSamplingProbabilities(A, b);

  //     else if (pi.ops.AOR.sampling == EXACT)
  // 	row_probabilities = exact_l2_norm_regression_sampling_probabilities(A, b);

  //     //! sampled according to the policy specified in double sampling_treshold
  //     row_sampling = sample_by_policy(pi.ops.AOR.sampling_threshold, row_probabilities, generator );
	
  //     //! update of a row sampled matrix is performed as a WEIGHTED UPDATE with 0 weight for the sampled-out rows
  //     //CORRECT
  //     for (int i = 0; i < m_j; i++)
  // 	row_sampling_weights(i) = (row_sampling[i] == false) ? 0 : 1;
  //   }


  LaVectorDouble x(n);		// x_0, ..., x_{n-2} contains the first n-1 elements of w[j]; w[n-1] is equal to 1; x_{n-1} contains the value of gamma;

  //if (flag_do_sampling)
  //  Maths::LaSVDLinearSolveIP(A, x, b, row_sampling_weights);
  //else if (flag_weighted_update)
  if (flag_weighted_update)
    Maths::LaSVDLinearSolveIP(A, x, b, *p);
  else
    Maths::LaSVDLinearSolveIP(A, x, b);

#ifdef __DEBUG__
  cout << "A" << endl << A << endl;
  cout << "b" << endl << b << endl;
  cout << "x" << endl << x << endl;
#endif

  for (int l = 0; l < n-1; l++)
    pi.w[j](l) = x(l);
  pi.w[j](n-1) = -1;
  pi.gamma[j] = x(n-1);
}

 */

  /*
void
exact::recalculate_l2_norm_hyperplane_parameters(Problem_instance & pi, int j, LaVectorDouble * p, bool flag_do_sampling, bool flag_worsening, base_generator_type & generator )
{
  assert (pi.problem_type.problem == HCLUSTERING);
#ifdef __DEBUG__
  cout << "exact::recalculate_l2_norm_hyperplane_parameters(Problem_instance & pi, int j, LaVectorDouble * p, bool flag_do_sampling, bool flag_worsening)" << endl;
#endif

  //! variables needed for the weighted case
  bool flag_weighted_update;

  int m_j = pi.M(j);
  int n = pi.N();
  if (m_j == 0)
    return;

  LaGenMatDouble A = pi.A(j);

  //! in case A is squared or less than squared (we mean: the system Aw - gamma = 0 is under determined), exact solution (so that Aw-gamma = 0, w != 0) is computed; in case an infinity of solutions exist, the minimal 2-norm one is returned (via SVD, since we already have code to compute it.. =)
  if (m_j <= n)
    {
      calculate_exact_fitting_parameters(pi, j);
      return;
    }
  //! nothing to be done if there are assigned points

  //! set the flag_weighted_update
  (p != NULL) ? flag_weighted_update = true : flag_weighted_update = false;

  vector<bool> row_sampling(m_j, true);
  vector<double> row_probabilities;
  LaVectorDouble row_sampling_weights(m_j);


  // if (flag_do_sampling)
  //   {
  //     if (pi.ops.AOR.sampling == SVD)	// clustering problem consiered as an affine regression problem
  // 	{
  // 	  LaGenMatDouble A = pi.A_aff_reg(j);
  // 	  LaVectorDouble b = pi.b_reg(j);
  // 	  row_probabilities = Maths::SVDSamplingProbabilities(A, b);
  // 	}
  //     else if (pi.ops.AOR.sampling == EXACT)
  // 	{
  // 	  LaGenMatDouble A = pi.A(j);
  // 	  row_probabilities = exact_l2_norm_hyperplane_sampling_probabilities(A);
  // 	}
  //     //! sampled according to the policy specified in double sampling_treshold
  //     row_sampling = sample_by_policy(pi.ops.AOR.sampling_threshold, row_probabilities, generator );
	
  //     //! update of a row sampled matrix is performed as a WEIGHTED UPDATE with 0 weight for the sampled-out rows
  //     for (int i = 0; i < m_j; i++)
  // 	row_sampling_weights = (row_sampling[i] == false) ? 0 : 1;
  //   }


  LaGenMatDouble B(n,n);
	
  //if (flag_do_sampling)
  //  B = exact::B(A, row_sampling_weights);
  //else if (flag_weighted_update)
  if (flag_weighted_update)
    B = exact::B(A, *p);
  else
    B = exact::B(A);

  double lambda;

  if (pi.ops.power_tolerance == 0)
    {
      if (!flag_worsening)
	Maths::SVDSmallestEigenVector(B, pi.w[j], lambda);
      else
	Maths::SVDLargestEigenVector(B, pi.w[j], lambda);
    }
  else
    {
      if (!flag_worsening)
	Maths::InversePowerMethod(B, pi.w[j], lambda, pi.ops.power_tolerance);
      else
	Maths::PowerMethod(B, pi.w[j], lambda, pi.ops.power_tolerance);
    }

  //if (flag_do_sampling)
  //  pi.gamma[j] = exact::gamma(A, row_sampling_weights);
  //else if (flag_weighted_update)
  if (flag_weighted_update)
    pi.gamma[j] = exact::gamma(A, pi.w[j], *p);
  else
    pi.gamma[j] = exact::gamma(A, pi.w[j]);
}

*/
