#include "algorithm.h"

using namespace std;


void
algorithm::bemporads_criterion ( Problem_instance & pi ) {


  int m = pi.M();
  int k = pi.K();

  for ( int passes = 0; passes < pi.ops.bemporad_passes; passes++ ){

    cout << "Bemporad, pass " << passes << endl;

    //BEMPORAD: begin
    // loop over points
    // loop of clusters
    // mark all clusters with distance within DELTA
    // if > 1:
    // mark all c = 10 points closest to the current point
    // for each point, read its cluster
    // count the popularity of the clusters
    // reassign point to most popular
    // if ( pi.ops.bemporad_delta != 0 ) { //apply criterion
      // cout << "Bemporad kicking in!" << endl;
      double delta;
      // if ( pi.ops.bemporad_delta < 0 ){
      // 	delta = pi.solution_measure()/m;
      // }
      // else {
      // 	delta = pi.ops.bemporad_delta;
      // }
      cout << "delta " << abs (pi.ops.bemporad_delta) << endl;
      cout << "delta " << pi.ops.bemporad_delta << endl;
      cout << "sol measure " << pi.solution_measure() << endl;
      delta = pi.solution_measure()/m * abs (pi.ops.bemporad_delta); //abs needed: it is NEGATIVE if the criterion is used as a step
      cout << "delta " << delta << endl;
    // if ( 1 ) { //apply criterion
    //   double delta;
    //   delta = pi.solution_measure()/m;
      int c;
      cout << "pi.ops.bemporad_c " << pi.ops.bemporad_c << endl;
      if ( pi.ops.bemporad_c > 0 ) {
	c = pi.ops.bemporad_c;
      }
      else if ( pi.ops.bemporad_c < 0 ) {
	c = int(m*abs(pi.ops.bemporad_c));
      }
      cout << "c " << c << endl;
      if ( c < 3 ) {
	c = 3;
      }
      cout << "c " << c << endl;
      int undecidable_ctr = 0;
      for (int ii = 0; ii < m; ii++) {
	bool already_found = false;
	for (int j = 0; j < k; j++) {
	  if ( pi.distance[ii][j] < delta && already_found == true) {
	    //undecidable point found!
	    undecidable_ctr++;
	    // cout << "point " << ii << ", assigned to cluster " << pi.J(ii) << " is undecidable w.r.t. cluster " << j << endl;
	    vector <int> votes (k, 0);
	    vector <bool> checked (m, false);	  
	    for (int passes = 0; passes < c; passes++) {
	      // cout << "pass " << passes << endl;
	      double min_dist = 1e300;
	      int min_dist_point = -1;
	      for (int iii = 0; iii < m; iii++) {
		// cout << "Checking distance between point " << ii << " and point " << iii << endl;
		if (iii == ii) {
		  continue;
		}
		double dist = Formulas::l2_norm_prototypal_distance(pi.a[ii], pi.a[iii]);
		// cout << "dist " << dist << endl;
		if ( dist < min_dist && checked[iii] == false ) {
		  min_dist = dist;
		  min_dist_point = iii;
		}
	      } //end of for iii
	      votes[pi.J(min_dist_point)]++;
	      checked[min_dist_point] = true; 
	    } //end of passes
	    int most_voted_cluster = -1;
	    int max_votes = 0;
	    for (int jj = 0; jj < k; jj++) {
	      if ( votes[jj] > max_votes ) {
		max_votes = votes[jj];
		most_voted_cluster = jj;
	      }
	    }
	    pi.assignment[ii][pi.J(ii)] = false;
	    pi.assignment[ii][most_voted_cluster] = true;
	    cout << "Bemporad: reassigned point " << ii << " from cluster " << pi.J(ii) << " to cluster " << most_voted_cluster << endl;
	    break; //not there in previous experiments; should only speed stuff up: the voting does not change even if repeated (at most) k-1 times
	  }
	  else if ( pi.distance[ii][j] < delta && already_found == false) {
	    already_found = true;
	  }
	} //end of for j 
      }
      plain_parameter_update( pi );
      pi.update_distances();
      //BEMPORAD: end
    // }
      cout << "LogBemporad: undecidable_ctr = " << undecidable_ctr << endl;
  }
}

bool
algorithm::everything_tabu( Problem_instance & pi ) {

  int m = pi.M();
  int k = pi.K();

  for ( int i = 0; i < m; i++ ) {
    for ( int j = 0; j < k; j++ ) {
      // if ( pi.tabu_list.can_be_reassigned ( i, pi.distance[i], pi.J( i ) ) == true )
      //   return false;
      if ( pi.tabu_list.does_aspire( i,  j, pi.distance[i][j] ) == true )
        return false; //a nontabu pair is found
    }
  }
  return true; //if we never returned before, we didn't find a nontabu pair
}

Return_value
algorithm::compute_trivial_one_cluster_solution(Problem_instance & pi )
{
  for (int i = 0; i < pi.M(); i++)
    pi.assignment[i][0] = 1;
  plain_combinatorial_reassignment(pi);
  pi.update_distances();
  Return_value ret;
  ret.obj = pi.solution_measure();
  return ret;
}

void
algorithm::plain_combinatorial_reassignment(Problem_instance & pi, vector<bool> * skip) {

  for (int i = 0; i < pi.M(); i++) {

    if (skip != NULL)
      if ( (*skip)[i] == true ) {
        continue;
      }

    int old_j = pi.J(i);

    int closest_j = find_closest_nontabu_cluster( pi, i );

    if ( closest_j == old_j || closest_j == -1 ) { //if closest_j == -1, all reassignments must be tabu
      continue;
    }

    pi.assignment[i][old_j] = false;
    pi.assignment[i][closest_j] = true;
    // cout << "[" << i << "," << old_j << "]-->[" << i << "," << closest_j << "]" << endl;
  }
  return;
}

void
algorithm::plain_parameter_update(Problem_instance & pi, const int j ) {

  int m_j = pi.M(j);
  int n = pi.N();
  if (m_j == 0)
    return;

  if (pi.problem_type.problem == PCLUSTERING) {
    pi.w[j] = pi.c(j);
  }
  else if (pi.problem_type.problem == HCLUSTERING) {
    //if (m_j < n) {
    //underdetermined system: exact or min 2-norm sol being computed
    //  exact::calculate_exact_fitting_parameters(pi, j);
    //  return;
    //}
    LaGenMatDouble A = pi.A(j);
    LaGenMatDouble B = exact::B(A);

    double lambda;

    if (pi.ops.power_tolerance == 0)
      Maths::SVDSmallestEigenVector(B, pi.w[j], lambda);
    else
      Maths::InversePowerMethod(B, pi.w[j], lambda, pi.ops.power_tolerance);
    pi.gamma[j] = exact::gamma(A, pi.w[j]);
  }
  else if (pi.problem_type.problem == LINEAR_REGRESSION) {
    LaGenMatDouble A = pi.A_lin_reg(j);
    LaVectorDouble b = pi.b_reg(j);

    LaVectorDouble x(n);		// x_0, ..., x_{n-2} contains the first n-1 elements of w[j]; w[n-1] is equal to 1; x_{n-1} contains the value of gamma, which SHOULD BE NULL!

    Maths::LaSVDLinearSolveIP(A, x, b);

    for (int l = 0; l < n-1; l++) {
      pi.w[j](l) = x(l);
    }

    pi.w[j](n-1) = -1;
    pi.gamma[j] = x(n-1);
    assert(pi.gamma[j] == 0);
  }
  else if (pi.problem_type.problem == AFFINE_REGRESSION) {
    LaGenMatDouble A = pi.A_aff_reg(j);
    LaVectorDouble b = pi.b_reg(j);

    LaVectorDouble x(n);		// x_0, ..., x_{n-2} contains the first n-1 elements of w[j]; w[n-1] is equal to 1; x_{n-1} contains the value of gamma;

    Maths::LaSVDLinearSolveIP(A, x, b);

    for (int l = 0; l < n-1; l++)
      pi.w[j](l) = x(l);

    pi.w[j](n-1) = -1;
    pi.gamma[j] = x(n-1);
  }
  else if ( pi.problem_type.problem == PAMF ) {
    pamf::Linear1normRegressionParameterUpdate( pi, j);
  }
}



void
algorithm::plain_parameter_update(Problem_instance & pi ) {

  for (int j = 0; j < pi.K(); j++)
    plain_parameter_update( pi, j );
}



int
algorithm::find_closest_nontabu_cluster( Problem_instance & pi, int i ) {

  // cout << " algorithm::find_closest_nontabu_cluster( Problem_instance & pi, int i ) " << endl;
  double min_dist = 1e300;
  int closest_cluster = -1;

  // finding the non tabu closest cluster
  for (int j = 0; j < pi.K(); j++) {
    // cout << "pi.ops.TS.flag_use_TS " << pi.ops.TS.flag_use_TS << endl;
    // cout << "pi.tabu_list.is_tabu(i,j) " << pi.tabu_list.is_tabu(i,j) << endl;
    // cout << "pi.tabu_list.does_aspire(i,j,pi.distance[i][j]) " << pi.tabu_list.does_aspire(i,j,pi.distance[i][j]) << endl;

    if (pi.ops.TS.flag_use_TS && pi.tabu_list.is_tabu(i,j) && pi.tabu_list.does_aspire(i,j,pi.distance[i][j]) == false) {
      // cout << "tabu pair found" << endl;
      continue;
    }
    else if (pi.distance[i][j] < min_dist) {
      min_dist = pi.distance[i][j];
      closest_cluster = j;
    }
  }

  // cout << "returning: point " << i << " will be reassigned to cluster " << closest_cluster << endl;
  return closest_cluster;
}


void
algorithm::degenerancy_avoider( Problem_instance & pi, int j, base_generator_type & generator ) {

  int n = pi.N();

  if (pi.M(j) == 0) {
    cout << "Empty cluster #" << j << " found: parameters randomly regenerated" << endl;

    random_generator::randomize(pi.w[j], generator );

    if (pi.problem_type.problem == AFFINE_REGRESSION || pi.problem_type.problem == LINEAR_REGRESSION)
      pi.w[j](n-1) = -1;

    pi.gamma[j] = random_generator::random( generator );

    if (pi.problem_type.problem == LINEAR_REGRESSION)
      pi.gamma[j] = 0;

    pi.update_distances( j );
  }
}
