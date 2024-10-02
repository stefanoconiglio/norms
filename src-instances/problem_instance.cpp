#include "problem_instance.h"

using namespace std;

int
Problem_instance::M(int j) {
  int m_j = 0;
  for (int i = 0; i < M(); i++) {
    if (assignment[i][j] == true)
      m_j++;
  }
  return m_j;
}

int
Problem_instance::J(int i)
{
  for (int jj = 0; jj < K(); jj++) {
    if (assignment[i][jj] == true)
      return jj;
  }
  return -1;
}

Problem_instance::Problem_instance(int m, int n, int k)
{
  file_name = "Empty";

  LaVectorDouble dummy_LaVectorDouble(n);
  dummy_LaVectorDouble = -1;
  vector<double> dummy_double_vector(k, -1);
  vector<bool> dummy_bool_vector(k,false);	
  a = vector<LaVectorDouble>(m, dummy_LaVectorDouble);
  w = vector<LaVectorDouble>(k, dummy_LaVectorDouble);
  gamma = vector<double>(k, -1);
  W = vector<LaVectorDouble>(k, dummy_LaVectorDouble);
  G = vector<double>(k, -1);
	
  assignment = vector<vector<bool> >(m, dummy_bool_vector);
  distance = vector<vector<double> >(m, dummy_double_vector);	

}

Problem_instance::Problem_instance(string & path, Problem_type & type, t_alg_options & ops, int external_k, int sample_ratio)
{
  problem_type = type;
  this->ops = ops;
  file_name = path;

  ifstream dataset(path.c_str());
  assert(dataset.good());

  string tmp;

  int m, n, k;
  int actual_m = 0;

  dataset >> tmp;    //! reads the title line
  dataset >> tmp;    //! reads "m: _"
  dataset >> m;

  dataset >> tmp;    //! reads "n: _"
  dataset >> n;

  LaVectorDouble dummy_point(n);

  dataset >> tmp;    //! reads "k: _"
  dataset >> k;

  do {
    dataset >> tmp;	//! reads "Data_set:"
  } while (strcmp(tmp.c_str(), "Data_set:"));

  // 	while (!dataset.eof())		// for some reason this reads ONE MORE LINE...
  for (int ctr = 0; ctr < m; ctr++)
    {	
      double current_coordinate;
      dataset >> tmp;	//! reads "p:_"
      for (int l = 0; l < n; l++)
	{
	  dataset >> current_coordinate;
	  dummy_point(l) = current_coordinate;
	}
      if ( (ctr % sample_ratio) == 0)
	{
	  a.push_back(dummy_point);
	  actual_m++;
	}
      else
	continue;
    }

  m = actual_m;

  dataset.close();

  // creating a set of k/external_k random clusters
  if (external_k != -1)
    k = external_k;

  LaVectorDouble dummy_vector(n);
  dummy_vector = -1;
  double dummy_scalar = -1;

  for (int j = 0; j < k; j++)
    {
      w.push_back(dummy_vector);
      gamma.push_back(dummy_scalar);
      W.push_back(dummy_vector);
      G.push_back(dummy_scalar);
    }

  // initializing distance and assignment data structures
  vector<bool> dummy_assignment_vector(k, false);
  dummy_assignment_vector[0] = true;		// by default all points are associated to cluster 0 (it's convenient to avoid the function int j(int i) to return inconsistent values...
  for (int i = 0; i < m; i++)
    assignment.push_back(dummy_assignment_vector);

  std::vector<double> dummy_distance_vector(k, 1e300);
  for (int i = 0; i < m; i++)
    distance.push_back(dummy_distance_vector);

  assert(int(a.size()) == M());
  assert(int(assignment.size()) == M());
  assert(int(distance.size()) == M());

  // initializing the tabu list
  tabu_list = full_tabu_list(ops.TS.length, M(), K(), ops.TS.flag_use_aspiration_criterion);

}

void
//Problem_instance::randomly_populate_clusters(int * ext_pred)
Problem_instance::randomly_populate_clusters( base_generator_type & generator )
{
  bool flag_do_regression = false;
  bool flag_gamma_fixed_to_0 = false;

  if (problem_type.problem == AFFINE_REGRESSION || problem_type.problem == LINEAR_REGRESSION || problem_type.problem == PAMF )
    flag_do_regression = true;

  if (problem_type.problem == LINEAR_REGRESSION)
    flag_gamma_fixed_to_0 = true;

  int n = N();
  int k = K();

  LaVectorDouble new_w(n);
  double new_gamma;

  for (int j = 0; j < k; j++) {
    //random_generator::randomize(w[j], ext_pred);
    random_generator::randomize( w[j], generator );
    
    //new_gamma = flag_gamma_fixed_to_0 ? 0 : random_generator::random(ext_pred);
    new_gamma = flag_gamma_fixed_to_0 ? 0 : random_generator::random( generator );
    
    if (flag_do_regression)
      w[j](n-1) = -1;
    
    gamma[j] = new_gamma;
    
    if (problem_type.problem == HCLUSTERING)
      utility::normalize_2norm(w[j], gamma[j]);
  }
}

ostream& operator<< (ostream& os, Problem_instance & ps)
{
  cout << "Objective " << ps.solution_measure() << endl;

  // cout << "Points" << endl;
  // for (int i = 0; i < ps.M(); i++) {		
  //   cout << "a[" << i << "] = " << ps.a[i];
  //   cout << endl;
  // }

  //! early returning if no clusters have been defined
  if (ps.K() == 0)
    {
      cout << "No clusters defined" << endl;
      return os;
    }

  cout << "Clusters" << endl;
  for (int i = 0; i < ps.K(); i++)
    {
      cout << "w[" << i << "] = " << ps.w[i] << " ";
      cout << "gamma[" << i << "] = " << ps.gamma[i] << " ";
      cout << "W[" << i << "] = " << ps.W[i] << " ";
      cout << "G[" << i << "] = " << ps.G[i] << endl;
    }

  cout << "Assignment" << endl;

  for (int i = 0; i < int(ps.assignment.size()); i++)
    cout << i << ": " << ps.assignment[i] << endl;

  cout << "Distances" << endl;

  for (int i = 0; i < int(ps.distance.size()); i++)
    cout << i << ": " << ps.distance[i] << endl;

  cout << "Objective " << ps.solution_measure() << endl;

  return os;
}

void
Problem_instance::output_solution_to_file(string & file_name)
{
  ofstream file( file_name.c_str() );

  for (int j = 0; j < K(); j++) {
    for (int l = 0; l < N(); l++)
      file << w[j](l) << ", ";
    file << gamma[j] << endl;
  }
  
  if ( problem_type.problem == PAMF )
    for (int j = 0; j < K(); j++) {
      for (int l = 0; l < N(); l++)
	file << W[j](l) << ", ";
      file << G[j] << endl;
    }

}


void
Problem_instance::output_solution_to_matlab(string & file_name, bool flag_print_lines)
{
  //! length variables
  int m = M();
  int n = N();
  int k = K();

  assert(n == 2);

  //! data structure variables
  int j_i;
	
  //! variables needed for the min/max computation
  vector<double> max(n,-1e300);
  vector<double> min(n,1e300);

  //! matlab output variables
  char str_representation_x0[256];
  char str_postprocessing_representation[2024];
  vector<string> assignment_string;
  char str_coord[256];
  char str_plot_command[256];
  string name;
  char str_j[256];

  //! opening the output file
  ofstream matlab_output(file_name.c_str());

  //! set the string vector
  for (int j = 0; j < k; j++ )
    {
      sprintf(str_j, "%i", j);
      assignment_string.push_back (string ("c") + string (str_j) + string (" = [") );
    }

  for (int i = 0; i < m; i++)
    {
      j_i = J(i);
      for (int l = 0; l < n; l++)
	{
	  sprintf (str_coord, "%f", a[i](l));
	  assignment_string[j_i] += string (str_coord) + string (" ");
	}
      assignment_string[j_i] += string(";\n");
    }

  //! close the string vector
  for (int j = 0; j < k; j++ )
    assignment_string[j] += string ("];");

  //! find the maximum and minimum per every dimension
  get_min_max_dimensions(min, max);

  //! print the presentation x
  sprintf(str_representation_x0, "x0 = %f:.01:%f;", min[0], max[0]);
  matlab_output << str_representation_x0 << endl;

  //! print the hyperplane equations
  for (int j = 0; j < k; j++ )
    {
      matlab_output << assignment_string[j] << endl;
      name = "y_c";
      sprintf(str_j, "%i", j);
      name += str_j;
      matlab_output << name << " = ";
      for (int l = 0; l < n-1; l++)
	matlab_output << " - " << w[j](l) << " / " << w[j](n-1) << " * x" << l << " ";

      matlab_output << " + " << gamma[j] << " / " << w[j](n-1) << ";" << endl;
    }

  matlab_output << "colors = jet(" << k << ");" << endl;

  //! plot the point sets, each one with a different color (thanks to octave, this is totally free)
  for (int j = 0; j < k; j++)
    {
      //sprintf(str_plot_command, "if (size(c%i) >0) plot (c%i(:,1), c%i(:,2), '@'); end", j, j, j);
      sprintf(str_plot_command, "if (size(c%i) >0) plot (c%i(:,1), c%i(:,2),'o', 'color', colors(%i,:)); end", j, j, j, j+1);
      matlab_output << str_plot_command << endl;
      matlab_output << "hold on;" << endl;
    }

  if (flag_print_lines)
    {
      //! postprocess the y_c vectors
      for (int j = 0; j < k; j++)
	{
	  sprintf(str_postprocessing_representation, "for i = 1 : size(y_c%i, 2) \n if y_c%i(i) > %f \n y_c%i(i) = %f; \n end \n if y_c%i(i) < %f \n y_c%i(i) = %f; \n end \n end", j, j, max[1], j, max[1], j, min[1], j, min[1]);
	  matlab_output << str_postprocessing_representation << endl;
	}
	
      if (problem_type.problem != PCLUSTERING)
	{
	  //! plot the hyperplanes
	  for (int j = 0; j < k; j++)
	    {
	      sprintf(str_plot_command, "plot(x0, y_c%i, 'color', colors(%i,:));", j, j+1);
	      matlab_output << str_plot_command << endl;
	      matlab_output << "hold on;" << endl;
	    }
	}
      else if (problem_type.problem == PCLUSTERING)
	{
	  //! plot the centroids
	  for (int j = 0; j < k; j++)
	    {
	      matlab_output << "centroid_" << j << " = [";
	      for (int l = 0; l < n; l++)
		matlab_output << w[j](l) << " ";
	      matlab_output << "];" << endl;

	      matlab_output << "plot(centroid_" << j << "(1,1), centroid_" << j << "(1,2), 'o', 'color', colors(" << j+1 << ",:));" << endl;
	      matlab_output << "hold on;" << endl;
	    }
	}
    }
}

void
Problem_instance::output_standard_ampl_data(string & file_name)
{
  ofstream ampl_data(file_name.c_str());

  int m = M();
  int n = N();
  int k = K();

  ampl_data << "set M := "; for (int i = 1; i <= m; i++) ampl_data << i << " "; ampl_data << ";" << endl;
  ampl_data << "set N := "; for (int i = 1; i <= n; i++) ampl_data << i << " "; ampl_data << ";" << endl;
  ampl_data << "set K := "; for (int i = 1; i <= k; i++) ampl_data << i << " "; ampl_data << ";" << endl;

  ampl_data << "param a : \n";
  for (int i = 1; i <= n; i++)
    ampl_data << " " << i;
  ampl_data << ":=" << endl;
  for (int i = 0; i < m; i++)
    {
      ampl_data << i+1 << " ";
      for (int l = 0; l < n; l++)
	ampl_data << " " << a[i](l);
      ampl_data << "\n";
    }
  ampl_data << ";" << endl;
}


void
Problem_instance::print_dat_file_assignment(string & instance_file_name)
{
  string file_name(instance_file_name);
  file_name += "_initial_pool.dat";
  ofstream file(file_name.c_str());

  file << "let k := " << K() << ";" << endl;

  /*
  file << "param IS(tr):";
  for (int j = 0; j < K(); j++)
    file << j+1 << " ";
  file << ":=" << endl;

  for (int i = 0; i < M(); i++)
    {
      file << i+1 << " ";
      for (int j = 0; j < K(); j++)
	{
	  if (assignment[i][j] == true)
	    file << "1 ";
	  else
	    file << "0 ";
	}
      file << endl;
    }
  file << ";" << endl;
  */
  
  for (int i = 0; i < M(); i++)
    {
      for (int j = 0; j < K(); j++)
	{
	  file << "let IS[" << i+1 << "," << j+1 << "] := ";
	    if (assignment[i][j] == true)
	      file << "1 ";
	    else
	      file << "0 ";
	  file << ";" << endl;
	}
    }
  
}


void
Problem_instance::print_dat_file_distances(string &instance_file_name)
{
  string file_name(instance_file_name);
  file_name += "_initial_solution_dist.dat";
  ofstream file(file_name.c_str());

  for (int i = 0; i < M(); i++)
    for (int j = 0; j < K(); j++)
      file << "let Dist[" << i+1 << "," << j+1 << "] := " << point_measure(i,j) << ";" << endl;
}

void
Problem_instance::print_dat_file_SNOPT_representation(string & file_path)
{
  for (int j = 0; j < K(); j++)
    {
      string cluster_file_name(file_path);
      char str_j[15];
      sprintf(str_j, "_cl%i", j);
      string cluster_idx(str_j);

      cluster_file_name += cluster_idx;

      ofstream cluster_file(cluster_file_name.c_str());

      for (int i = 0; i < M(); i++)
	{
	  cluster_file << "let d[" << i+1 << "] := ";
	  if (assignment[i][j]  == true)
	    cluster_file << "1;" << endl;
	  else
	    cluster_file << "0;" << endl;

	  cluster_file << "let Dist[" << i+1 << "] := " << point_measure(i,j) << ";" << endl;
	}
      for (int l = 0; l < N(); l++)
	{
	  cluster_file << "let w[" << l+1 << "] := " << w[j](l) << ";" << endl;
	}
      cluster_file << "let w0 := " << gamma[j] << ";" << endl;
    }
}


void
Problem_instance::print_binary_representation(string & file_path)
{
  ofstream binary_representation_file(file_path.c_str());

  for (int i = 0; i < M(); i++)
    {
      for (int j = 0; j < K(); j++)
	{
	  binary_representation_file << "x [" << i+1 << "," << j+1 << "] = ";
	  if (assignment[i][j] == true)
	    binary_representation_file << "1;" << endl;
	  else
	    binary_representation_file << "0;" << endl;
	}
    }
  double max_of_max_of_l2_norm_squared_distances = 0;
  for (int j = 0; j < K(); j++)
    {
      binary_representation_file << "Eps(" << j << "):";
      double max_of_l2_norm_squared_distance = 0;
      for (int i = 0; i < M(); i++)
	{
	  if (assignment[i][j] != true)
	    continue;
	  if (distance[i][j] > max_of_l2_norm_squared_distance)
	    max_of_l2_norm_squared_distance = pow(distance[i][j],2);
	}
      if (max_of_l2_norm_squared_distance > max_of_max_of_l2_norm_squared_distances)
	max_of_max_of_l2_norm_squared_distances = max_of_l2_norm_squared_distance;
      binary_representation_file << max_of_l2_norm_squared_distance << ";" << endl;
    }
  binary_representation_file << "Max Eps: " << max_of_max_of_l2_norm_squared_distances << ";" << endl;
}


LaGenMatDouble
Problem_instance::A(int j) {
  int m = M();
  int m_j = M(j);
  int n = N();

  assert (m_j > 0);

  LaGenMatDouble A(m_j, n);

  for (int i = 0, ii = 0; i < m; i++, ii++) {
    if (assignment[i][j] == false) {
      ii--;
      continue;
    }
    for (int l = 0; l < n; l++)
      A(ii,l) = a[i](l);
  }
  return A;
}

LaGenMatDouble
Problem_instance::A_aff_reg(int j)
{
  LaGenMatDouble A = Problem_instance::A(j);
  int m_j = A.rows();
  int n = A.cols();	
  for (int i = 0; i < m_j; i++)
    A(i,n-1) = -1;
  return A;
}


LaGenMatDouble
Problem_instance::A_lin_reg(int j)
{
  LaGenMatDouble A = Problem_instance::A(j);
  int m = A.rows();
  int n = A.cols();	
  for (int i = 0; i < m; i++)
    A(i,n-1) = 0;
  // NOTE: since we use SVD to solve Ax = b --> x = pinv(A)*b, we always get the minimal norm solution.
  // Hence, even though, with the last column of zeros, gamma could have ANY value,
  // we will always get a NULL gamma. That's perfect for Panda
  return A;
}

LaVectorDouble
Problem_instance::b_reg(int j)
{
  assert (problem_type.problem != HCLUSTERING);

  int m = M();
  int m_j = M(j);
  int n = N();

  assert (m_j > 0);

  LaVectorDouble b(m_j,1);

  for (int i = 0, ii = 0; i < m; i++, ii++)	{
    if (assignment[i][j] == false ) {
      ii--;
      continue;
    }
    b(ii) = a[i](n-1);
  }
  return b;
}

LaVectorDouble
Problem_instance::c(int j)
{
  int m = M();
  int n = N();
  int m_j = M(j);
  LaVectorDouble c(n);
  for (int l = 0; l < n; l++)
    c(l) = 0;

  if (m_j == 0)
    return c;
  
  for (int i = 0; i < m; i++)
    {
      if (assignment[i][j] == false)
	    continue;
      for (int l = 0; l < n; l++)
	{
	  c(l) += a[i](l);
	}
    }
  for (int l = 0; l < n; l++)
    {
      c(l) /= m_j;
    }
  return c;
}

void
Problem_instance::update_distances(const int i, const int j) {

  // cout << "Problem_instance::update_distances" << endl;
  // cout << "w[j] " << w[j] << endl;
  // cout << "gamma[j] " << gamma[j] << endl;

  if (problem_type.problem == HCLUSTERING)
    distance[i][j] = Formulas::l2_norm_distance(a[i], w[j], gamma[j]);

  else if ((problem_type.problem == AFFINE_REGRESSION || problem_type.problem == LINEAR_REGRESSION || problem_type.problem == PAMF ))
    distance[i][j] = Formulas::residual(a[i], w[j], gamma[j]);

  else if (problem_type.problem == PCLUSTERING)
    distance[i][j] = Formulas::l2_norm_prototypal_distance(a[i], w[j]);

  else
    abort();
}

void
Problem_instance::update_distances(const int j)
{
  for (int i = 0; i < M(); i++)
    update_distances(i, j);
}

void
Problem_instance::update_distances()
{
  for (int j = 0; j < K(); j++)
    update_distances(j);
}

double
Problem_instance::point_measure(int i, int j) const
{
  return distance[i][j];
}
double
Problem_instance::cluster_measure(int j) const
{
  double measure = 0;
  for (int i = 0; i < M(); i++)
    // if (assignment[i][j] == true && )
    //   measure += pow(distance[i][j], 2);
    if (assignment[i][j] == true ) {
    //if ( !( problem_type.problem == AFFINE_REGRESSION || problem_type.problem == LINEAR_REGRESSION ) )
      if ( !( problem_type.problem == PAMF ) )
	measure += pow(distance[i][j], 2);
      else
	measure += distance[i][j];
    }
  return measure;
}

double
Problem_instance::solution_measure() const
{
  double measure = 0;
  for (int j = 0; j < K(); j++)
    measure += cluster_measure(j);
  return measure;
}

void
Problem_instance::get_min_max_dimensions(vector<double> & min_coords, vector<double> & max_coords) const
{
  //! variables
  int m = M();
  int n = N();

  for (int i = 0; i < m; i++) {
    for (int l = 0; l < n; l++) {
      if (a[i](l) < min_coords[l])
	min_coords[l] = a[i](l);
      
      if (a[i](l) > max_coords[l])
	max_coords[l] = a[i](l);
    }
  }
}

void
Problem_instance::normalize(double factor) {
  int m = M();
  int n = N();
  vector<double> min(m,100);
  vector<double> max(m,-100);
  get_min_max_dimensions(min, max);
  for (int i = 0; i < m; i++)
    for (int l = 0; l < n; l++)
      a[i](l) -= min[l];

  for (int i = 0; i < m; i++)
    for (int l = 0; l < n; l++)
      a[i](l) /= (max[l]-min[l]);

  if ( factor != 1 ) {
    for (int i = 0; i < m; i++)
      for (int l = 0; l < n; l++)
	a[i](l) *= factor;
  }
}

void
Problem_instance::cutoffprecision(double precision) {
  int m = M();
  int n = N();
  for (int i = 0; i < m; i++)
    for (int l = 0; l < n; l++)
      a[i](l) = floor(a[i](l) * precision) / precision;
}


void
Problem_instance::create_random_solution_and_clean_TS( base_generator_type & generator ) {

  this->randomly_populate_clusters( generator );
  this->update_distances();
  algorithm::plain_combinatorial_reassignment( *this );
  if ( this->ops.TS.flag_use_TS )
    this->tabu_list.empty();
  algorithm::plain_parameter_update( *this ); //if all my algs start by moving points, this is fine
  this->update_distances();

  //Later addition: this way, any randomly generated solution is RLP-feasible

  //if ( (this->problem_type.problem == AFFINE_REGRESSION || this->problem_type.problem == LINEAR_REGRESSION) ) {
  if ( this->ops.flag_generate_RLP_feasible_random_solutions ) {
    // pamf::reassign_with_RLP( *this );
    pamf::plain_RLP_parameter_update ( *this );
    pamf::plain_RLP_induced_combinatorial_reassignment( *this );
    algorithm::plain_parameter_update( *this );
    this->update_distances();
  }


}


bool
Problem_instance::assert_correctness () {
    
  int m = M();
  int k = K();


  //assert assignment correctness
  for (int i = 0; i < m; i++) {
    bool assigned = false;
    for (int j = 0; j < k; j++) {
      if ( assignment[i][j] == true )
	assigned = true;
      }
    if ( !assigned ) {
      cout << "assert_correctness = false --assignment failed" << endl;
      return false;
    }
  }

  //assert correctness of distances
  for (int i = 0; i < m; i++) {
    for (int j = 0; j < k; j++) {
      
      double actualDistance;

      switch (problem_type.problem ) {
      case HCLUSTERING:
	actualDistance = Formulas::l2_norm_distance(a[i], w[j], gamma[j]);
	break;
      case AFFINE_REGRESSION:
      case LINEAR_REGRESSION:
      case PAMF:
	actualDistance = Formulas::residual(a[i], w[j], gamma[j] ) ;
	break;
      case PCLUSTERING:
	actualDistance = Formulas::l2_norm_prototypal_distance(a[i], w[j]);
	break;
      }

      if ( distance[i][j] != actualDistance ) {
	  cout << "assert_correctness = false --distances failed: "
               << distance[i][j] << " != " << actualDistance
	       << endl;
	  cout << " i = " << i << endl;
	  cout << " j = " << i << endl;
	  cout << "a[i] = " << a[i] << endl;
	  cout << "w[j] = " << w[j] << endl;
	  cout << "gamma[j] = " << gamma[j] << endl;
	  cout << assignment << endl;
	  return false;
      }
    }
  }

  if ( problem_type.problem == PAMF ) {
    // double misclass_error = pamf::check_RLP( *this );
    double misclass_error = pamf::plain_RLP_parameter_update( *this );
    if ( misclass_error > 0.01 ) {
      cout << "Misclassification error = " << misclass_error << endl;
      return false;
    }
  }

  return true;

}



