#include "inclusions.h"

using namespace std;




Return_value
algorithm::PW(Problem_instance & pi, bool warmStar, double maxTime, base_generator_type & generator )
{
  const int k = pi.K();
  assert(k > 1);
  const int m = pi.M();

  Chronometer timer;
  timer.start();
  double init_time = timer.elapsed_time();

  double min_dist;
  int cl = -1;
  int old_cl = -1;
  int pt;			//! idx of the point chosen for reallocation

  vector<double> distances(m,-1);
  vector<bool> already_reassigned(m,false);
  vector<int> ordered_vector(m,-1);

  if ( warmStar == false)
    pi.create_random_solution_and_clean_TS ( generator );

  cout << "LogInitPW\t"
       << " initObj "      << pi.solution_measure()
       << endl;


  double last_obj = 1e300;
  double obj = pi.solution_measure();

  int iteration = 0;
  while ( obj < last_obj && timer.elapsed_time() < maxTime - TIMETOL ) {
    last_obj = obj;
    
    //! compiling the distances vector
    for (int i = 0; i < m; i++) {
      distances[i] = pi.distance[i][pi.J(i)];
      //distances[i] = random_generator::integer(5);
      ordered_vector[i] = i;
      already_reassigned[i] = false;
    }
    
    //! cycling over all the points
    for (int i = 0; i < m; i++) {
      //! find the maximum distance point
      pointwise_distance_order pdo;
      pdo.set_vectors(distances, already_reassigned);
      vector<int>::iterator mdpt = max_element(ordered_vector.begin(), ordered_vector.end(), pdo);
      pt = *mdpt;
      
      already_reassigned[pt] = true;
      
      min_dist = 1e300;
      old_cl = pi.J(pt);
      
      for (int j = 0; j < k; j++) {
	if (pi.distance[pt][j] < min_dist) {
	  min_dist = pi.distance[pt][j];
	  cl = j;
	}
      }
      
      cout << "i = " << i << " on cl = " << old_cl << " -- dist " << pi.distance[pt][old_cl];
      cout << "i = " << i << " to cl = " << cl << " -- dist " << pi.distance[pt][cl];

      cout << pi << endl;

      if (old_cl == cl)
	continue;
      else {
	cout << "[" << pt << "," << old_cl << "] --> [" << pt << "," << cl << "]" << endl;
	
	pi.assignment[pt][pi.J(pt)] = false;
	pi.assignment[pt][cl] = true;
	
	// exact::plain_recalculate_l2_norm_hyperplane_parameters(pi, old_cl);
	// exact::plain_recalculate_l2_norm_hyperplane_parameters(pi, cl);
	plain_parameter_update ( pi, old_cl );
	plain_parameter_update ( pi, cl );
	
	pi.update_distances();

	cout << "LogPointPW\t"
	     << " currObj "      << obj
	     << endl;

      cout << pi << endl;
      }
    }
    
    obj = pi.solution_measure();
    //cout << "It: " << iteration << "\t" << obj << endl;
    cout << "LogCurrPW\t"
	 << " currObj "      << obj
	 << endl;
    
    iteration++;
  } // end of the iterations
  
  cout << "LogFinalPW\t"
       << " currObj "      << obj
       << endl;


  Return_value retval;
  retval.obj = obj;
  retval.time_to_best = timer.elapsed_time() - init_time;
  retval.total_time = retval.time_to_best; //the same (last sol = best sol)
  retval.iterations_to_best = iteration;
  retval.total_iterations = retval.iterations_to_best;  //the same (last sol = best sol)

  return retval;
}

