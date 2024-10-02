#include "inclusions.h"

using namespace std;

//---------------------------------------------------------------------------------------
// ALGORITHMS
//---------------------------------------------------------------------------------------



Return_value
algorithm::papavero( Problem_instance & pi, bool warmStart, double maxTime, base_generator_type & generator ) {

  cout << "papavero " << endl;

  const int m = pi.M();
  const int k = pi.K();
  assert(k > 1);

  Chronometer timer;
  timer.start();

  Return_value retval;

  //initial random solution
  if ( warmStart == false ) {
    pi.create_random_solution_and_clean_TS ( generator );
  }


  Problem_instance best_pi(pi);

  retval.obj = pi.solution_measure();


  cout << "LogInitPAPAVERO\t"
      << " initObj "      << pi.solution_measure()
      << endl;

  int i; //this way, I can use it outside
  for ( i = 0; timer.elapsed_time() <= maxTime - TIMETOL ; i++ ) {

    if (pi.ops.TS.flag_use_TS)  //this way, its starts from 0, rather than from -1
      pi.tabu_list++;

    
    // plain_combinatorial_reassignment( pi );
    // papavero_combinatorial_reassignment( pi );
    papavero_combinatorial_reassignment( pi, generator );
    
    plain_parameter_update( pi );
    pi.update_distances();

    // pamf::perform_RLP( pi, true );
    pamf::plain_RLP_parameter_update ( pi );


    if ( pi.ops.flag_intertwine_RLP ) {
      cout << "Intertwining RLP" << endl;
      //pamf::reassign_with_RLP( pi );
      pamf::plain_RLP_induced_combinatorial_reassignment( pi );

      plain_parameter_update( pi );
      pi.update_distances();

#ifdef __DEBUG__
      double misclass_error = pamf::check_RLP ( pi );

      if ( pi.N() == 2 && misclass_error > 0 ) {
      	cout << "Assertion of separability failed: " << misclass_error << endl;
      	string matlab_data_file_name("pino");
      	matlab_data_file_name += string("_sol.m");
      	cout << "Outputting to matlab on file " << matlab_data_file_name << endl;
      	pi.output_solution_to_matlab(matlab_data_file_name, true); //last parameter (true) means: print cluster lines = {on, off}
      	cout << pi << endl;
      }	

      assert ( misclass_error == 0 );
#endif
    }

    double curr_obj = pi.solution_measure();

    if ( curr_obj < best_pi.solution_measure() ) {
      best_pi = pi;
      retval.obj = curr_obj;
    }

    cout << "LogCurrPAPAVERO\t"
        << " bestObj "      << retval.obj
        << " currObj "      << curr_obj
        << " timeLeft "     << maxTime - timer.elapsed_time()
        << endl;


    retval.total_iterations = i;

  } // end of iterations loop

  //update pi to best_pi (it is used in APR) -- also, if a better sol is found via polishing, it must be stored
  if ( pi.solution_measure() < best_pi.solution_measure() )
    best_pi = pi;
  else
    pi = best_pi; //when terminating, pi must be best sol that we return

  //applied to the BEST SOLUTION FOUND, before returning
  if ( pi.ops.flag_RLP_at_the_end ) {
    cout << "RLP was not used in intertwining not directly in the combinatorial reassignment. Reassigning with RLP now ON THE BEST SOLUTION FOUND, then updating the parameters and the distances" << endl;
    // cout << "Solution of value " << pi.solution_measure() << " before calling pamf::reassign_with_RLP( pi ); " << pi.assignment << endl;
 
    //pamf::reassign_with_RLP( pi );
    pamf::plain_RLP_parameter_update ( pi );
    pamf::plain_RLP_induced_combinatorial_reassignment( pi );

    // cout << "Solution of value " << pi.solution_measure() << " after calling pamf::reassign_with_RLP( pi ); " << pi.assignment << endl;
    algorithm::plain_parameter_update( pi );
    pi.update_distances();
    retval.obj = pi.solution_measure();
    cout << "LogCurrPAPAVERO\t"
        << " bestObj "      << retval.obj
	<< " currObj "      << pi.solution_measure()
        << " timeLeft "     << maxTime - timer.elapsed_time()
        << endl;
  }
  

  //iterations_to_best and time_to_best are NO MORE SUPPORTED, due to S.I.M.P.L.I.F.I.C.A.T.I.O.N.
  retval.obj = pi.solution_measure();
  retval.time_to_best = timer.elapsed_time();
  retval.total_time = retval.time_to_best; //the same (last sol = best sol)
  //  retval.iterations_to_best = i + BM_on_pi.total_iterations + BM_on_best_pi.total_iterations;
  retval.iterations_to_best = i;
  retval.total_iterations = retval.iterations_to_best;  //the same (last sol = best sol)


  cout << "LogFinalPAPAVERO\t"
      // << " bestObj "      << best_pi.solution_measure()
      << " bestobj "      << pi.solution_measure()
      << " timeLeft "     << maxTime - timer.elapsed_time()
      << endl;

  return retval;
}//algorithm::PAPAVERO()



Return_value
algorithm::BM(Problem_instance & pi, bool warmSTart, double maxTime, base_generator_type & generator ) {

  const int k = pi.K();
  assert(k > 1);

  Chronometer timer;
  timer.start();

  if ( warmSTart == false ) {
    pi.create_random_solution_and_clean_TS( generator );
  }


  cout << "LogInitBM\t"
      << " initObj "      << pi.solution_measure()
      << endl;

  double last_obj = 1e300;
  double obj = pi.solution_measure();

  int i;
  //used in ejor2011
  for ( i = 0; obj < last_obj && timer.elapsed_time() <= maxTime - TIMETOL; i++ ) {
  //for ( i = 0; obj != last_obj && timer.elapsed_time() <= maxTime - TIMETOL; i++ ) {

    if (pi.ops.TS.flag_use_TS)
      pi.tabu_list++;

    last_obj = pi.solution_measure();

    // if ( pi.ops.flag_intertwine_RLP == false ) { // if true, they are done in the random solution (first sol) and before exiting the cycle (any sol)
    //   plain_parameter_update( pi );
    //   pi.update_distances();
    // }

    if ( pi.ops.flag_RLP_based_combinatorial_reassignment )
      pamf::RLP_based_combinatorial_reassignment( pi, maxTime - timer.elapsed_time() );
    else
      plain_combinatorial_reassignment( pi );


    if ( pi.ops.flag_intertwine_RLP ) {
      cout << "LogCurrBM\t Intertwining RLP" << endl;
      //pamf::reassign_with_RLP( pi );
      pamf::plain_RLP_parameter_update ( pi );
      pamf::plain_RLP_induced_combinatorial_reassignment( pi );

      // plain_parameter_update( pi );
      // pi.update_distances();
    }

    plain_parameter_update( pi );
    pi.update_distances();

    // cout << "sol measure " << pi.solution_measure() << endl;

    if ( pi.ops.bemporad_delta < 0 ) { //apply criterion
      cout << "Bemporad kicking in as a step!!" << endl;
      bemporads_criterion( pi );
    }

    obj = pi.solution_measure();

    cout << "LogCurrBM\t"
        << " currObj "      << obj
        << endl;

  } // end of the iterations

  if ( pi.ops.bemporad_delta > 0 ) { //apply criterion
    cout << "Bemporad kicking in at the end!" << endl;
    bemporads_criterion( pi );
    obj = pi.solution_measure();
  }


  if ( pi.ops.flag_RLP_at_the_end ) {
    cout << "RLP was not used in intertwining not directly in the combinatorial reassignment. Reassigning with RLP now, then updating the parameters and the distances" << endl;
    //cout << "Assignment before calling pamf::reassign_with_RLP( pi ); " << pi.assignment << endl;
    //pamf::reassign_with_RLP( pi );
    pamf::plain_RLP_parameter_update ( pi );
    pamf::plain_RLP_induced_combinatorial_reassignment( pi );

    //cout << "Assignment after calling pamf::reassign_with_RLP( pi ); " << pi.assignment << endl;
    //algorithm::plain_parameter_update( pi ); //removed on 9.oct.2014: this way, when comparing to the "two phase" algorithm, it does NOT perform a "Domain partitioning" step which also recomputes the submodels; rather, it just reassigns the points to the submodel corresponding to the subdomain containing them and returns; this way, WORSE objective function values should be found!
    pi.update_distances();
    obj = pi.solution_measure();
    cout << "LogCurrBM\t"
        << " currObj "      << obj
        << endl;
  }
 
  cout << "LogFinalBM\t"
      << " currObj "      << obj
      << endl;


  Return_value retval;

  retval.obj = obj;
  retval.time_to_best = timer.elapsed_time();
  retval.total_time = retval.time_to_best; //the same (last sol = best sol)
  retval.iterations_to_best = i;
  retval.total_iterations = retval.iterations_to_best;  //the same (last sol = best sol)
  

  return retval;
}

Return_value
algorithm::TSBM(Problem_instance & pi, bool warmStar, double maxTime, base_generator_type & generator ) {
  
  cout << "Return_value algorithm::TSBM(Problem_instance & pi, bool warmStar, double maxTime, base_generator_type & generator )" << endl;

  const int m = pi.M();
  const int n = pi.N();
  const int k = pi.K();
  assert(k > 1);

  Chronometer timer;
  timer.start();

  if ( warmStar == false)
    pi.create_random_solution_and_clean_TS ( generator );

  cout << "LogInitTSBM\t"
      << " initObj "      << pi.solution_measure()
      << endl;

  Problem_instance best_pi = pi;

  cout << "TSBM:"
      << " conservative = " << pi.ops.APR.conservative
      << " samples = "<< pi.ops.TS.neighborhood_samples
      << " swaps = " << pi.ops.TS.swaps
      << endl;


  int iteration = 0;
  while ( timer.elapsed_time() < maxTime - TIMETOL ) {

    if (pi.ops.TS.flag_use_TS)
      pi.tabu_list++;

    // pi.tabu_list.full_print_tabu_list();

    //copy of the matrix of distances before performing reassignments; used to put stuff in the tabu list
    vector< double > distances_before( m, -1 );
    vector< int > assignments_before( m, -1 );
    for ( int i = 0; i < m; i++ ) {
      int curr_j = pi.J( i );
      distances_before[i] = pi.distance[i][curr_j];
      assignments_before[i] = pi.J( i );
    }

    //begin copied from BM

    if ( pi.ops.flag_RLP_based_combinatorial_reassignment )
      pamf::RLP_based_combinatorial_reassignment( pi, maxTime - timer.elapsed_time() );
    else
      plain_combinatorial_reassignment( pi );


    if ( pi.ops.flag_intertwine_RLP ) {
      cout << "Intertwining RLP" << endl;
      //pamf::reassign_with_RLP( pi );
      pamf::plain_RLP_parameter_update ( pi );
      pamf::plain_RLP_induced_combinatorial_reassignment( pi );

      // plain_parameter_update( pi );
      // pi.update_distances();
    }

    plain_parameter_update( pi );
    pi.update_distances();

    //end copied from BM

    bool worsening = false;
    if ( pi.solution_measure() < best_pi.solution_measure() ) { //store best
      best_pi = pi;
    } else {
      worsening = true;
      //add the whole new solution to the tabu list
      for ( int i = 0; i < m; i++ ) {
	int old_j = assignments_before[i];
	int new_j = pi.J( i );
	// pi.tabu_list.add_old_and_remove_new( i, old_j, distances_before[i], new_j );
	if ( pi.tabu_list.use_aspiration_criterion() == false )
	  pi.tabu_list.add( i, old_j );
	else
	  pi.tabu_list.add( i, old_j, distances_before[i] );
	
	// cout << "adding pair " << i << " " << old_j << " of distance " << distances_before[i] << endl;

      }
    }

    iteration++;

    //cout << "It: " << iteration << "\t" << obj << endl;
    cout << "LogCurrTSBM\t"
        << " currObj "      << pi.solution_measure()
        << " worsening "    << worsening
        << " bestObj "      << best_pi.solution_measure()
        << " iterations "   << iteration
        << " timeLeft "     << maxTime - timer.elapsed_time()
        << endl;

    // pi.tabu_list.full_print_tabu_list();

  } // end of the iterations

  cout << "LogFinalTSBM\t"
      // << " currObj "      << obj
      << " currObj "      << pi.solution_measure()
      << " bestObj "      << best_pi.solution_measure()
      << " iterations "   << iteration
      << endl;


  Return_value retval;
  retval.obj = best_pi.solution_measure();
  retval.time_to_best = timer.elapsed_time();
  retval.total_time = retval.time_to_best; //the same (last sol = best sol)
  retval.iterations_to_best = iteration;
  retval.total_iterations = retval.iterations_to_best;  //the same (last sol = best sol)

  return retval;
}



// Return_value
// algorithm::BM(Problem_instance & pi, bool warmSTart, double maxTime, base_generator_type & generator ) {

//   const int k = pi.K();
//   assert(k > 1);

//   Chronometer timer;
//   timer.start();

//   if ( warmSTart == false ) {
//     pi.create_random_solution_and_clean_TS( generator );
//   }


//   cout << "LogInitBM\t"
//       << " initObj "      << pi.solution_measure()
//       << endl;

//   double last_obj = 1e300;
//   double obj = pi.solution_measure();

//   int i;
//   //used in ejor2011
//   //r ( i = 0; obj < last_obj && timer.elapsed_time() <= maxTime - TIMETOL; i++ ) {
//   for ( i = 0; obj != last_obj && timer.elapsed_time() <= maxTime - TIMETOL; i++ ) {

//     if (pi.ops.TS.flag_use_TS)
//       pi.tabu_list++;

//     last_obj = pi.solution_measure();

//     //cout << pi << endl;

//     plain_parameter_update( pi );

//     pi.update_distances();

//     plain_combinatorial_reassignment( pi );


//     cout << "pi.ops.flag_intertwine_RLP = " << pi.ops.flag_intertwine_RLP << endl;

//     if ( pi.ops.flag_intertwine_RLP && (pi.problem_type.problem == AFFINE_REGRESSION || pi.problem_type.problem == LINEAR_REGRESSION) ) {
//       // if ( i > 1 )
//       cout << "Before       pamf::reassign_with_RLP( pi ); " << endl;
//       pamf::reassign_with_RLP( pi );
//       cout << "After       pamf::reassign_with_RLP( pi ); " << endl;
//       cout << pi.assignment << endl;
//       plain_parameter_update( pi );
//       pi.update_distances();
//     }

//     obj = pi.solution_measure();

//     cout << "LogCurrBM\t"
//         << " currObj "      << obj
//         << endl;

//   } // end of the iterations
 
//   cout << "Pre log final" << endl;

//   cout << pi.assignment << endl;

//   cout << "LogFinalBM\t"
//       << " currObj "      << obj
//       << endl;


//   Return_value retval;

//   retval.obj = obj;
//   retval.time_to_best = timer.elapsed_time();
//   retval.total_time = retval.time_to_best; //the same (last sol = best sol)
//   retval.iterations_to_best = i;
//   retval.total_iterations = retval.iterations_to_best;  //the same (last sol = best sol)

//   cout << "pre return" << endl;
//   cout << pi.assignment << endl;


//   return retval;
// }

Return_value
algorithm::PR( Problem_instance & pi, double alpha, bool warmStart, double maxTime, base_generator_type & generator ) {

  const int m = pi.M();
  const int k = pi.K();
  assert(k > 1);

  Chronometer timer;
  timer.start();

  double init_alpha = alpha; //needed for linear decrease

  vector<bool> dummy( m, false );
  vector<vector<bool> > ill_ass_vector( k, dummy );

  Return_value retval;

  double ill_assigned_num = 1e300;

  //initial random solution
  if ( warmStart == false ) {
    pi.create_random_solution_and_clean_TS ( generator );
  }

  // //We assume that the starting solution --must be-- RLP-feasible. This is needed, since PR might find worsened solutions and must keep track of the best one that is found. If we allow to evaluate any solution which is infeasbile, the algorith will likely find it better than any feasible one and return it.
  // if ( (pi.problem_type.problem == AFFINE_REGRESSION || pi.problem_type.problem == LINEAR_REGRESSION) ) {
  //   cout << "Before       pamf::reassign_with_RLP( pi ); " << endl;
  //   pamf::reassign_with_RLP( pi );
  //   cout << "After       pamf::reassign_with_RLP( pi ); " << endl;
  //   cout << pi.assignment << endl;
  //   plain_parameter_update( pi );
  //   pi.update_distances();
  // }

  Problem_instance best_pi(pi);

  retval.obj = pi.solution_measure();


  cout << "LogInitPR\t"
      << " initObj "      << pi.solution_measure()
      << " alpha "        << alpha
      << endl;

  int i; //this way, I can use it outside
  for ( i = 0; ( timer.elapsed_time() <= maxTime - TIMETOL ) &&
  ( alpha > 0 )                                      &&
  ( ill_assigned_num > 0 ) ; i++ ) {

    if (pi.ops.TS.flag_use_TS)  //this way, its starts from 0, rather than from -1
      pi.tabu_list++;

    //find ill-assigned points
    for (int j = 0; j < k; j++)
      ill_ass_vector[j] = find_ill_assigned( pi, j, alpha, generator );

    //merge "local" (per cluster) ill_ass_vectors into a global one
    vector<bool>global_ill_ass(m, false);
    for (int ii = 0; ii < m; ii++)
      for (int j = 0; j < k; j++)
        if ( ill_ass_vector[j][ii] == true )
          global_ill_ass[ii] = true;

    //reassign ill-assigned
    ill_assigned_num = 0;

    // for (int j = 0; j < k; j++) {
    //   degenerancy_avoider( pi, j, generator );
    // }

    if ( pi.ops.flag_RLP_based_combinatorial_reassignment ) { //math heuristic, using CPLEX
      pamf::RLP_based_combinatorial_reassignment( pi, maxTime - timer.elapsed_time(), &global_ill_ass );
    }
    else {

      // cout << "APR combinatorial reassignment" << endl;
      for (int j = 0; j < k; j++) {
	ill_assigned_num += APR_combinatorial_reassignment(pi, j, ill_ass_vector[j], generator );
	// cout << "after APR_comb_reass on cluster " << j << " " << pi.assignment << endl;
	
	if ( pi.ops.APR.extra_update == true ) {
	  plain_parameter_update( pi );
	  pi.update_distances();
	}

      }

      // ill_assigned = full_APR_combinatorial_reassignment(pi, global_ill_ass, generator );
      //for some reason, it was giving different results...
      
      //reassign nonill-assigned
      // cout << "Plain combinatorial reassignment" << endl;

      plain_combinatorial_reassignment(pi, &global_ill_ass);

    }


    // //BEMPORAD: begin
    // // loop over points
    // // loop of clusters
    // // mark all clusters with distance within DELTA
    // // if > 1:
    // // mark all c = 10 points closest to the current point
    // // for each point, read its cluster
    // // count the popularity of the clusters
    // // reassign point to most popular
    // if ( pi.ops.bemporad_delta != 0 ) { //apply criterion
    //   cout << "Bemporad kicking in!" << endl;
    //   double delta;
    //   if ( pi.ops.bemporad_delta < 0 ){
    // 	delta = pi.solution_measure()/m;
    //   }
    //   else {
    // 	delta = pi.ops.bemporad_delta;
    //   }
    // // if ( 1 ) { //apply criterion
    // //   double delta;
    // //   delta = pi.solution_measure()/m;
    //   int c = 10;
    //   for (int ii = 0; ii < m; ii++) {
    // 	bool already_found = false;
    // 	for (int j = 0; j < k; j++) {
    // 	  if ( pi.distance[ii][j] < delta && already_found == true) {
    // 	    //undecidable point found!
    // 	    // cout << "point " << ii << ", assigned to cluster " << pi.J(ii) << " is undecidable w.r.t. cluster " << j << endl;
    // 	    vector <int> votes (k, 0);
    // 	    vector <bool> checked (m, false);	  
    // 	    for (int passes = 0; passes < c; passes++) {
    // 	      // cout << "pass " << passes << endl;
    // 	      double min_dist = 1e300;
    // 	      int min_dist_point = -1;
    // 	      for (int iii = 0; iii < m; iii++) {
    // 		// cout << "Checking distance between point " << ii << " and point " << iii << endl;
    // 		if (iii == ii) {
    // 		  continue;
    // 		}
    // 		double dist = Formulas::l2_norm_prototypal_distance(pi.a[ii], pi.a[iii]);
    // 		// cout << "dist " << dist << endl;
    // 		if ( dist < min_dist && checked[iii] == false ) {
    // 		  min_dist = dist;
    // 		  min_dist_point = iii;
    // 		}
    // 	      }
    // 	      votes[pi.J(min_dist_point)]++;
    // 	      checked[min_dist_point] = true; 
    // 	    } //end of passes
    // 	    int most_voted_cluster = -1;
    // 	    int max_votes = 0;
    // 	    for (int jj = 0; jj < k; jj++) {
    // 	      if ( votes[jj] > max_votes ) {
    // 		max_votes = votes[jj];
    // 		most_voted_cluster = jj;
    // 	      }
    // 	    }
    // 	    pi.assignment[ii][pi.J(ii)] = false;
    // 	    pi.assignment[ii][most_voted_cluster] = true;
    // 	    cout << "Bemporad: reassigned point " << ii << " from cluster " << pi.J(ii) << " to cluster " << most_voted_cluster << endl;
    // 	  }
    // 	  else if ( pi.distance[ii][j] < delta && already_found == false) {
    // 	    already_found = true;
    // 	  }
    // 	}
    //   }
    //   plain_parameter_update( pi );
    //   pi.update_distances();
    //   //BEMPORAD: end
    // }
    if ( pi.ops.bemporad_delta < 0 ) { //apply criterion
      cout << "Bemporad kicking in as a step!!" << endl;
      bemporads_criterion( pi );
    }
    

    if ( pi.ops.flag_intertwine_RLP ) {
      cout << "Intertwining RLP" << endl;
      //pamf::reassign_with_RLP( pi );
      pamf::plain_RLP_parameter_update ( pi );
      pamf::plain_RLP_induced_combinatorial_reassignment( pi );

      // plain_parameter_update( pi );
      // pi.update_distances();

#ifdef __DEBUG__
      double misclass_error = pamf::check_RLP ( pi );

      if ( pi.N() == 2 && misclass_error > 0 ) {
      	cout << "Assertion of separability failed: " << misclass_error << endl;
      	string matlab_data_file_name("pino");
      	matlab_data_file_name += string("_sol.m");
      	cout << "Outputting to matlab on file " << matlab_data_file_name << endl;
      	pi.output_solution_to_matlab(matlab_data_file_name, true); //last parameter (true) means: print cluster lines = {on, off}
      	cout << pi << endl;
      }	

      assert ( misclass_error == 0 );
#endif
    }

    //update parameters
    // for (int j = 0; j < k; j gg++)
    //   plain_parameter_update(pi, j );
    // cout << "Plain parameter update + distances" << endl;
    plain_parameter_update( pi );
    pi.update_distances();


    double curr_obj = pi.solution_measure();

    if ( curr_obj < best_pi.solution_measure() ) {
      best_pi = pi;
      retval.obj = curr_obj;
    }

    cout << "LogCurrPR\t"
        << " bestObj "      << retval.obj
        << " currObj "      << curr_obj
        << " alpha "        << alpha
        << " ill_assigned " << ill_assigned_num
        << " timeLeft "     << maxTime - timer.elapsed_time()
        << endl;

    // decrease alpha
    if (i % pi.ops.APR.alpha_decrease_frequency == 0) {
      if (pi.ops.APR.flag_use_linear_threshold == false)
        alpha = alpha * pi.ops.APR.rho_down;
      else
        //alpha = ((pi.ops.APR.init_alpha - i*pi.ops.APR.rho_down) >= 0) ? (pi.ops.APR.init_alpha - i*pi.ops.APR.rho_down) : 0;
        alpha = ((init_alpha - i*pi.ops.APR.rho_down) >= 0) ? (init_alpha - i*pi.ops.APR.rho_down) : 0;

      //! zero-threshold alpha
      if (alpha < pi.ops.APR.min_alpha)
        alpha = 0;
    }

    retval.total_iterations = i;

  } // end of iterations loop


  //update pi to best_pi (it is used in APR) -- also, if a better sol is found via polishing, it must be stored
  if ( pi.solution_measure() < best_pi.solution_measure() )
    best_pi = pi;
  else
    pi = best_pi; //when terminating, pi must be best sol that we return

  if ( pi.ops.bemporad_delta > 0 ) { //apply criterion
    cout << "Bemporad kicking in at the end!" << endl;
    bemporads_criterion( pi );
    // obj = pi.solution_measure();
  }

  //applied to the BEST SOLUTION FOUND, before returning
  if ( pi.ops.flag_RLP_at_the_end ) {
    cout << "RLP was not used in intertwining not directly in the combinatorial reassignment. Reassigning with RLP now ON THE BEST SOLUTION FOUND, then updating the parameters and the distances" << endl;
    // cout << "Solution of value " << pi.solution_measure() << " before calling pamf::reassign_with_RLP( pi ); " << pi.assignment << endl;
    // pamf::reassign_with_RLP( pi );
    pamf::plain_RLP_parameter_update ( pi );
    pamf::plain_RLP_induced_combinatorial_reassignment( pi );

    // cout << "Solution of value " << pi.solution_measure() << " after calling pamf::reassign_with_RLP( pi ); " << pi.assignment << endl;
    //algorithm::plain_parameter_update( pi ); //removed on 9.oct.2014: this way, we can compare to the k-HC heuristic with domain partitioning at the end, so to show the impact of the idea of doing the domain partitioning step at each iteration EVEN against an algorithm using our criterion, but which carries out such step only once at the end
    pi.update_distances();
    retval.obj = pi.solution_measure();
    cout << "LogCurrPR\t"
        << " bestObj "      << retval.obj
	<< " currObj "      << pi.solution_measure()
        << " alpha "        << alpha
        << " ill_assigned " << ill_assigned_num
        << " timeLeft "     << maxTime - timer.elapsed_time()
        << endl;
  }
  

  //iterations_to_best and time_to_best are NO MORE SUPPORTED, due to S.I.M.P.L.I.F.I.C.A.T.I.O.N.
  retval.obj = pi.solution_measure();
  retval.time_to_best = timer.elapsed_time();
  retval.total_time = retval.time_to_best; //the same (last sol = best sol)
  //  retval.iterations_to_best = i + BM_on_pi.total_iterations + BM_on_best_pi.total_iterations;
  retval.iterations_to_best = i;
  retval.total_iterations = retval.iterations_to_best;  //the same (last sol = best sol)


  cout << "LogFinalPR\t"
      // << " bestObj "      << best_pi.solution_measure()
      << " bestobj "      << pi.solution_measure()
      << " alpha "        << alpha
      << " ill_assigned " << ill_assigned_num
      << " timeLeft "     << maxTime - timer.elapsed_time()
      << endl;

  return retval;
}//algorithm::PR()


Return_value
algorithm::anytimer( Problem_instance & pi, t_algorithm_choice algorithm, double maxTime, base_generator_type & generator ) {

  Chronometer timer;
  timer.start();

//  pi.create_random_solution_and_clean_TS( generator );

  Return_value ret_global;
  Return_value ret_local;

  Problem_instance best_pi( pi );

  ret_global.obj = pi.solution_measure();
  ret_global.time_to_best = timer.elapsed_time();
  ret_global.iterations_to_best = 0;

  cout << "LogInitANY\t"
       << " initObj "      << pi.solution_measure()
       << endl;

  uint i;
  //for ( i = 0; /**/; i++ ) {
  for ( i = 0; timer.elapsed_time() <= maxTime - TIMETOL; i++ ) {
    //for ( i = 0; timer.elapsed_time() <= maxTime ; i++ ) {

    cout << "LogCurrANY\t" 
        << " start "          << i+1
        << " remaining time " << maxTime - timer.elapsed_time()
        << endl;

    //base_generator_type generator( i );
    //Given generator must be used, otherwise this algorithm could not be multistarted (needed for statistical comparison purposes)

    pi.create_random_solution_and_clean_TS( generator );

    switch ( algorithm ) {
    case t_aBM:
      ret_local = algorithm::BM( pi, false,
          std::max( maxTime - timer.elapsed_time(), 0.0 ),
          generator );
      break;
    case t_aPR:
      ret_local = algorithm::PR( pi, pi.ops.APR.init_alpha, false,
          std::max( maxTime - timer.elapsed_time(), 0.0 ),
          generator );
      break;
    case t_aPW:
      ret_local = algorithm::PW( pi, false,
          std::max( maxTime - timer.elapsed_time(), 0.0 ),
          generator );
      break;
    case t_aPAPAVERO:
      ret_local = algorithm::papavero( pi, false,
          std::max( maxTime - timer.elapsed_time(), 0.0 ),
          generator );
      break;
    }    

    // if ( timer.elapsed_time() >= maxTime ) {
    //   cout << "Anytime: timeout" << endl;
    //   break;
    // }

    if ( ret_local.obj < ret_global.obj ) {
      best_pi = pi;
      ret_global.obj = ret_local.obj;
      ret_global.time_to_best = ret_global.total_time + ret_local.time_to_best;
      ret_global.iterations_to_best = ret_global.total_iterations + ret_local.iterations_to_best;
    }

    ret_global.total_iterations += ret_local.total_iterations;
    ret_global.total_time = timer.elapsed_time();

    cout << "LogCurrANY\t"
        << " bestObj "      << ret_global.obj
        << " currObj "      << ret_local.obj
        << " totIters "     << ret_global.total_iterations+1
        << " bestIters "    << ret_global.iterations_to_best+1
        << " time "         << ret_global.total_time
        << " bestTime "     << ret_global.time_to_best
        << " starts "       << i+1
        << endl;

  } //end of multi-starts



  cout << "LogGlobal\t"
      << " bestObj "      << ret_global.obj
      << " totIters "     << ret_global.total_iterations+1
      << " bestIters "    << ret_global.iterations_to_best+1
      << " time "         << ret_global.total_time
      << " bestTime "     << ret_global.time_to_best
      << " starts "       << i+1
      << endl;

  pi = best_pi;

  ret_global.total_starts = i+1;

  return ret_global;
}

//Return
Problem_instance
algorithm::multistarter( Problem_instance & pi, t_algorithm_choice algorithm, double maxTime, base_generator_type & generator ) {

  cout << "algorithm::multistarter - " << algorithm << endl;

  Return_value ret;

  Problem_instance best_pi;

  //variables for statistics
  double best_obj = 1e300;
  double best_occurr = 0;
  vector<double> obj_per_start;
  vector<double> iter_per_start;
  vector<double> time_per_start;
  vector<double> start_per_start;

  for ( uint i = 0; i < pi.ops.multistarts; i++ ) {

    cout << "--------------------Start #" << i << endl;

    base_generator_type generator( i ); //each start uses the same generator (recreated): this way, all the algoritms that are multistarted start from the same random solution

    pi.create_random_solution_and_clean_TS( generator );

    cout << "LogInitMS\t"
        << " initObj "      << pi.solution_measure()
        << endl;


    switch ( algorithm ) {
    case t_mBM:
      ret = algorithm::BM( pi, true, 1e300, generator );
      break;
    case t_mPR:
      ret = algorithm::PR( pi, pi.ops.APR.init_alpha, true, 1e300, generator );
      break;
    case t_mPAPAVERO:
      ret = algorithm::papavero( pi, true, 1e300, generator );
      break;
    case t_maBM:
      ret = anytimer( pi, t_aBM, maxTime, generator );
      break;
    case t_maPR:
      ret = anytimer( pi, t_aPR, maxTime, generator );
      break;
    case t_maPAPAVERO:
      ret = anytimer( pi, t_aPAPAVERO, maxTime, generator );
      break;      
    case t_mAPR:
      ret = algorithm::APR( pi, true, maxTime, generator );
      break;
    case t_mPW:
      ret = algorithm::PW( pi, true, maxTime, generator );
      break;
    case t_maPW:
      ret = anytimer( pi, t_aPW, maxTime, generator );
      break;
    case t_msBM:
      ret = algorithm::anytimeShakenBM( pi, maxTime, generator );
      break;
    case t_mTS:
      ret = algorithm::TS( pi, false, maxTime, generator );
      break;
    case t_mTSBM:
      ret = algorithm::TSBM( pi, false, maxTime, generator );
      break;
    case t_mTSpamf:
      ret = algorithm::TSpamf( pi, true, maxTime, generator );
      break;
    case t_mJustRandom:
      ret = algorithm::justRandom( pi, true, maxTime, generator );
      break;
    default:
      cout << "t_algorithm incompatible with algorithm::multistarter" << endl;
      abort();
    }

    obj_per_start.push_back( ret.obj );
    time_per_start.push_back( ret.total_time );
    iter_per_start.push_back( ret.total_iterations );
    start_per_start.push_back ( ret.total_starts );

    if (ret.obj < best_obj) {
      best_obj = ret.obj;
      best_pi = pi;
      best_occurr = 1;
    }
    else if (best_obj - 1e-06 <= ret.obj && ret.obj <= best_obj + 1e-06 )
      best_occurr ++;

    cout << "LogCurrMS\t"
	 << " obj "          << ret.obj
	 << " starts "       << i+1
	 << " iters "        << ret.total_iterations+1
	 << " time "         << ret.total_time
	 << endl;

    cout << "LogGlobalMS\t"
	 << " bestObj "        << best_obj
	 << " starts "         << i+1
	 << " avgObj "         << utility::average( obj_per_start )
	 << " devStdObj "      << utility::dev_std( obj_per_start )
	 << " minObj "         << utility::min( obj_per_start )
	 << " maxObj "         << utility::max( obj_per_start )
	 << " avgIters "       << utility::average( iter_per_start )
	 << " devStdIter "     << utility::dev_std( iter_per_start )
	 << " minIters "       << utility::min( iter_per_start )
	 << " maxIter "        << utility::max( iter_per_start )
	 << " avgTime "        << utility::average( time_per_start )
	 << " devStdTime "     << utility::dev_std( time_per_start )
	 << " avgStarts "      << utility::average( start_per_start )
	 << " devStdStarts "   << utility::dev_std( start_per_start )
	 << " minStarts "      << utility::min( start_per_start )
	 << " maxStarts "      << utility::max( start_per_start )
	 << " bestFreq "       << best_occurr / obj_per_start.size()
	 << endl;

    cout << "LogObjGlobalMS\t"
         << obj_per_start << endl;

    cout << "LogIterGlobalMS\t"
         << iter_per_start << endl;

    cout << "LogStartsGlobalMS\t"
         << start_per_start << endl;

  } //end of multi-starts

  pi = best_pi;  //to return the best sol: not used in EJOR2011 code

  return best_pi;
} //algorithm::multistart()




Problem_instance
algorithm::pts_neighboring_solution( const Problem_instance & pi, vector<double> & probs, base_generator_type & generator ) {

  int m = pi.M();

  Problem_instance new_pi = pi;

  vector<int> is;
  vector<int> oldjs;
  vector<int> newjs;
  vector<double> ds;

  vector<bool>ill_ass(m, false);

  for (int g = 0; g < pi.ops.TS.swaps; g++) {

    int i, new_j;

    i = random_generator::sample_discrete_distribution( probs, generator );
    int old_j = new_pi.J( i );

    if ( pi.ops.APR.conservative == false )
      new_j = suggest_nontabu_reassignment( new_pi, i, generator );
    else
      new_j = suggest_conservative_nontabu_reassignment( new_pi, i, generator );

    // cout << " i = " << i << " --> " << " j = " << new_j << endl; //PATANJALI

    if ( new_j == -1)
      continue;
    else
      ill_ass[i] = true; //needed if we perform a reassignment of non ill-assigned points afterwards

    new_pi.assignment[i][old_j] = false;
    new_pi.assignment[i][new_j] = true;

    if ( pi.ops.TS.flag_use_TS ) {    // 2Feb2012: skip ahead if not TS (regular one, not the TSpamf one) is used
      is.push_back( i );
      oldjs.push_back( old_j );
      newjs.push_back(new_j);
      ds.push_back( pi.distance[i][old_j] );
    }

  } //all swaps performed: a neighboring solution has been created

  // //DEBATABLE!!!
  // if ( pi.ops.APR.conservative == true )
  //   plain_combinatorial_reassignment(pi, &ill_ass);


  plain_parameter_update( new_pi ); //many points are swapped: must update everything

  new_pi.update_distances( );

  cout << "New neighboring solution value = " << new_pi.solution_measure() << endl;

  // if ( pi.ops.TS.cooper )
  if ( pi.ops.TS.cooper == 1 )
    algorithm::BM( new_pi, true, 1e300, generator ); //issue: the moves done in cooper's method are not added to the tabu list...

  if ( pi.ops.TS.flag_use_TS ) {   //2Feb2012: made faster: if no TS is used, just skip ahead
    if ( new_pi.solution_measure() > pi.solution_measure() ) { //worsening solution
      // cout << "Worsening sol. found: adding stuff to tabu list"  << endl;
      for ( int aaa = 0; aaa < is.size(); aaa++ )
	new_pi.tabu_list.add_old_and_remove_new( is[aaa], oldjs[aaa], ds[aaa], newjs[aaa] );
    }
  }

  return new_pi;
}


Return_value
algorithm::TS(Problem_instance & pi, bool warmStar, double maxTime, base_generator_type & generator ) {

  const int m = pi.M();
  const int n = pi.N();
  const int k = pi.K();
  assert(k > 1);

  Chronometer timer;
  timer.start();

  if ( warmStar == false)
    pi.create_random_solution_and_clean_TS ( generator );

  cout << "LogInitTS\t"
      << " initObj "      << pi.solution_measure()
      << endl;

  Problem_instance best_pi = pi;

  cout << "TS:"
      << " conservative = " << pi.ops.APR.conservative
      << " samples = "<< pi.ops.TS.neighborhood_samples
      << " swaps = " << pi.ops.TS.swaps
      << endl;


  int iteration = 0;
  while ( timer.elapsed_time() < maxTime - TIMETOL ) {

    if (pi.ops.TS.flag_use_TS)
      pi.tabu_list++;

    vector<double> probs = ill_assigned_probabilities( pi, generator );

    Problem_instance best_neighboring_pi = pi;

    bool first_sol = true;

    cout << "LogIntTS\t"
        << " origNeighSol "      << pi.solution_measure()
        << endl;


    if ( pi.ops.TS.exhaustiveEnumeration == false ) {

      for ( int s = 0; s < pi.ops.TS.neighborhood_samples; s++) {  //build a new sol

        Problem_instance new_pi = pts_neighboring_solution( pi, probs, generator );

        if ( new_pi.solution_measure() < best_neighboring_pi.solution_measure() || first_sol == true ) { //store best neighboring sol
          best_neighboring_pi = new_pi;
        }
        first_sol = false;

        cout << "LogIntTS\t"
            << " newNeig "      << new_pi.solution_measure()
            << " bestNeig "      << best_neighboring_pi.solution_measure()
            << endl;


      } //neighborhood explored

    }
    else if ( pi.ops.TS.exhaustiveEnumeration == true ) {

      if ( pi.ops.TS.flag_use_TS ) {
        while ( everything_tabu ( pi ) == true ) {
          pi.tabu_list++;
        }
      }
      for ( int i = 0; i < m; i++) { //complete 1-exchange neighborhood: for all points, for all clusters

        for ( int j = 0; j < k; j++) {

          if ( j == pi.J ( i ) ||
               pi.tabu_list.does_aspire( i, j, pi.distance[i][j] ) == false ) {
            continue;
          }

          Problem_instance new_pi = ts_exchange_neighboring_solution( pi, i, j ); //the solution being created

          cout << "LogIntTS\t"
               << " currNeig "      << new_pi.solution_measure()
               << endl;

          if ( new_pi.solution_measure() < best_neighboring_pi.solution_measure() || first_sol == true ) { //store best neighboring sol
            best_neighboring_pi = new_pi;
          }
          first_sol = false;
        }
      } //neighborhood explored

    }//end of " if ( pi.ops.TS.exhaustiveEnumeration == true ) "

    bool worsening = false;
    if ( best_neighboring_pi.solution_measure() < best_pi.solution_measure() ) { //store best
      best_pi = best_neighboring_pi;
    } else
      worsening = true;

    pi = best_neighboring_pi; //next solution is build AROUND the current one...

    iteration++;

    //cout << "It: " << iteration << "\t" << obj << endl;
    cout << "LogCurrTS\t"
        << " currObj "      << pi.solution_measure()
        << " worsening "    << worsening
        << " bestObj "      << best_pi.solution_measure()
        << " iterations "   << iteration
        << " timeLeft "     << maxTime - timer.elapsed_time()
        << endl;

    // pi.tabu_list.full_print_tabu_list();

  } // end of the iterations

  cout << "LogFinalTS\t"
      // << " currObj "      << obj
      << " currObj "      << pi.solution_measure()
      << " bestObj "      << best_pi.solution_measure()
      << " iterations "   << iteration
      << endl;


  Return_value retval;
  retval.obj = best_pi.solution_measure();
  retval.time_to_best = timer.elapsed_time();
  retval.total_time = retval.time_to_best; //the same (last sol = best sol)
  retval.iterations_to_best = iteration;
  retval.total_iterations = retval.iterations_to_best;  //the same (last sol = best sol)

  return retval;
}




//---------------------------------------------------------------------------------------
// ILL-ASSIGNED CRITERION
//---------------------------------------------------------------------------------------

double
algorithm::ill_assigned_value ( int i, Problem_instance & pi, base_generator_type & generator ) {

  int k = pi.K();
  int j = pi.J( i );

  double value;

  if ( pi.ops.APR.reassignment_type == AMM_DISTANCES || pi.ops.APR.reassignment_type == DISTANCES) {

    value = pi.distance[i][j];

    if ( pi.ops.APR.reassignment_type == AMM_DISTANCES ) {
      double second_best_dist = 1e300;
      // int second_best_j;
      for (int jj = 0; jj < k; jj++) {
        if (jj == j)
          continue;
        if (pi.distance[i][jj] < second_best_dist) {
          second_best_dist = pi.distance[i][jj];
          // second_best_j = jj;
        }
      }
      //ammortizing
      value /= ( second_best_dist + 0.001 );    //0.001 acts as a shift factor
    }
  }
  else if ( pi.ops.APR.reassignment_type == RANDIA  )  {
    value = random_generator::random( generator );
  }

  return value;  
}


vector<bool>
algorithm::find_ill_assigned(Problem_instance & pi, int j, double alpha, base_generator_type & generator ) {

  const int m = pi.M();
  const int k = pi.K();

  vector<bool>ill_ass(m, false);
  vector<bool>ass(m, false);

  int m_j = pi.M(j);
  if (m_j == 0)
    return ill_ass;

  // calculate maximum rank per cluster, i.e., the maximum rank above which a point is dubbed 'ill-assigned'
  // double max_ranking = (m_j-1) * (1-alpha);
  double max_ranking = floor( m_j * alpha );

  if ( max_ranking == 0 ) { //all points are ill-assigned: skip the computations
    // cout << "Max_ranking is 0, nothing is ill-assigned: skipping the computations" << endl;
    return ill_ass; //nothing is ill-assigned
  }

  for (int i = 0; i < m; i++) { //ass is initialized as a copy of pi.assignment[*][j]
    if (pi.assignment[i][j] == true)
      ass[i] = true;
  }

  if ( max_ranking == m_j ) { //all points are ill-assigned: skip the computations
    // cout << "Max_ranking = m_j: skipping the computations" << endl;
    return ass; //everything is ill-assigned
  }

  vector<double> values(m, -1);

  for (int i = 0; i < m; i++) {
    if (pi.assignment[i][j] == false)
      values[i] = -1; //this way, they are easily recognized when debugging
    else
      values[i] = ill_assigned_value ( i, pi, generator );
  }

  vector<int> podium(m, -1); //podium[i] is the index of the point that classified in position i
  for (int i = 0; i < m; i++) {
    podium[i] = i;
  }

  // less_then lto;
  // lto.set_vectors(values);
  // sort(podium.begin(), podium.end(), lto);

  // for (int i = 0; i < m; i++) {
  //   ranking[podium[i]] = i -m + m_j; //remove 0 ranking (not assigned to j) points
  //   // ranking[podium[i]] = i; //remove 0 ranking (not assigned to j) points
  // }

  // // calculate maximum rank per cluster, i.e., the maximum rank above which a point is dubbed 'ill-assigned'
  // // double max_ranking = (m_j-1) * (1-alpha);
  // double max_ranking = ceil( (m_j-1) * (1-alpha) );


  // for (int i = 0; i < m; i++) {
  //   if (pi.assignment[i][j] == false)
  //     continue;
  //   // if (ranking[i] > max_ranking)
  //   if (ranking[i] >= max_ranking) //to allow for 100% ill-ass when alpha = 1
  //     ill_ass[i] = true;
  // }


  greater_then gto;
  gto.set_vectors(values);
  sort(podium.begin(), podium.end(), gto);


  for ( int i = 0; i < max_ranking; i++ ) {
    ill_ass[podium[i]] = true;
  }

  // cout << "max_ranking " << max_ranking << endl;
  // cout << "m(j)\t " << m_j << endl;
  // cout << "ass\t " << ass << endl;
  // cout << "values\t " << values << endl;
  // cout << "podiumm\t " << podium << endl;
  // cout << "ill_ass\t " << ill_ass << endl;  

  return ill_ass;
}

vector<double>
algorithm::ill_assigned_probabilities( Problem_instance & pi, base_generator_type & generator ) {

  const int m = pi.M();
  const int k = pi.K();

  vector<double> values (m, -1);

  for (int i = 0; i < m; i++) {
    values[i] = ill_assigned_value ( i, pi, generator );
  }

  //normalize
  utility::normalize_1norm( values );

  return values;
}

//---------------------------------------------------------------------------------------
// ILL-ASSIGNED REASSIGNMENT
//---------------------------------------------------------------------------------------

//CAVEAT: It does NOT consider the tabu list. If the suggested reassignment is tabu, it must be prevented in the caller.
// int
// algorithm::suggest_reassignment( Problem_instance & pi, int i, base_generator_type & generator ) {

//   int new_j;

//   int k = pi.K();

//   if ( pi.ops.APR.reassignment_target == CLOSEST ) { //CLOSEST different from current one
//     double best_dist = 1e300;

//     for (int jj = 0; jj < k; jj++) {
//       if (jj == pi.J( i ) )
// 	continue;
//       if (pi.distance[i][jj] < best_dist) {
// 	best_dist = pi.distance[i][jj];
// 	new_j = jj;
//       }
//     }
//   }
//   else if ( pi.ops.APR.reassignment_target == RAND_DIST ) {
//     vector<double> cluster_probs( k, -1);

//     for (int jj = 0; jj < k; jj++) {
//       if (jj == pi.J( i ) )
// 	cluster_probs[jj] = 0;
//       else
// 	cluster_probs[jj] = 1.0/pi.distance[i][jj]; //1 was amiss
//     }

//     utility::normalize_1norm( cluster_probs );

//     // do {
//     //   new_j = random_generator::sample_discrete_distribution( cluster_probs, generator );
//     // } while ( new_j == pi.J( i ) );
//     //find SEPTEMBER_ISSUE to read a comment on the same issue fixed at another place
//     new_j = random_generator::sample_discrete_distribution( cluster_probs, generator );

//   }
//   else if ( pi.ops.APR.reassignment_target == RAND ) {
//     new_j = random_generator::integer( k, generator );
//   }

//   // assert ( new_j != -1 );

//   return new_j;

// }

int
algorithm::suggest_nontabu_reassignment( Problem_instance & pi, int i, base_generator_type & generator ) {

  int new_j;

  int k = pi.K();

  int jjj = pi.J( i );

  int closest_cluster = -1; // -1 is returned if no suitable nontabu (and noncurrent) cluster is found

  if (pi.ops.APR.reassignment_target == CLOSEST) {
    double min_dist = 1e300;
    for (int j = 0; j < k; j++) {
      if ( jjj == j )
        continue;
      if ( pi.ops.TS.flag_use_TS && pi.tabu_list.is_tabu(i,j)
              && pi.tabu_list.does_aspire(i,j,pi.distance[i][j]) == false ) {
        continue;
      }
      else if (pi.distance[i][j] < min_dist) {
        min_dist = pi.distance[i][j];
        closest_cluster = j;
      }
    }
  }
  else if (pi.ops.APR.reassignment_target == RAND_DIST || pi.ops.APR.reassignment_target == RAND ) {
    vector <double> probs(k, 0);
    for (int j = 0; j < k; j++) {
      if ( jjj == j ||
          ( pi.ops.TS.flag_use_TS && pi.tabu_list.is_tabu(i,j)
              && pi.tabu_list.does_aspire(i,j,pi.distance[i][j]) == false )
      )
        probs[j] = 0;
      else {
        if ( pi.ops.APR.reassignment_target == RAND_DIST )
          probs[j] = 1.0/pi.distance[i][j];
        else if ( pi.ops.APR.reassignment_target == RAND )
          probs[j] = 1.0; //uniform distribution
      }
    }

    if ( utility::sum( probs ) > 1e-6 ) { //the probs are not null: an assignment can be performed

      utility::normalize_1norm( probs );

      //SEPTEMBER_ISSUE: modified on sept 1 2011: I allow to sample jjj again. Suppose k=2 and jjj-2; assigning point i to cluster j=1 is tabu; then, it must be assigned to cluster jjj=2. The other option is recognizing all options are tabu and moving the tabu iteration counter on. It amounts to just sample the cluster jjj=2 and "tanti saluti". Yes but, with this, the point might be reassigned to the current cluster with high probability
      //CHECK-NEO
      closest_cluster = random_generator::sample_discrete_distribution( probs, generator );

    }
  }


  return closest_cluster;
}

int
algorithm::suggest_conservative_nontabu_reassignment( Problem_instance & pi, int i, base_generator_type & generator ) {

  int new_j;
  int old_j = pi.J( i );

  assert ( pi.ops.APR.conservative == true );

  new_j = find_closest_nontabu_cluster( pi, i );
  if ( new_j == old_j )
    new_j = suggest_nontabu_reassignment( pi, i, generator );

  return new_j;
}



int
algorithm::APR_combinatorial_reassignment(Problem_instance & pi, int curr_j, vector<bool> & ill_ass, base_generator_type & generator ) {

  int new_j;
  int ill_assigned_found = 0;

  const int m = pi.M();

  for (int i = 0; i < m; i++) {
    if ( pi.assignment[i][curr_j] == false || ill_ass[i] == false )
      continue;
    else
      ill_assigned_found++;

    if ( pi.ops.APR.conservative == false ) {
      new_j = suggest_nontabu_reassignment( pi, i, generator );
      // cout << new_j << endl; 
    }
    else {
      new_j = suggest_conservative_nontabu_reassignment( pi, i, generator );
    }

    //the old pair is made tabu regardless of the objective function in PR!
    if (new_j != -1 && new_j != curr_j) {

      if ( pi.ops.TS.flag_use_TS )
        pi.tabu_list.add_old_and_remove_new( i, curr_j, pi.distance[i][curr_j], new_j );

      pi.assignment[i][curr_j] = false;
      pi.assignment[i][new_j] = true;
    }
  }

  return ill_assigned_found;
}



//----------------------------------------------------------
// unused in cor2011
//----------------------------------------------------------


void
algorithm::random_shaker(Problem_instance & pi, base_generator_type & generator ) {

  int m = pi.M();
  int k = pi.K();

  // vector<double> probs = ill_assigned_probabilities( pi ); //12-09-2011
  vector<double> probs = ill_assigned_probabilities( pi, generator );
  cout << "ill-assignment probs " << probs << endl;

  for ( int i = 0; i < m; i++) {

    if ( random_generator::random( generator ) <= probs[i] ) {

      cout << "Shaking point " << i << endl;

      int new_j = -1;

      if ( pi.ops.APR.reassignment_target == CLOSEST ) {
        double best_dist = 1e300;

        for (int jj = 0; jj < k; jj++) {
          if (jj == pi.J( i ) )
            continue;
          if (pi.distance[i][jj] < best_dist) {
            best_dist = pi.distance[i][jj];
            new_j = jj;
          }
        }
      }
      else if ( pi.ops.APR.reassignment_target == RAND_DIST ) {
        vector<double> cluster_probs( k, -1);

        for (int jj = 0; jj < k; jj++) {
          if (jj == pi.J( i ) )
            cluster_probs[jj] = 0;
          else
            cluster_probs[jj] = 1.0/pi.distance[i][jj]; //1 was amiss
        }

        utility::normalize_1norm( cluster_probs );

        // do {
        //   new_j = random_generator::sample_discrete_distribution( cluster_probs, generator );
        // } while ( new_j == pi.J( i ) );
        //find SEPTEMBER_ISSUE in this file to get a comment on the error that was fixed
        new_j = random_generator::sample_discrete_distribution( cluster_probs, generator );

      }
      else if ( pi.ops.APR.reassignment_target == RAND ) {
        new_j = random_generator::integer( k, generator );
      }

      assert ( new_j != -1 );

      // cout << i << "from " << pi.J(i) << " to " << new_j << endl;

      int old_j = pi.J(i);
      pi.assignment[i][pi.J(i)] = false;
      pi.assignment[i][new_j] = true;

      //! the GOOD assignment becomes tabu ....
      // cout << "TS before" << endl;
      // pi.tabu_list.full_print_tabu_list(); cout << endl;

      if ( old_j != new_j )
        pi.tabu_list.add( i, pi.J(i), pi.distance[i][pi.J(i)] );

      // cout << "TS after" << endl;
      // pi.tabu_list.full_print_tabu_list(); cout << endl;

    }
  }

  plain_parameter_update( pi );

  pi.update_distances( );
}

Return_value
algorithm::anytimeShakenBM( Problem_instance & pi, double maxTime, base_generator_type & generator ) {

  Chronometer timer;
  timer.start();

  pi.create_random_solution_and_clean_TS( generator );

  Return_value ret_global;
  Return_value ret_local;

  Problem_instance best_pi( pi );

  ret_global.obj = pi.solution_measure();
  ret_global.time_to_best = timer.elapsed_time();
  ret_global.iterations_to_best = 0;

  pi.create_random_solution_and_clean_TS( generator );

  cout << "LogInitSBM\t" 
      << " obj "          << pi.solution_measure()
      << endl;


  uint i;
  for ( i = 0; timer.elapsed_time() <= maxTime - TIMETOL; i++ ) {

    cout << "LogCurrSBM\t" 
        << " start "          << i+1
        << " remaining time " << maxTime - timer.elapsed_time()
        << endl;

    //base_generator_type generator( i );
    //Given generator must be used, otherwise this algorithm could not be multistarted (needed for statistical comparison purposes)


    ret_local = algorithm::BM( pi, true,
        std::max( maxTime - timer.elapsed_time(), 0.0 ),
        generator );

    if ( ret_local.obj < ret_global.obj ) {
      best_pi = pi;
      ret_global.obj = ret_local.obj;
      ret_global.time_to_best = ret_global.total_time + ret_local.time_to_best;
      ret_global.iterations_to_best = ret_global.total_iterations + ret_local.iterations_to_best;
    }

    ret_global.total_iterations += ret_local.total_iterations;
    ret_global.total_time = timer.elapsed_time();

    cout << "LogCurrSBM\t"
        << " bestObj "      << ret_global.obj
        << " currObj "      << ret_local.obj
        << " totIters "     << ret_global.total_iterations+1
        << " bestIters "    << ret_global.iterations_to_best+1
        << " time "         << ret_global.total_time
        << " bestTime "     << ret_global.time_to_best
        << " starts "       << i+1
        << endl;

    //apply shaker

    cout << "pre shaker solution value = " << pi.solution_measure() << endl;

    random_shaker ( pi, generator );

    cout << "post-shaker solution value = " << pi.solution_measure()  << endl;


  } //end of multi-starts



  cout << "LogGlobalSBM\t"
      << " bestObj "      << ret_global.obj
      << " totIters "     << ret_global.total_iterations+1
      << " bestIters "    << ret_global.iterations_to_best+1
      << " time "         << ret_global.total_time
      << " bestTime "     << ret_global.time_to_best
      << " starts "       << i+1
      << endl;

  pi = best_pi;

  ret_global.total_starts = i+1;

  return ret_global;
}
