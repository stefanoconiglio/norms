#include "algorithm.h"

using namespace std;


Return_value
algorithm::APR( Problem_instance & pi, bool warmStart, double maxTime, base_generator_type & generator ) {

  const int k = pi.K();
  assert(k > 1);
  const int m = pi.M();

  Problem_instance best_pi(pi);

//  int noobjectivechange = 0;
//  double old_obj;

  //double alpha = pi.ops.APR.init_alpha;
  double current_restart_alpha = pi.ops.APR.restart_alpha;

  Chronometer timer;
  timer.start();
  //double init_time = timer.elapsed_time();

  bool alpha_peak_already_reached = false;
  double v_at_peak_alpha = 1e300;

  //variables for statistics
  //vector<double> obj_per_start;
  //vector<double> iter_per_start;
  //vector<double> time_per_start;
  //double best_occurr = 0;

  //double timeOfLastStart = init_time;

  Return_value retval;

  //initial random solution
  if ( warmStart == false ) {
    pi.create_random_solution_and_clean_TS ( generator );
  }

  best_pi = pi;

  retval.obj = pi.solution_measure();
  retval.time_to_best = timer.elapsed_time();
  retval.total_time = retval.time_to_best; //the same (init sol = best sol)
  retval.iterations_to_best = 0;
  retval.total_iterations = retval.iterations_to_best;  //the same (init sol = best sol)

  cout << "LogInitAPR\t"
      << " initObj "      << pi.solution_measure()
      << endl;

  // PR iterations (increasing-decreasing alpha)
  int i; //needed to print it outside of the cycle
  //for ( i = 0; /*i < pi.ops.multistarts*/; i++ ) {
  for ( i = 0; timer.elapsed_time() <= maxTime - TIMETOL; i++ ) {

    //if ( pi.chrono.elapsed_time()-init_time >= pi.ops.max_time ) {
    //if ( timer.elapsed_time() - init_time >= maxTime  ) {
    // if ( timer.elapsed_time() >= maxTime  ) {
    //   cout << "APR: Time out" << endl;
    //   break;
    // }

    //anti-cycle at peak_alpha

    if ( current_restart_alpha == pi.ops.APR.init_alpha ) {
      if (alpha_peak_already_reached == false) {
        v_at_peak_alpha = pi.solution_measure();
        alpha_peak_already_reached = true;
        cout << "APR: Alpha-peak reached. Solution = " << pi.solution_measure() << endl;
      }
      else if (alpha_peak_already_reached == true) {
        if (v_at_peak_alpha == pi.solution_measure() ) {
          // pi.randomly_populate_clusters( generator );
          // pi.update_distances();
          // plain_combinatorial_reassignment(pi);
          pi.create_random_solution_and_clean_TS ( generator ); //TS was not emptied
          cout << "APR: Alpha-peak loop. Random solution set = " << pi.solution_measure() << endl;
        }
        alpha_peak_already_reached = false;
        v_at_peak_alpha = 1e300;
      }
    }



    //calling PR
    Return_value PR_ret = algorithm::PR( pi, current_restart_alpha, true, std::max ( maxTime - timer.elapsed_time(), 0.0 ), generator );

    //obj_per_start.push_back( pi.solution_measure() );
    //time_per_start.push_back( timer.elapsed_time() - timeOfLastStart );
    //iter_per_start.push_back( PR_ret.total_iterations );

    if (pi.solution_measure() < best_pi.solution_measure() ) {
      best_pi = pi;
      retval.obj = best_pi.solution_measure();
      retval.time_to_best = retval.total_time + PR_ret.time_to_best;
      retval.iterations_to_best = retval.total_iterations + PR_ret.iterations_to_best;
      current_restart_alpha = pi.ops.APR.restart_alpha;
      //alpha = current_restart_alpha;
      //best_occurr = 1;
    }
    else {
      current_restart_alpha = current_restart_alpha/pi.ops.APR.rho_up < pi.ops.APR.init_alpha ? current_restart_alpha/pi.ops.APR.rho_up : pi.ops.APR.init_alpha;
      //alpha = current_restart_alpha;
      //best_occurr++;
    }

    //retval.total_time = timer.elapsed_time() - init_time;
    retval.total_time = timer.elapsed_time();
    //retval.total_iterations = i;
    retval.total_iterations += PR_ret.total_iterations;

    cout << "LogLocalAPR\t"
        << " bestObj "      << best_pi.solution_measure()
        << " currObj "      << pi.solution_measure()
        << " restart_alpha "<< current_restart_alpha
        << " PRstarts "     << i+1
        << endl;

//    if ( pi.solution_measure() == old_obj )
//      noobjectivechange++;
//    else
//      noobjectivechange = 0;

//    old_obj = pi.solution_measure();

    //if ( noobjectivechange >= 1 ) {
    //  cout << "APR: Objective has stalled! Terminating" << endl;
    //  break;
    //}


  } //end of meta-iterations

  cout << "LogFinalAPR\t"
      << " bestObj "      << best_pi.solution_measure()
      << " currObj "      << pi.solution_measure()
      << " restart_alpha "<< current_restart_alpha
      << " PRstarts "     << i+1
      << endl;

  retval.total_starts = i+1;

  return retval;
} //end of algorithm::APR()

Problem_instance
algorithm::ts_exchange_neighboring_solution( const Problem_instance & pi, int i, int j ) {

  Problem_instance new_pi = pi;

  int old_j = new_pi.J( i );
  double old_distance = new_pi.distance[i][old_j];

  new_pi.assignment[i][old_j] = false;
  new_pi.assignment[i][j] = true;

  plain_parameter_update( new_pi, old_j );
  plain_parameter_update( new_pi, j );

  new_pi.update_distances( old_j );
  new_pi.update_distances( j );

  if ( new_pi.solution_measure() > pi.solution_measure() ) //worsening solution
    new_pi.tabu_list.add_old_and_remove_new( i, old_j, old_distance, j );

  return new_pi;
}

