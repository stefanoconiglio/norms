#ifndef __INCLUSIONS__
#define __INCLUSIONS__


#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
	
#include <fstream>
	
#include <iostream>
#include <string.h>
	
//#include <gsl/gsl_vector.h>
//#include <gsl/gsl_matrix.h>

#ifdef __ORLAB__
#include <lapackpp.h>
#else
#include <lapackpp/lapackpp.h>
#endif
	
#include <cmath>
	
#include <string>
#include <vector>
#include <deque>

#include <algorithm>

#include <sys/times.h>

#include "timer.h"

#include <boost/random.hpp>
// This is a typedef for a random number generator.
using namespace boost;
typedef boost::mt19937 base_generator_type;

#include <boost/random/discrete_distribution.hpp>

#define TIMETOL 0.001


//! problem enums
enum t_problem {HCLUSTERING, AFFINE_REGRESSION, LINEAR_REGRESSION, PCLUSTERING, PAMF};

// makeinstance enums
enum t_error {GAUSSIAN, UNIFORM};

// NOTE: ALL the fields MUST be filled in with something. Most important, i.e., for the tabu list.
// EVen if it is NOT GOING TO BE USED, all its fields must be coherent!!

struct t_TS_options
{	
  bool flag_use_TS;
  int length;
  bool flag_use_aspiration_criterion;
  bool exhaustiveEnumeration;
  /* double neighborhood_samples; */
  /* double fraction_of_swaps; */
  int neighborhood_samples;
  /* int swaps; */
  double swaps;  //2Feb2012
  int cooper;
  bool diversify;
  double diversification_swaps_percent;
};

//! reassignment type (in the combinatorial move)
enum t_reassignment_type {DISTANCES, AMM_DISTANCES, RANDIA};

//! target of the ill_assigned points
//! - RAND (but different from the current one)
//! - RAND_DIST: assigned on a random cluster (different fro the current one) on the basis of the inverse of the distance
//! - CLOSEST: reassigned to the closest cluster; to the second closest one if the point is ALREADY on the closest one
enum t_reassignment_target {RAND, RAND_DIST, CLOSEST};

struct t_APR_options
{
  t_reassignment_type reassignment_type;
  t_reassignment_target reassignment_target;
  double init_alpha;
  double rho_down;
  double rho_up;
  double min_alpha;
  int alpha_decrease_frequency;
  double restart_alpha;
  bool flag_use_linear_threshold;
  bool extra_update;
  bool conservative;
};

struct t_alg_options
{
  int multistarts;
  double max_time;
  double power_tolerance;		// hysteresis always uses a fixed tolerance for the whole run of the algorithm
  t_TS_options TS;
  bool flag_intertwine_RLP;
  bool flag_RLP_based_combinatorial_reassignment;
  bool flag_RLP_at_the_end;
  t_APR_options APR;
  bool flag_use_shakers;
  bool flag_generate_RLP_feasible_random_solutions;
  int RLP_based_combinatorial_reassignment_sol_limit;
  double bemporad_delta; //REV1: if 0, Bemporad's criterion is not applied; if > 0: it is applied at the end; if < 0, at each step
  int bemporad_passes; //REV1: number of passes for the criterion
  double bemporad_c; //REV1: number of points
};


struct Problem_type {			// information that DO NOT CHANGE for the whole execution of the alg are stored here
  t_problem problem;
};


//! available data construction types
typedef enum { 
  RANDOM,
  SEMI_RANDOM,
  SEMI_RANDOM_CONTINUOS_REGRESSION,
  SEMI_RANDOM_NON_CONTINUOS_REGRESSION,
  WAVE,
  RWAVE,
  PCLUSTERING_INSTANCE
} t_data_type;

struct Instance_type
{
  t_data_type data_type;
  double variance;
  //double fraction_of_misclassified;
  double fraction_of_random_points; //old add_noise
  double sigma_level_threshold;     //if 0, not used
  t_error error;
  /* bool minimum_population_guarantee; */
};

//! available algorithms
typedef enum { 
  t_PAPAVERO,
  t_mPAPAVERO,
  t_aPAPAVERO,
  t_maPAPAVERO,
  t_BM,
  t_mBM,
  t_aBM, //cannot be called directly from main
  t_maBM,
  t_PR,
  t_mPR,
  t_aPR, //cannot be called directly from main
  t_maPR,
  t_APR,
  t_mAPR,
  t_PW,
  t_mPW,
  t_aPW, //cannot be called directly from main
  t_maPW,
  t_sBM,
  t_msBM,
  t_TS,
  t_mTS,
  t_TSBM,
  t_mTSBM,
  t_TSpamf,
  t_mTSpamf,
  t_mJustRandom
} t_algorithm_choice;


//! instance generation options: stored inside the problem instance for printing purposes
struct Instance_generation_printing_data
{
  double variance;
  std::vector<double> variances;
  std::string name;
  double seed;
};


class Return_value {
 public:
  double obj;
  double time_to_best;
  int iterations_to_best;
  double total_time;
  int total_iterations;
  int total_starts;

  Return_value () {
    obj = 1e300;
    time_to_best = 0;
    iterations_to_best = 0;
    total_time = 0;
    total_iterations = 0;
    total_starts = 0;
  }
};

	
#include "ampl_output_reader.h"

#include "random.h"
	
#include "utility.h"

#include "formulas.h"
#include "maths.h"

#include "tabu_list.h"

#include "problem_instance.h"
#include "instance_generator.h"
	
#include "exact.h"

#include "pamf.h"
	
#include "algorithm.h"

#endif

