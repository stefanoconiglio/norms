#include"inclusions.h"

// Note: most of the stuff on random number in the Boost Greph Library is taken from the example random_example.cpp (see Boost documentation)

// extern int nLeastEigPowersCalled;
// extern int nLeastEigPowersIterationsCalled;


using namespace boost;


int
main(int argc, char *argv[]) {
  // nLeastEigPowersCalled = 0;
  // nLeastEigPowersIterationsCalled = 0;

  //! usage output
  if (argc == 1) {
      cerr << "Example of use:" << endl;
      cerr << "/argv[0] -f ../data/cor2011/random-medium/sr_750_6_8_ins -m mBM -m 10" << endl;
      return 0;
  }

  string file_name;
  // string solution_file_name;

  int precision = 6;

  //! default setting of various parameters
  t_algorithm_choice algorithm;

  Problem_type problem_type;

  problem_type.problem = HCLUSTERING;

  t_alg_options ops;

  ops.power_tolerance = 0;
  ops.max_time = 1e300;  //caveat: cplex default; if given larger, Cplex throws an exception

  ops.flag_intertwine_RLP = false;
  ops.flag_RLP_based_combinatorial_reassignment = false;
  ops.flag_RLP_at_the_end = false;
  ops.RLP_based_combinatorial_reassignment_sol_limit = 2100000000; //CPLEX default

  ops.flag_generate_RLP_feasible_random_solutions = false;

  ops.flag_use_shakers = false;

  //REV1
  ops.bemporad_delta = 0; //if 0, Bemporad's criterion is not applied
  ops.bemporad_passes = 1; //if 0, Bemporad's criterion is not applied
  ops.bemporad_c = 10; //if 0, Bemporad's criterion is not applied

  ops.TS.flag_use_TS = false;
  //  ops.TS.flag_use_TS = true;
  ops.TS.length = 2;
  ops.TS.flag_use_aspiration_criterion = true;

  ops.TS.exhaustiveEnumeration = false;
  // ops.TS.neighborhood_samples = 1;
  // ops.TS.fraction_of_swaps = 0.05;
  ops.TS.neighborhood_samples = 1;
  ops.TS.swaps = 1;
  ops.TS.cooper = 0;
  ops.TS.diversify = false;
  ops.TS.diversification_swaps_percent = 0.5;


  ops.APR.flag_use_linear_threshold = false;
  ops.APR.reassignment_type = AMM_DISTANCES;
  // ops.APR.reassignment_type = DISTANCES;
  ops.APR.reassignment_target = CLOSEST;
  ops.APR.init_alpha = 1;
  ops.APR.restart_alpha = 0.6;
  ops.APR.min_alpha = 0.001;					//! it is CRUCIAL that is is nonzero, otherwise the obscillations in the objective function will NOT be identifies and the algorithm will not stop!
  //ops.APR.restart_alpha = ops.APR.init_alpha;
  //ops.APR.restart_alpha = 0.2;
  ops.APR.rho_down = 0.99; //8-7-2011
  ops.APR.rho_down = 0.999;
  //ops.APR.rho_up = pow(0.99,2);
  ops.APR.rho_up = 0.9;
  ops.APR.alpha_decrease_frequency = 1;
  ops.APR.extra_update = false;
  ops.APR.conservative = true;

  int o; //! current option
  bool flag_output_matlab_data = false;
  bool flag_output_solution = false;

  int given_k = -1;

  //bool flag_normalize_points = false;
  //bool flag_statistical_normalization;

  string algorithm_choice;

  //! awful reading var needed in the switch that Mr Nice-Guy C++ doesn't allow me to declare therein...
  int divisor_1, divisor_2, divisor_3, divisor_4, divisor_5, divisor_6, divisor_7, divisor_8, divisor_9, divisor_10;
  char tmp;
  int str_length;
  string str_tabu_choice;
  string str_aor_choice;
  string str_bemporad_choice;

  ops.multistarts = int(1e300);

  while ((o = getopt(argc, argv, "d:f:P:k:a:m:p:z:t:C:c:T:I:R:M::D:U:E:L:X:N:W:O:r:e::w::h::F::S:o::B:")) != -1) {
    switch (o) {
    case 'f':
      file_name = optarg;
      break;
    case 'P':
      if (!strcmp(optarg, "hclustering"))
	problem_type.problem = HCLUSTERING;
      else if (!strcmp(optarg, "affine_regression"))
	problem_type.problem = AFFINE_REGRESSION;
      else if (!strcmp(optarg, "linear_regression"))
	problem_type.problem = LINEAR_REGRESSION;
      else if (!strcmp(optarg, "pclustering"))
	problem_type.problem = PCLUSTERING;
      else if (!strcmp(optarg, "pamf"))
	problem_type.problem = PAMF;
      else {
	cout << "Unrecognized -P option" << endl;
	exit(-1);
      }
      break;
    case 'w':
      ops.flag_intertwine_RLP = true;
      cout << "ops.flag_intertwine_RLP set to " << ops.flag_intertwine_RLP << endl;
      break;
    case 'h':
      ops.flag_RLP_based_combinatorial_reassignment = true;
      cout << "ops.flag_RLP_based_combinatorial_reassignment to " << ops.flag_RLP_based_combinatorial_reassignment << endl;
      break;
    case 'e':
      ops.flag_RLP_at_the_end = true;
      cout << "ops.flag_RLP_at_the_end " << ops.flag_RLP_at_the_end << endl;
      break;
    case 'S':
      ops.RLP_based_combinatorial_reassignment_sol_limit = atoi( optarg );
      break;
    case 'k':
      given_k = atoi(optarg);
      break;
    case 'a':
      algorithm_choice = optarg;
      if (!strcmp(optarg, "PAPAVERO"))
	algorithm = t_PAPAVERO;
      else if (!strcmp(optarg, "mPAPAVERO"))
	algorithm = t_mPAPAVERO;
      else if (!strcmp(optarg, "maPAPAVERO"))
	algorithm = t_maPAPAVERO;
      else if (!strcmp(optarg, "BM"))
	algorithm = t_BM;
      else if (!strcmp(optarg, "PR"))
	algorithm = t_PR;
      else if (!strcmp(optarg, "mBM"))
	algorithm = t_mBM;
      else if (!strcmp(optarg, "mPR"))
	algorithm = t_mPR;
      else if (!strcmp(optarg, "aBM"))
	algorithm = t_APR;
      else if (!strcmp(optarg, "maBM"))
	algorithm = t_maBM;
      else if (!strcmp(optarg, "maPR"))
	algorithm = t_maPR;
      else if (!strcmp(optarg, "mAPR"))
	algorithm = t_mAPR;
      else if (!strcmp(optarg, "PW"))
	algorithm = t_PW;
      else if (!strcmp(optarg, "mPW"))
	algorithm = t_mPW;
      else if (!strcmp(optarg, "aPW"))
	algorithm = t_maPW;
      else if (!strcmp(optarg, "sBM"))
	algorithm = t_sBM;
      else if (!strcmp(optarg, "msBM"))
	algorithm = t_msBM;
      // else if (!strcmp(optarg, "TS"))
      else if (!strcmp(optarg, "TS"))
	algorithm = t_TS;
      ////// else if (!strcmp(optarg, "mTS"))
      else if (!strcmp(optarg, "mTS"))
	algorithm = t_mTS; 
      else if (!strcmp(optarg, "TSBM"))
	algorithm = t_TSBM;
      else if (!strcmp(optarg, "mTSBM"))
	algorithm = t_mTSBM; 
      else if (!strcmp(optarg, "TSpamf"))
	algorithm = t_TSpamf; 
      else if (!strcmp(optarg, "mTSpamf"))
	algorithm = t_mTSpamf; 
      else if (!strcmp(optarg, "mJustRandom"))
	algorithm = t_mJustRandom; 
      else {
	cout << "Wrong algorithm parameter specified" << endl;
	exit(-1);
      }
      break;
    case 'm':
      ops.multistarts = atoi(optarg);
      break;
    case 'p':
      precision = atoi(optarg);
      break;
    case 'z':
      ops.power_tolerance = atof(optarg);
      break;
    case 't':
      ops.TS.flag_use_TS = true;
      str_tabu_choice = optarg;
      str_length = str_tabu_choice.length();
      
      divisor_1 = str_tabu_choice.find_first_of(":");
      ops.TS.length = atoi(str_tabu_choice.substr(0, divisor_1).c_str());
      
      if (ops.TS.length == -1)
        ops.TS.flag_use_TS = false;


      divisor_2 = str_tabu_choice.find_first_of(":", divisor_1+1);
      if (str_tabu_choice.substr(divisor_1+1, divisor_2-divisor_1-1) == "aspire")
	ops.TS.flag_use_aspiration_criterion = true;
      else if  (str_tabu_choice.substr(divisor_1+1, divisor_2-divisor_1-1) == "dont_aspire")
	ops.TS.flag_use_aspiration_criterion = false;
      else {
	cout << "Wrong -t option specified" << endl;
	exit(-1);
      }
      // cout << "Tabu search options: " << ops.TS.length << ":" << (ops.TS.flag_use_aspiration_criterion ? "aspire" : "don't aspire") << endl;
      break;
    case 'C':
      if (!strcmp(optarg, "distances"))
	ops.APR.reassignment_type = DISTANCES;
      else if (!strcmp(optarg, "amm_distances"))
	ops.APR.reassignment_type = AMM_DISTANCES;
      else if (!strcmp(optarg, "random"))
	ops.APR.reassignment_type = RANDIA;
      else {
	cout << "Unrecognized -C option" << endl;
	exit(-1);
      }
      break;
    case 'c':
      if (!strcmp(optarg, "random"))
	ops.APR.reassignment_target = RAND;
      else if (!strcmp(optarg, "random_distances"))
	ops.APR.reassignment_target = RAND_DIST;
      else if (!strcmp(optarg, "closest"))
	ops.APR.reassignment_target = CLOSEST;
      else {
	cout << "Unrecognized -c option" << endl;
	exit(-1);
      }
      break;
    case 'T':
      ops.max_time = atof(optarg);
      break;
      //APR parameters
    case 'I':
      ops.APR.init_alpha = atof( optarg );
      break;
    case 'R':
      ops.APR.restart_alpha = atof( optarg );
      break;
    // case 'M':
    //   ops.APR.min_alpha = atof( optarg );
    //   break;
    case'M':
      flag_output_matlab_data = true;
      break;
    case 'o':
      flag_output_solution = true;
      // solution_file_name = optarg;
      break;
    case 'D':
      ops.APR.rho_down = atof( optarg );
      break;
    case 'U':
      ops.APR.rho_up = atof( optarg );
      break;
    // case 'F':
    //   ops.APR.alpha_decrease_frequency = atoi( optarg );
    //   break;
    case 'F':
      ops.flag_generate_RLP_feasible_random_solutions = true;
      break;
    case 'E':
      ops.APR.extra_update = atoi( optarg );
      break;
    case 'L':
      ops.APR.flag_use_linear_threshold = atoi ( optarg );
      break;
    case 'X':
      ops.TS.exhaustiveEnumeration = atoi( optarg );
      break;
    case 'N':
      ops.TS.neighborhood_samples = atoi( optarg );
      break;
    case 'W':
      ops.TS.swaps = atof( optarg );
      break;
    case 'O':
      ops.APR.conservative = atoi( optarg );
      break;
    case 'r':
      ops.TS.cooper = atoi( optarg );
      break;
    case 'd':
      ops.TS.diversify = true;
      ops.TS.diversification_swaps_percent  = atof( optarg );
      break;
    case 'B':
      ops.bemporad_delta = atof ( optarg );
      str_bemporad_choice = optarg;
      
      divisor_1 = str_bemporad_choice.find_first_of(":");
      ops.bemporad_delta = atof ( str_bemporad_choice.substr(0, divisor_1).c_str());
      
      divisor_2 = str_bemporad_choice.find_first_of(":", divisor_1+1);
      ops.bemporad_passes = atof ( str_bemporad_choice.substr(divisor_1+1, divisor_2-divisor_1-1).c_str() );

      divisor_3 = str_bemporad_choice.find_first_of(":", divisor_2+1);
      ops.bemporad_c = atof ( str_bemporad_choice.substr(divisor_2+1, divisor_3-divisor_2-1).c_str() );
    }
  }
  
  
  cout << "--------------PARAMETERS-----------------------------------------" << endl;

  cout << "problem type "                         << problem_type.problem << endl;

  cout << "T use = "                              << ops.TS.flag_use_TS << endl;
  cout << "T tabu list length  = "                << ops.TS.length << endl;
  cout << "T tabu list aspire  = "                << ops.TS.flag_use_aspiration_criterion << endl;


  cout << "I init_alpha = "                       << ops.APR.init_alpha << endl;
  cout << "R restart_alpha = "                    << ops.APR.restart_alpha << endl;
  cout << "M min_alpha (fixed) = "                << ops.APR.min_alpha << endl;
  cout << "D rho_down = "                         << ops.APR.rho_down << endl;
  cout << "U rho_up = "                           << ops.APR.rho_up << endl;
  cout << "F alpha_decrease_frequency (fixed) = " << ops.APR.alpha_decrease_frequency << endl;
  cout << "E extra_update = "                     << ops.APR.extra_update << endl;

  cout << "L linear_threshold = "                 << ops.APR.flag_use_linear_threshold << endl;
  
  cout << "Bemporad's criterion delta = "         << ops.bemporad_delta << endl;
  cout << "Bemporad's criterion passes = "        << ops.bemporad_passes << endl;
  cout << "Bemporad's criterion c = "             << ops.bemporad_c << endl;

  cout << "-----------------------------------------------------------------" << endl;

  //random number generator with a GIVEN SEED of 42u
  base_generator_type generator(42u);

  //! set the given precision
  cout.precision(precision);

  //! populates a new problem_instance from file
  cout << "Populating a new problem solution from file " << file_name << "...";
  Problem_instance pi(file_name, problem_type, ops, given_k);

//  if (ops.TS.length == -1) {
//    pi.ops.TS.length = int(double(pi.M())/100.0);
//    if (pi.ops.TS.length < 2)
//      pi.ops.TS.length = 2;
//    cout << "TS length set to " << pi.ops.TS.length << endl;
//  }
  
  cout << "... done!" << endl;

  //for TS and TSpamf: if ops.TS.swaps < 0, it uses it as a percentage of pi.M();
  if ( ops.TS.swaps < 0 && ops.TS.swaps > -1 ) { //the -1 case is considered as: "sample all"; this is computationally lighter and is done in a special way
    cout << "Converting ops.TS.swaps from " << ops.TS.swaps;
    ops.TS.swaps = floor( pi.M() * (-ops.TS.swaps) );
    cout << " to " << ops.TS.swaps << endl;
    pi.ops.TS.swaps = ops.TS.swaps;
  }


  if ( pi.K() == 1) { //given_k ?
    cout << "Only one cluster has been required: closed-form solution being returned" << endl;
    algorithm::compute_trivial_one_cluster_solution(pi );
    exit(0);
  }
  
  switch ( algorithm ) {
      // algorithms that naturally terminate (ATNT)
    case t_BM:
      algorithm::BM( pi, false, ops.max_time, generator );
      break;
    case t_PR:
      algorithm::PR( pi, 1, false, ops.max_time, generator );
      break;
      // ATNT, multistarted to evaluate, on average, the quality of the solutions they find
    case t_PAPAVERO:
      algorithm::papavero( pi, false, ops.max_time, generator );
      break;
    case t_mBM:
      algorithm::multistarter( pi, t_mBM, ops.max_time, generator );
      break;
    case t_mPR:
      algorithm::multistarter( pi, t_mPR, ops.max_time, generator );
      break;
      // anytime versions of the ATNT (within a given timelimit, they go on multistarting, then return the best solution)
      // APR, which is an anytime algorithm
    case t_mPAPAVERO:
      algorithm::multistarter( pi, t_mPAPAVERO, ops.max_time, generator );
      break;
    case t_APR:
      algorithm::APR( pi, false, ops.max_time, generator );
      break;
      // anytime versions of the ATNT, multistarted to evaluate, on average, the quality of the solutions they find
    case t_maBM:
      algorithm::multistarter( pi, t_maBM, ops.max_time, generator );
      break;
    case t_maPR:
      algorithm::multistarter( pi, t_maPR, ops.max_time, generator );
      break;
    case t_maPAPAVERO:
      algorithm::multistarter( pi, t_maPAPAVERO, ops.max_time, generator );
      break;
    case t_mAPR:
      algorithm::multistarter( pi, t_mAPR, ops.max_time, generator );
      break;
    case t_PW:
      algorithm::PW( pi, false, ops.max_time, generator );
      break;
    case t_mPW:
      algorithm::multistarter( pi, t_mPW, ops.max_time, generator );
      break;
    case t_maPW:
      algorithm::multistarter( pi, t_maPW, ops.max_time, generator );
      break;
    case t_sBM:
      algorithm::anytimeShakenBM( pi, ops.max_time, generator );
      break;
    case t_msBM:
      algorithm::multistarter( pi, t_msBM, ops.max_time, generator );
      break;
    case t_TS:
      algorithm::TS( pi, false, ops.max_time, generator );
      break;
    case t_mTS:
      algorithm::multistarter( pi, t_mTS, ops.max_time, generator );
      break;
    case t_TSBM:
      algorithm::TSBM( pi, false, ops.max_time, generator );
      break;
    case t_mTSBM:
      algorithm::multistarter( pi, t_mTSBM, ops.max_time, generator );
      break;
    case t_TSpamf:
      algorithm::TSpamf( pi, false, ops.max_time, generator );
      break;
    case t_mTSpamf:
      algorithm::multistarter( pi, t_mTSpamf, ops.max_time, generator );
      break;
    case t_mJustRandom:
      algorithm::multistarter( pi, t_mJustRandom, ops.max_time, generator );
      break;
    default:
      cout << "t_algorithm incompatible with main" << endl;
      exit(-1);
    }

  // cout << "Solution after the algorithms are done " << pi.solution_measure() << endl << pi.assignment << endl;

  // if ( !pi.ops.flag_intertwine_RLP && !pi.ops.flag_RLP_based_combinatorial_reassignment && (pi.problem_type.problem == AFFINE_REGRESSION || pi.problem_type.problem == LINEAR_REGRESSION) ) { 
  //   cout << "RLP was not used in intertwining not directly in the combinatorial reassignment. Reassigning with RLP now, then updating the parameters and the distances" << endl;
  //   cout << "Assignment before calling pamf::reassign_with_RLP( pi ); " << pi.assignment << endl;
  //   pamf::reassign_with_RLP( pi );
  //   cout << "Assignment after calling pamf::reassign_with_RLP( pi ); " << pi.assignment << endl;
  //   algorithm::plain_parameter_update( pi );
  //   pi.update_distances();
  // }


  //cout << pi << endl;
  // cout << "Final solution before returning " << pi.solution_measure() << endl << pi.assignment << endl;



  // output to matlab and to file
  if (flag_output_matlab_data && pi.N() == 2) {
    string matlab_data_file_name(file_name);
    matlab_data_file_name += string("_sol.m");
    cout << "Outputting to matlab on file " << matlab_data_file_name << endl;
    pi.output_solution_to_matlab(matlab_data_file_name, true); //last parameter (true) means: print cluster lines = {on, off}

    //cout << pi << endl;
  }

  if ( pi.assert_correctness() == false ) {
    cout << "pi.assert_correctness() == false" << endl;
    cout.flush();
    exit(-1);
  }

  if ( flag_output_solution ) {
    cout << pi.K() << endl;
    string solution_file_name( file_name );
    solution_file_name += "_sol_k";
    char tmp[2];
    sprintf( tmp, "%i", pi.K() );
    solution_file_name += tmp;
    pi.output_solution_to_file( solution_file_name );
    cout << "Outputting solution to file " << solution_file_name << endl;
  }



}
