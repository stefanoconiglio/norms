#include"inclusions.h"

#include<getopt.h>

int
main ( int argc, char *argv[] ) {

  if (argc == 1) {
      cerr << "Example of use:" << endl;
      cerr << "/argv[0] -t sr -m 750 -n 6 -k 8 -s 666 [name is automatically given]" << endl;
      return 0;
  }

  
  string str_file_name = "";
  string string_type = "";
  int precision = 6;	//! cout precision
  int m = 20;
  int n = 2;
  int k = 3;

  Instance_type instance_type;
  //! default initialization
  instance_type.data_type = RANDOM;
  instance_type.variance = 0.003;
  instance_type.fraction_of_random_points = 0;
  // instance_type.fraction_of_misclassified = 0;
  instance_type.sigma_level_threshold = 0; //if 0: not used
  instance_type.error = GAUSSIAN;
  // instance_type.minimum_population_guarantee = false;

  int o;					//! current option
	
  bool flag_file_name_given       = false;
  bool flag_output_matlab_data    = false;
  bool flag_output_ampl_data      = false;
  bool b_variance_given           = false;		//randomizes a fraction of the data set
  // bool flag_add_misclassified_points = false;
  bool flag_add_outsiders         = false;
  bool flag_normalize             = true;
  double normalization_factor     = 10;
  bool flag_cutoffprecision       = true;
  double cutoffprecision          = 1e02;

  int seed = 0;

  string file_name;

  while ((o = getopt(argc, argv, "f:t:s:p:v:m:M::m:n:k:N:P:e:U::G::A::")) != -1)
    {
      switch (o)
	{
	case 'f':
          str_file_name = optarg;
	  flag_file_name_given = true;
	  break;
	case 't':
          string_type = optarg;
	  if (!strcmp(optarg, "r"))
	    instance_type.data_type = RANDOM;
	  else if (!strcmp(optarg, "sr"))
	    instance_type.data_type = SEMI_RANDOM;
	  else if (!strcmp(optarg, "srcr"))
	    instance_type.data_type = SEMI_RANDOM_CONTINUOS_REGRESSION;
	  else if (!strcmp(optarg, "srncr"))
	    instance_type.data_type = SEMI_RANDOM_NON_CONTINUOS_REGRESSION;
	  else if (!strcmp(optarg, "wave"))
	    instance_type.data_type = WAVE;
	  else if (!strcmp(optarg, "rwave"))
	    instance_type.data_type = RWAVE;
	  else if (!strcmp(optarg, "pclustering"))
	    instance_type.data_type = PCLUSTERING_INSTANCE;
	  else {
	    cout << "Unrecognized instance type" << endl;
	    abort();
	  }
	  break;
	case 'v':
          b_variance_given = true;
	  instance_type.variance = atof(optarg);
	  break;
	case 'e':
	  instance_type.sigma_level_threshold = atof(optarg);
	  break;
	case 's':
	  cout << "seed = " << optarg << "(as read), " << atoi(optarg) << " (atoi)" << endl;
          // random_generator::set_seed(atoi(optarg));
	  // cout << "Seed set to " << random_generator::get_seed() << endl;
	  seed = atoi( optarg );
	  break;
	case 'p':
          precision = atoi(optarg);
	  break;
	case 'm':
	  m = atoi(optarg);
	  break;
	case 'n':
	  n = atoi(optarg);
	  break;
	case 'k':
	  k = atoi(optarg);
	  break;
	case 'M':
          flag_output_matlab_data = true;
	  break;
	// case 'c':
        //   flag_add_misclassified_points = true;
	//   instance_type.fraction_of_misclassified = atoi(optarg);
	//   break;	
	case 'a':
          flag_add_outsiders = true;
	  instance_type.fraction_of_random_points = atoi(optarg);
	  break;
	case 'U':  //uses uniform error distrib. instead of guassian one.
	  instance_type.error = UNIFORM;
	  break;
	// case 'G':  //guarantees a minimum population per each cluster
	//   instance_type.minimum_population_guarantee = true;
	//   break;
       case 'A':
          flag_output_ampl_data = true;
          break;
	case 'N':
	  flag_normalize = true;
	  normalization_factor = atof(optarg);
	  break;
	case 'P':
	  flag_cutoffprecision = true;
	  cutoffprecision = atof(optarg);
	}
    }


  base_generator_type generator( seed );

  //default variance
  if (b_variance_given == false) {
    double dev_std = double(1)/(2*4 + 6*4 + 8-1); // 1/(2*2 + (k-2)*4 + k-1), with k=max_k=8
    instance_type.variance = pow(dev_std,2);
    cout << "Using a variance of " << instance_type.variance << endl;
  }
  
  //! set the given precision
  cout.precision(precision);

  //! creating an empty Problem_instance with the right data structure dimensions
  Problem_instance pi;
  Instance_generation_printing_data igpd;

  pi = Problem_instance(m, n, k);
	
  //! creating the default file name		
  if (flag_file_name_given) {
    file_name = string(str_file_name);
  }
  else {//automatic name generation 	
    char auto_file_name[256];
    sprintf(auto_file_name, "%s_%i_%i_%i_%i", string_type.c_str(), m, n, k, seed);
    
    if (flag_add_outsiders)
      sprintf(auto_file_name, "%s_o%f", auto_file_name, instance_type.fraction_of_random_points);
    
    // if (flag_add_misclassified_points)
    //   sprintf(auto_file_name, "%s_m%f", auto_file_name, instance_type.fraction_of_random_points);
    file_name = string(auto_file_name);
  }
  
  igpd.variances = vector<double>(k, 0);
  igpd.name = string_type;
  igpd.seed = seed;


  Instance_generator::populate( pi, instance_type, igpd, generator );

  if (flag_normalize)
    pi.normalize(normalization_factor);

  if (flag_cutoffprecision)
    pi.cutoffprecision(cutoffprecision);

  cout << "Instance generated with: m=" << m << ", n=" << n << ", k=" << k << endl;

  file_name = file_name + "_ins";

  //! output to file
  Instance_generator::write_points_to_file( pi, file_name, igpd );

  cout << "File written to: " << file_name << endl;

  //! output to ampl
  if (flag_output_ampl_data) {
    cout << "Outputting to ampl" << endl;
    string ampl_data_file_name(file_name);
    ampl_data_file_name += string(".dat");
    pi.output_standard_ampl_data(ampl_data_file_name);
  }

  //! output to matlab
  if (flag_output_matlab_data && n == 2) {
    cout << "Outputting to matlab" << endl;
    string matlab_data_file_name(file_name);
    matlab_data_file_name += string(".m");
    pi.output_solution_to_matlab(matlab_data_file_name, true); //last parameter (true) means: print cluster lines = {on, off}
  }
}

