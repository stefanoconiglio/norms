#include "inclusions.h"

using namespace std;

#ifndef __PROBLEM_SOLUTION__
#define __PROBLEM_SOLUTION__


class Problem_instance {
  friend class cross_validation;
  friend class ampl;

  // DATA STRUCTURES
 public:
  // data points
  vector<LaVectorDouble> a;

  // clusters
  vector<LaVectorDouble> w;
  vector<double> gamma;
  vector<LaVectorDouble> W;
  vector<double> G;
	
  // settings --norm, sos or non sos obj function, clustering/regression
  Problem_type problem_type;
  t_alg_options ops;

  // solution data: assignment, distances...
  vector<vector<double> >	distance;
  vector<vector<bool> >	assignment;

  full_tabu_list tabu_list;	// no always used but who cares? :-)

  string file_name;

  // PUBLIC INTERFACE: constructors and internal data structure compilers
 public:
  //! empty constructor
  Problem_instance() {};

  //! constructor from data file; if k is not provided, it ia assumed the one read from the data file is valid; if sample_ratio is provided, the given % of points only are  actually loaded from file (useful for digital images with a huge amount of points)
  Problem_instance(string & file_name, Problem_type & type, t_alg_options & ops, int k = -1, int sample_ratio = 1);

  //! 'emmpty' constructor; it only gives the data structures the correct dimensions; used by Instance_generator
  Problem_instance(int m, int n, int k);

  // PUBLIC INTERFACE: getters
 public:
  inline int M() const {return a.size();}		
  inline int N() const {return a[0].size();}
  inline int K() const {return w.size();} // the capital letter is a workaround the fact that no index m can be defined when a function m() exists...

  //! returns the number of points covered by cluster j
  int M(int j);

  //! gets the cluster j(i) covering point i
  int J(int i);

  //! returns the m_j times n matrix containing the points covered by cluster j
  LaGenMatDouble A(int j);

  //! returns A(j) with the last column filled with ones, for affine regression problems
  LaGenMatDouble A_aff_reg(int j);

  //! returns A(j) with the last column filled with zeroes (equivalent to zero gamma for min norm LS problems)
  LaGenMatDouble A_lin_reg(int j);

  //! returns the b vector containing the last entry of each datapoint (useful for only for REGRESSION problems)
  LaVectorDouble b_reg(int j);

  //! returns the centroid of cluster j
  LaVectorDouble c(int j);

  // PUBLIC INTERFACE: setters
 public:
  //! updates a single i-j pair
  void update_distances(const int i, const int j);

  //! updated the distances from all points to a specified cluster
  void update_distances(const int j);

  //! updated the distances from all points to all clusters
  void update_distances();

  // PUBLIC INTERFACE: computational setters
 public:

  //returns the index of the cluster the domain of which contains the point according to RLP
  int index_of_RLP_induced_cluster( const int i );

  // creates a new random solution (in w, gamma, and x) and a cleans the Tabu List (if used)
  void create_random_solution_and_clean_TS( base_generator_type & generator );

  //! randomly populates the (already created!!) set of clusters
  //void randomly_populate_clusters(int * ext_pred = NULL);
  void randomly_populate_clusters( base_generator_type & generator );
	
  // PUBLIC INTERFACE: printers
 public:
  //! outputs a matlab-readable output JUST FOR TWO dimensional (n=2) problems
  void output_solution_to_matlab(string & file_name, bool flag_print_lines);

  void output_solution_to_file(string & file_name);

  //old stuff
  void print_binary_representation(string & file_path);

  void print_dat_file_IS_representation(string & file_path);

  void print_dat_file_SNOPT_representation(string & file_path);


  void print_dat_file_assignment(string & instance_file_name);
  void print_dat_file_distances(string &instance_file_name);

  // PUBLIC INTERFACE: rep invariant
  bool assert_correctness ();  

  // GENERAL FORMULAS
 public:
  //! caution: point_measure ONLY returns the distance stores in the distances matrix!
  double point_measure(int i, int j) const;
  double cluster_measure(int j) const;
  double solution_measure() const;

  // PUBLIC INTERFACE: various utilities
 public:
  //! find min and max coordinates in the set of points
  void get_min_max_dimensions(vector<double> & v_min, vector<double> & v_max) const;

  //! normalizes the dataset so as to be contained into the [0,1]^n box
  void normalize(double factor = 1);

  void cutoffprecision(double precision = 1000);

  // INTERNAL PRINTING METHODS
 public:
  //! prints M =, N=, K= and the dataset A.
  void output_standard_ampl_data(string & file_name);

  // FRIEND PRINT METHODS
 public:
  friend std::ostream& operator<< (ostream& os, Problem_instance& ps);
};



#endif
