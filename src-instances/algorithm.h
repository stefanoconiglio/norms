#include "inclusions.h"


using namespace std;

#ifndef __ALGORITHM__
#define __ALGORITHM__

class algorithm
{

  //-----------------------
  //! COMPUTATIONAL SETTERS
  //-----------------------
 public:

  //REV1 
  //My implementation of Bemporad's criterion; does pi.ops.bemporad_passes iterations of the criterion; if pi.ops.bemporad_delta > 0, the method SHOULD be called as a PHASE; if pi.ops.bemporad_delta < 0, it SHOULD be called as a move; in any case, it uses "delta = pi.solution_measure()/m * abs (pi.ops.bemporad_delta)" for the parameter delta (names as in Bemporad's paper); if pi.ops.bemporad_c > 0, it considers the bemporad_c  2norm closest points for voting; if  pi.ops.bemporad_c < 0, it considers the bempoard_c fraction "c = int(m*abs(pi.ops.bemporad_c))"; in any cae, if c < 3, it is by force set to 3;
"

  static void
    bemporads_criterion(Problem_instance & pi);

  static Return_value
    papavero(Problem_instance & pi, bool warmStart, double maxTime, base_generator_type & generator );

  //CHECKED 25-09-2011
  // returns the ill-assigned value of point i
  static double
    ill_assigned_value ( int i, Problem_instance & pi, base_generator_type & generator );

  //CHECKED 25-09-2011
  //if cluster j is empty, randomly chooses a new set of hyperplane parameters for it, and updates the distances
  static void
    degenerancy_avoider( Problem_instance & pi, int j, base_generator_type & generator );

  //CHECKED 25-09-2011
  //returns a bool vector indicating which of the points assigned to j are ill-assigned
  static vector<bool>
    find_ill_assigned(Problem_instance & pi, int j, double alpha, base_generator_type & generator );

  //CHECKED 25-09-2011
  //Bradley and Mangasarian's algorithm; takes the tabu list into account, if used
  static Return_value
    BM(Problem_instance & pi, bool warmStart, double maxTime, base_generator_type & generator );
  
  //an anytime BM that stores the best solution found and makes tabu the whole set of moves that lead to a worsening one
  static Return_value
    TSBM(Problem_instance & pi, bool warmStar, double maxTime, base_generator_type & generator );

  //CHECKED 25-09-2011
  //reassigns each point to its closest nontabu cluster; it skips points with (*skip)[i] == true;
  static void
    plain_combinatorial_reassignment(Problem_instance & pi, vector<bool>* skip = NULL);

  static void
    papavero_combinatorial_reassignment( Problem_instance & pi, base_generator_type & generator );

  //CHECKED 25-09-2011
  //recomputes the hyperplane parameters according to the problem type: clustering, lin. regress, aff. regress, for cluster j
  static void
    plain_parameter_update(Problem_instance & pi, const int j );

  //CHECKED 25-09-2011
  //as before, cycling over all the clusters
  static void
    plain_parameter_update(Problem_instance & pi );

  //CHECKED 25-09-2011
  //returns the index of the closest, nontabu cluster for i
  static int
    find_closest_nontabu_cluster( Problem_instance & pi, int i );

  //CHECKED 25-09-2011
  //PR --COR 2011
  static Return_value
    PR( Problem_instance & pi, double alpha, bool warmStart, double maxTime, base_generator_type & generator );

  //CHECKED 25-09-2011
  //reassigns all the ill_ass points assigned to cluster j
  static int
    APR_combinatorial_reassignment( Problem_instance & pi, int j, vector<bool> & ill_ass, base_generator_type & generator );

  //CHECKED 25-09-2011
  //returns the index of a target cluster, possibly skipping the current one; when reassigning to the closest cluster, the current one is skipped, as well
  //as tabu reassignments; if no reassignment is possible, it returns -1. If reassigning randomly (both options), returning the current cluster or a tabu one
  //is a 0-measure event, thus being possible
  static int
    suggest_nontabu_reassignment( Problem_instance & pi, int i, base_generator_type & generator );

  //CHECKED 25-09-2011
  //first considers the closest cluster; if it coincides with the current one, it calls suggest_nontabu_reassignment
  static int
    suggest_conservative_nontabu_reassignment( Problem_instance & pi, int i, base_generator_type & generator );


  //CHECKED 25-09-2011
  //returns ill-assigned probabilities obtained by calling ill_assigned_value() for all points, then normalizing the vector
  static
    vector<double> ill_assigned_probabilities( Problem_instance & pi, base_generator_type & generator );

  //CHECKED 25-09-2011
  //returns true if no point-cluster pair is nontabu, false otherwise
  static bool
    everything_tabu( Problem_instance & pi );

  //CHECKED 26-09-2011
  //calls either PTS or TS, as detailed in the paper
  //note: it does not call APR_combinatorial_reassignment, because, here, we add the reassignments of ill-ass. points to the tabu list
  //only when they give rise to a worsening solution; APR_combinatorial_reassignment adds them regardless of the value of the new solution
  static Return_value
    TS(Problem_instance & pi, bool warmStar, double maxTime, base_generator_type & generator );

  // TS algorithm for pamf
  static Return_value
    TSpamf(Problem_instance & pi, bool warmStar, double maxTime, base_generator_type & generator );

  static Problem_instance
    pts_neighboring_solution_sampling_all( const Problem_instance & pi, base_generator_type & generator );


  static Return_value
    justRandom( Problem_instance & pi, bool warmstart, double maxTime, base_generator_type & generator );


  //CHECKED 26-09-2011
  //samples a solution in the PTS neighborhood, also refining it with BM if pi.ops.TS.cooper == true
  static Problem_instance
    pts_neighboring_solution( const Problem_instance & pi, vector<double> & probs, base_generator_type & generator );

  //CHECKED 26-09-2011
  //creates a solution by moving point i to cluster j
  static Problem_instance
    ts_exchange_neighboring_solution( const Problem_instance & pi, int i, int j );


  //CHECKED 26-09-2011
  //multistart embedded --COR 2011
  static Problem_instance
    multistarter( Problem_instance & pi, t_algorithm_choice algorithm, double maxTime, base_generator_type & generator );


  //CHECKED 26-09-2011
  //Needed for comparisons with APR and TS and PTS, which are anytime algorithms
  //anytime version of BM or PR (and PW and PAPAVERO): keeps on multistarting BM or PR (or PW) until maxTime is exceeded, returning the best solution found.
  //can only be run in its multistarted version (maBM or maPR), letting m=1 to do a single run; This is because the initial solution to this method MUST BE PROVIDED (it is not randomly generated at start)
  static Return_value
    anytimer( Problem_instance & pi, t_algorithm_choice algorithm, double maxTime, base_generator_type & generator );

  //CHECKED 26-09-2011
  //APR --COR 2011
  static Return_value APR( Problem_instance & pi, bool warmStart, double maxTime, base_generator_type & generator );


  //CHECKED 26-09-2011
  // solves the orthgonal regression problem where all the points are assigned to the same cluster
  static Return_value compute_trivial_one_cluster_solution(Problem_instance & pi );


  //UNUSED IN COR2011-------------------------------------------------------------------------------------------------------

  //! point-wise mangasarian move: relocates every point, but one by one
  static Return_value PW(Problem_instance & pi, bool warmStar, double maxTime, base_generator_type & generator );

  //shakes the solution via ill-assignment stuff
  static void
    random_shaker(Problem_instance & pi, base_generator_type & generator );

  static Return_value
    anytimeShakenBM( Problem_instance & pi, double maxTime, base_generator_type & generator );


};

class ordering_class
{
 public:
  virtual bool operator()(const int & a, const int & b) {return false;}
  virtual ~ordering_class(){}
};

class less_then : ordering_class
{
 public:
  less_then() {}
  bool operator()(const int & a, const int & b)
  {
    return (values[a] < values[b]);
  }

  void set_vectors(vector <double> & values)
  {
    this->values = values;
  }

 private:
  vector<double> values;
};


class greater_then : ordering_class
{
 public:
  greater_then() {}
  bool operator()(const int & a, const int & b)
  {
    return (values[a] > values[b]);
  }

  void set_vectors(vector <double> & values)
  {
    this->values = values;
  }

 private:
  vector<double> values;
};

 class pointwise_distance_order : ordering_class
 {
  public:
   pointwise_distance_order() {}
   bool operator()(const int& a, const int& b)
   {
     if (disabled[a])
       return true;
     else if (disabled[b])
       return false;
     else
       return (distances[a] < distances[b]);
   }

   void set_vectors(const vector <double> & distances , const vector<bool> & disabled)
   {
     this->distances = distances;
     this->disabled = disabled;
   }

  private:
   vector<double> distances;
   vector<bool> disabled;
 };

/* class max_pointwise_distance_order : ordering_class */
/* { */
/*  public: */
/*   max_pointwise_distance_order() {} */
/*   bool operator()(const int& a, const int& b) */
/*   { */
/*     if (disabled[a]) */
/*       return false; */
/*     else if (disabled[b]) */
/*       return true; */
/*     else */
/*       return (distances[a] < distances[b]); */
/*   } */

/*   void set_vectors(const vector <double>& distances , const vector<bool>& disabled) */
/*   { */
/*     this->distances = distances; */
/*     this->disabled = disabled; */
/*   } */

/*  private: */
/*   vector<double> distances; */
/*   vector<bool> disabled; */
/* }; */

#endif
