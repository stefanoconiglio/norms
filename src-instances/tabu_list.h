
#include "inclusions.h"

using namespace std;

#ifndef __TABU_LIST__
#define __TABU_LIST__

class full_tabu_list {
 private:
  vector<vector<int> > table; 	//! representation: -1: element NOT in tabu list
  vector<vector<double> > tabu_value;	//! point-to-cluster distance between point i and cluster j when the pair (i,j) is added to the tabu list
  deque <pair<int, int> > list;
  int length;
  int iteration;
  bool aspire;
  
 public:
  full_tabu_list() {};		// empty constructor: needed because of full_tabu_list tabu_list belonging to Problem_instance
  
  full_tabu_list(int length, int m, int k, bool aspire = false);

  //CHECKED 26-09-2011
  //sets tabu status and aspiration distance
  void
    add(int m, int k, double value = -1);

  //CHECKED 26-09-2011
  //when changing (i, old_j) into (i, new_j), it adds (i, old_j) to the tabu list and removes (i, new_j) from it (if it was there)
  void
    add_old_and_remove_new( int i, int old_j, double old_distance, int new_j );

  bool
    is_tabu(const int i, const int j) const;
  
  bool
    inline use_aspiration_criterion() {return aspire;}

  //CHECKED 25-09-2011
  //returns true if either the pair (i,j) is not tabu, or if, though being tabu, their distance is strictly smaller than the one
  //they had when put in the tabu list
  bool
    does_aspire(const int i, const int j, double value) const;
  
  //returns false if reassigning point i to any cluster (different from the current one) is tabu, true otherwise
  bool
    can_be_reassigned(const int i, const vector<double> distances, int current_j );

  //CHECKED 26-09-2011
  // unsets the pair (i,j) as a tabu pair
  void remove_tabu_pair(const int i, const int j);

  //call it as the beginning of the loop of your alg, so it is initialized to 0, rather than -1
  void operator++(int);
  
  void
    empty();
    
  void
    full_print_tabu_list() const; 
};

#endif

