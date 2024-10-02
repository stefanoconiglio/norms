#include "inclusions.h"

#include <ext/hash_map>
using __gnu_cxx::hash_map;

using namespace std;

#ifndef __SOLUTION_TABU_LIST__
#define __SOLUTION_TABU_LIST__

class solution_tabu_list {
 private:
  hash_map<vector< int >, int> map;
  int length;
  int iteration;
  
 public:
  solution_tabu_list() {};		// empty constructor: needed because of solution_tabu_list tabu_list belonging to Problem_instance
  
  solution_tabu_list(int length, int m, int k, bool aspire = false);

  //CHECKED 26-09-2011
  //sets tabu status and aspiration distance
  void
    add(int m, int k);

  bool
    is_tabu( vector<int> solution ) const;
  
  void remove_tabu_pair(const int i, const int j);

  //call it as the beginning of the loop of your alg, so it is initialized to 0, rather than -1
  void operator++(int);
  
  void
    empty();
    
};

#endif

