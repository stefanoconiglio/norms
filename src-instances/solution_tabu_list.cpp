#include "solution_tabu_list.h"

solution_tabu_list::solution_tabu_list(int length, int m, int k, bool aspire) {

  this->length = length;
  
  vector <int> dummy_vector(k, -1);
  table = vector<vector<int> >(m, dummy_vector);
  
  vector <double> dummy_vector_double(k, -1);
  tabu_value = vector<vector<double> >(m, dummy_vector_double);
  
  iteration = 0;
  this->aspire = aspire;
}


void
solution_tabu_list::operator++(int) {

  if ( iteration < std::numeric_limits<int>::max() )
    iteration++;
  else { //ugly shit: it RESETS the TS whenever the counter goes out of the blue
    iteration = 0;
    empty();
  }

  
}

void
solution_tabu_list::add(int i, int j, double value) {

  table[i][j] = iteration;
  
  if ( aspire ) {
    assert (value != -1);		//! when using the aspiration criterion a value must ALWAYS be provided
    tabu_value[i][j] = value;
  }
}

void
solution_tabu_list::add_old_and_remove_new( int i, int old_j, double old_distance, int new_j ) {

  // add the old assignment to the tabu list, which we want to prevent
  if ( use_aspiration_criterion() == false )
    add( i, old_j );
  else
    add( i, old_j, old_distance );

  // since we are performing a (i, new_j) assignment, we make it nontabu anymore (if it was tabu; it was not, the method just puts a harmless -1)
  remove_tabu_pair( i, new_j );
}

bool
solution_tabu_list::is_tabu(const int i, const int j) const {

  if ( table[i][j] == -1 )
    return false;
  else
    return ( iteration <= table[i][j] + length );
}

bool
solution_tabu_list::does_aspire(const int i, const int j, double value) const {

  if ( is_tabu( i, j ) == false )
    return true;

  return ( value < tabu_value[i][j] );
  
}

bool
solution_tabu_list::can_be_reassigned(const int i, const vector<double> distances, int current_j ) {
  
  int k = distances.size();
  
  for ( int j = 0; j < k; j++ ) {


    // cout << "--------------" << endl;

    // cout << "testing i = " << i << " " << "j = " << j << " dist = " << distances[j] << "-- stored = " << tabu_value[i][j] << endl;
    // cout << "assigned = " << (j == current_j) << endl;
    // cout << "tabu = " << is_tabu(i, j) << endl;
    // cout << "aspire = " << does_aspire(i, j, distances[j] ) << endl;

    // cout << "--------------" << endl;

    if ( j == current_j )
      continue;
    if ( is_tabu( i, j ) == false )
      return true;
    else  {
      if ( use_aspiration_criterion() == true )       
	if ( does_aspire( i, j, distances[j] ) ) 
	  return true;
    }
      
  }

  // cout << "Not reassignable" << endl;
  return false;


}

void
solution_tabu_list::remove_tabu_pair(const int i, const int j) {
  table[i][j] = -1;		//! the pair is no more tabu
  tabu_value[i][j] = -1;
}

