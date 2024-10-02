#include"inclusions.h"

#ifndef __RANDOM_GENERATOR__
#define __RANDOM_GENERATOR__

using namespace std;

class random_generator {
 public:
  // returns a random number in the real interval [0,1).
  static double random( base_generator_type & generator );

  // returns a random number from a Gaussian distribution 
  static double nrandom( base_generator_type & generator );

  // returns a random integer in the interval [0, limit)
  static int integer(int limit, base_generator_type & generator );

  //! randomizes a LaVectorDouble by continous [0,1] random numbers, following a given chain (if specified), appending a random + or - sign before each component
  static void randomize(LaVectorDouble & vec, base_generator_type & generator );

  //! randomizes a LaVectorDouble by continous [0,1] random numbers, following a given chain (if specified)
  static void randomize01(LaVectorDouble & vec, base_generator_type & generator );

  //! randomizes a vector<double> by continous [0,1] random numbers, following a given chain (if specified)
  static void randomize01( vector<double> & vec, base_generator_type & generator );

  //! randomizes a LaVectorDouble by continous [0,1] random numbers, following a given chain (if specified)
  static void randomize01(LaGenMatDouble & matrix, base_generator_type & generator );

  static int sample_discrete_distribution( vector<double> & probs, base_generator_type & generator );

  //! returns a random number in the real interval (min,max).
  static double random_in_box(double min, double max, base_generator_type & generator );

  //! returns a vector contained the occurrencies of a samples by a cumulative distribution
  static vector<int> sampleByCumulative( vector<double> & probabilities, int number_of_samples, base_generator_type & generator );

};
#endif
