#include "random.h"

using namespace std;


int
random_generator::sample_discrete_distribution( vector<double> & probs, base_generator_type & generator ) {

  double rand = random( generator );

  double sum = 0;
  uint i;
  for ( i = 0; i < probs.size(); i++ ) {
    sum += probs[i];
    
    if ( rand <= sum )
      // break; //read below
      return i;
  }
  
  // return i; //if from "i=0" to "i=probs.size()-1" we don't meet success, we get "i = probs.size()", the for loop terminates and "return i" returns "probs.size()"; fixed
  return probs.size()-1;
}


double
random_generator::random( base_generator_type & generator ) {

  // Define a uniform random number distribution which produces "double"
  // values between 0 and 1 (0 inclusive, 1 exclusive).
  boost::uniform_real<> uni_dist(0,1);
  boost::variate_generator<base_generator_type&, boost::uniform_real<> > uni(generator, uni_dist);

  return uni( );
}

double
random_generator::nrandom( base_generator_type & generator ) {

  boost::normal_distribution<> nd(0.0, 1.0);
  boost::variate_generator<base_generator_type&, boost::normal_distribution<> > var_nor(generator, nd);

  return var_nor();
}

int
random_generator::integer( const int limit, base_generator_type & generator ) {

  typedef boost::uniform_int<> distribution_type;
  typedef boost::variate_generator<base_generator_type&, distribution_type> gen_type;
  gen_type die_gen( generator, distribution_type(0, limit-1) );

  return die_gen();
}


void 
random_generator::randomize01( vector<double> & vec, base_generator_type & generator ) {
  for (int i = 0; i < vec.size(); i++)
    vec[i] = random( generator );
}


void 
random_generator::randomize01(LaVectorDouble & vec, base_generator_type & generator ) {
  for (int i = 0; i < vec.size(); i++)
    vec(i) = random( generator );
}

void 
random_generator::randomize01(LaGenMatDouble & matrix, base_generator_type & generator)
{	
  for (int i = 0; i < matrix.rows(); i++)
    for (int j = 0; j < matrix.cols(); j++)
      matrix(i,j) = random( generator );
}


void 
random_generator::randomize(LaVectorDouble & vec, base_generator_type & generator) {
  randomize01(vec,  generator );
  
  for (int i = 0; i < vec.size(); i++)
    if (random_generator::integer(2,  generator ) == 0)
      vec(i) = -vec(i);
}


double
random_generator::random_in_box(double min, double max, base_generator_type & generator)
{
  return random_generator::random( generator ) * fabs(max-min) + min;
}

vector<int> 
random_generator::sampleByCumulative( vector<double> & probabilities, int number_of_samples, base_generator_type & generator ) {

  //compute cumulative distribution
  int N = probabilities.size();
  vector<int> samples(N, 0);
  vector<double> cumulative(N,0);
  cumulative[0] = probabilities[0];
  for (int i = 1; i < N; i++)
    cumulative[i] = cumulative[i-1] + probabilities[i];

  //sample
  int element_to_samples = -1;

  //! do the cumulative sampling!
  for (int n = 0; n < number_of_samples; n++) {
    //! throws a 'normalized dice' [0-1)
    double rand_normalized = random_generator::random( generator );
    //! find the element of the cumulative that has the maximum value smaller (or equal) to the rand_normalized
    for (element_to_samples = N-1; element_to_samples >= 0; element_to_samples--) {
      if (cumulative[element_to_samples] <= rand_normalized)
	break;
    }
    
    //! updates the samples counter
    assert(element_to_samples >= 0);
    (samples[element_to_samples])++;
  }
  //! returns the samples counter
  return samples;
}
