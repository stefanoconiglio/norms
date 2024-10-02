#include "inclusions.h"

#ifndef __PAMF__
#define __PAMF__


namespace pamf {

  //assuming pi.W, pi.G are consistent (due to calling plain_RLP_parameter_update(Problem_instance & pi ) ), it reassigns each point i to the cluster the domain of which containts i
  void plain_RLP_induced_combinatorial_reassignment( Problem_instance & pi );

  //computers RLP parameters for the current partition, minimizing the misclassification error; returns the total misclassification error
  double plain_RLP_parameter_update(Problem_instance & pi );

    // computers the parameter update in linear l1 norm using cplex on cluster j
  void Linear1normRegressionParameterUpdate( Problem_instance & pi, int j );
  
  /* //! calls AMPL and CPLEX with mangasarian's RLP model, returning a set of missclassified points with their 'closest' cluster */
  /* vector<pair<int,int> > perform_RLP(Problem_instance & pi, bool doRLPonly = false); */
  
  /* //! reassignment of missclassified points with RLP */
  /* bool reassign_with_RLP(Problem_instance & pi); */
  
  // reassing the points to the hyperplanes (for given distances), guaranteeing that the assignment is linearly partitionable, so as to minimize the total sum of distances
  void RLP_based_combinatorial_reassignment( Problem_instance & pi, double time, vector<bool> * skip = NULL );

  
}

#endif
