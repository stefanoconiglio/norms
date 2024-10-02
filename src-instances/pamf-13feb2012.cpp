#include "pamf.h"

#define ILOUSESTL
#include <ilcplex/ilocplex.h>



void
pamf::Linear1normRegressionParameterUpdate( Problem_instance & pi, int j ) {

  IloInt m_j = IloInt( pi.M( j ) );
  LaGenMatDouble A_j = pi.A( j );
  IloInt n = IloInt( pi.N() );

  //cout << "CPLEX-BEGIN void exact::Linear1normRegressionParameterUpdate( Problem_instance & pi, int j ) " << "for j = " << j << ", with m_j = " << m_j << endl;
  
  assert (m_j > 0);

  IloEnv env;

  IloModel regressionModel = IloModel( env );
  IloCplex regression( regressionModel );
  regression.setOut(env.getNullStream()); //seems to be not working... 
  regression.setError(env.getNullStream()); //seems to be not working... 
  regression.setWarning(env.getNullStream()); //seems to be not working... 
  regression.setParam ( IloCplex::Threads, 1 );
 

  IloNumVarArray w( env, n, -IloInfinity, IloInfinity, IloNumVar::Float );
  IloNumVar gamma( env, -IloInfinity, IloInfinity, IloNumVar::Float );
  IloNumVarArray d( env, m_j, 0, IloInfinity, IloNumVar::Float );
  
  IloNumArray w_sol( env, n);
  IloNum gamma_sol;

  IloObjective obj;
  obj = IloAdd ( regressionModel, IloMinimize( env, IloSum( d ) ) );
  //regression = IloCplex ( regressionModel );
  
  for ( int i = 0; i < m_j; i++ ) {
    IloNumArray a( env, n );
    for ( int l = 0; l < n; l++)
      a[l] = A_j(i, l);
    regressionModel.add ( d[i] - IloScalProd(a, w) + gamma >= 0);
    regressionModel.add ( d[i] + IloScalProd(a, w) - gamma >= 0);
  }
  regressionModel.add( w[n-1] == -1 );
		   
  regression.extract( regressionModel );
  //cout << regressionModel << endl;
  //regression.exportModel( "/home/aa23/eclipse_workspace/hyperplane_clustering/bin/pluto.lp" );
  regression.solve();
  
  gamma_sol = regression.getValue ( gamma );  
  
  for ( int l = 0; l < n; l++ ) {
    if ( regression.isExtracted( w[l] ) )
      w_sol[l] = regression.getValue( w[l]);
    // if all the points have a zero component l, the corresponding w[l] variable will not be extracted in the model, leading to an error when reading its value after solving
    else
      w_sol[l] = 0;
    pi.w[j](l) = w_sol[l];
  }
  pi.gamma[j] = gamma_sol;

  //cout << "CPLEX-END void exact::Linear1normRegressionParameterUpdate( Problem_instance & pi, int j ) " << "for j = " << j << ", with m_j = " << m_j << endl;

  env.end();
}



double
pamf::plain_RLP_parameter_update(Problem_instance & pi )  {
  
  //cout << "CPLEX-BEGIN vector<pair<int,int> > pamf::perform_RLP(Problem_instance & pi)" << endl;

  //! general variables
  vector<pair<int,int> >	reassignments;

  IloInt m = IloInt( pi.M() );
  IloInt k = IloInt( pi.K() );
  IloInt n = IloInt (pi.N() );		//! we are working on the independent dimensions ONLY

  IloEnv env;

  IloModel RLPModel = IloModel( env );
  IloCplex RLP( RLPModel );
  RLP.setOut(env.getNullStream()); //seems to be not working... 
  RLP.setError(env.getNullStream()); //seems to be not working... 
  RLP.setWarning(env.getNullStream()); //seems to be not working... 
  RLP.setParam ( IloCplex::Threads, 1 );

  typedef IloArray<IloNumVarArray> IloNumVarArrayArray;

  IloNumVarArrayArray W( env, k ); 
  for ( IloInt j = 0; j < k; j++) {
    W[j] = IloNumVarArray ( env, n-1, -IloInfinity, IloInfinity, IloNumVar::Float );
  }
  
  IloNumVarArray G( env, k, -IloInfinity, IloInfinity, IloNumVar::Float );
  IloNumVarArray e( env, m, 0, IloInfinity, IloNumVar::Float );
  
  IloObjective obj;
  obj = IloAdd ( RLPModel, IloMinimize( env, IloSum( e ) ) );
  //RLP = IloCplex ( RLPModel );


  for ( IloInt i = 0; i < m; i++ ) {
    IloNumArray a( env, n-1 );
    for ( int l = 0; l < n-1; l++)
      a[l] = pi.a[i](l);
    
    for ( IloInt j = 0; j < k; j++ ) {
      if ( pi.assignment[i][j] == false )
	continue; //! no constraint if the ith point is not assigned to the jth cluster
      for ( IloInt jj = 0; jj < k; jj++ ) {
	if ( j == jj )
	  continue;
	RLPModel.add ( e[i] + IloScalProd(a, W[j])  - IloScalProd(a, W[jj])  - G[j] + G[jj] >= 1 );
      }
    }
  }
  
  RLP.extract( RLPModel );
  // cout << RLPModel << endl;
  RLP.solve();

  for ( IloInt j = 0; j < k; j++ ) {
    for ( IloInt l = 0; l < n-1; l++ ) {
      if ( RLP.isExtracted( W[j][l] ) )
	pi.W[j](l) = RLP.getValue( W[j][l] );
      else
	pi.W[j](l) = 0;
    }
    pi.G[j] = RLP.getValue( G[j] );
  }

  // for ( IloInt i = 0; i < m; i++ )
  //   if ( RLP.isExtracted( e[i] ) )
  //     cout << "e[ " << i << "]" << RLP.getValue( e[i] ) << endl;

  IloNum ret = RLP.getValue( obj );

  env.end();

  return double( ret );

}



// vector<pair<int,int> >
// pamf::perform_RLP(Problem_instance & pi, bool doRLPonly )  {
  
//   //cout << "CPLEX-BEGIN vector<pair<int,int> > pamf::perform_RLP(Problem_instance & pi)" << endl;

//   //! general variables
//   vector<pair<int,int> >	reassignments;

//   IloInt m = IloInt( pi.M() );
//   IloInt k = IloInt( pi.K() );
//   IloInt n = IloInt (pi.N() );		//! we are working on the independent dimensions ONLY

//   IloEnv env;

//   IloModel RLPModel = IloModel( env );
//   IloCplex RLP( RLPModel );
//   RLP.setOut(env.getNullStream()); //seems to be not working... 
//   RLP.setError(env.getNullStream()); //seems to be not working... 
//   RLP.setWarning(env.getNullStream()); //seems to be not working... 
//   RLP.setParam ( IloCplex::Threads, 1 );

//   typedef IloArray<IloNumVarArray> IloNumVarArrayArray;

//   IloNumVarArrayArray W( env, k ); 
//   for ( IloInt j = 0; j < k; j++) {
//     W[j] = IloNumVarArray ( env, n-1, -IloInfinity, IloInfinity, IloNumVar::Float );
//   }
  
//   IloNumVarArray G( env, k, -IloInfinity, IloInfinity, IloNumVar::Float );
//   IloNumVarArray e( env, m, 0, IloInfinity, IloNumVar::Float );
  
//   IloObjective obj;
//   obj = IloAdd ( RLPModel, IloMinimize( env, IloSum( e ) ) );
//   //RLP = IloCplex ( RLPModel );


//   for ( IloInt i = 0; i < m; i++ ) {
//     IloNumArray a( env, n-1 );
//     for ( int l = 0; l < n-1; l++)
//       a[l] = pi.a[i](l);
    
//     for ( IloInt j = 0; j < k; j++ ) {
//       if ( pi.assignment[i][j] == false )
// 	continue; //! no constraint if the ith point is not assigned to the jth cluster
//       for ( IloInt jj = 0; jj < k; jj++ ) {
// 	if ( j == jj )
// 	  continue;
// 	RLPModel.add ( e[i] + IloScalProd(a, W[j])  - IloScalProd(a, W[jj])  - G[j] + G[jj] >= 1 );
//       }
//     }
//   }
  
//   RLP.extract( RLPModel );
//   // cout << RLPModel << endl;
//   RLP.solve();

//   for ( IloInt j = 0; j < k; j++ ) {
//     for ( IloInt l = 0; l < n-1; l++ ) {
//       if ( RLP.isExtracted( W[j][l] ) )
// 	pi.W[j](l) = RLP.getValue( W[j][l] );
//       else
// 	pi.W[j](l) = 0;
//     }
//     pi.G[j] = RLP.getValue( G[j] );
//   }

//   // for ( IloInt i = 0; i < m; i++ )
//   //   if ( RLP.isExtracted( e[i] ) )
//   //     cout << "e[ " << i << "]" << RLP.getValue( e[i] ) << endl;




//   //--------------------
//   //! compile the reassignments vector
//   //--------------------

//   // if ( RLP.getValue( obj ) > 0 ) {	//nothing has to be reassigned if nothing is misclassified
//   if ( !doRLPonly && RLP.getValue( obj ) > -1 ) {	//always do it

 
//     double	  max;
//     double	  curr; 
//     int		  max_cluster = -1;
//     pair<int,int> tmp_pair;
    
//     for ( int i = 0; i < m; i++) {

//       // double value_of_current_cluster = pi.W[pi.J(i)] * pi.a[i] + pi.a[i](n-1) - pi.G[pi.J(i)];
//       // max = -1e300;
//       // for ( int j = 0; j < k; j++) {
//       // 	curr = pi.W[j] * pi.a[i] + pi.a[i](n-1) - pi.G[j]; // the last component of W is set to -1 (the last component of 'a' is the dependent variable y, which is not involved in the partitioning). In the scalar product, it is summed up multiplied by (-1); adding the last component again we remove it.
//       // 	// cout << "curr " << curr << " of j = " << j << endl;
//       // 	if (curr >= max - 1e-06) {
//       // 	  max = curr;
//       // 	  max_cluster = j;
//       // 	}
//       // }
//       // assert (max > -1e300 && max_cluster != -1);

//       // int max_cluster2 = pi.index_of_RLP_induced_cluster( i );
//       // cout << "max_cluster " << max_cluster << endl;
//       // cout << "max_cluster2 " << max_cluster2 << endl;
//       // assert (max_cluster == max_cluster2);

//       int max_cluster = pi.index_of_RLP_induced_cluster( i );
      

//       // if (pi.assignment[i][max_cluster] == false &&
//       // 	  max >= value_of_current_cluster ) { 
//       // 	tmp_pair.first = i;
//       // 	tmp_pair.second = max_cluster;
//       // 	reassignments.push_back(tmp_pair);
//       // }
//       tmp_pair.first = i;
//       tmp_pair.second = max_cluster;
//       reassignments.push_back(tmp_pair);

//       //Note that ">=" is needed both in
//       // "if (curr >= max) {"
//       // and in
//       // "	  max >= value_of_current_cluster ) { "
//       // Assume ">". Pick srncr_100_2_2_10022 problem and run "./k-hc-new -f ../data/dam2011/pool1/medium/srncr_100_2_2_10022_ins -a maPR -P affine_regression -m 10 -T 0.1 -O 1 -D 0.99 -I 0.99 -t 2:aspire -C amm_distances -c closest -M -w -F" on it. You will get a solution where two points with the same first component are assigned to two different clusters. The RLP will yield a separating hyperplane containing both of them. When checking separability, we get a solution:
//       // Variable Name           Solution Value
//       // id42                          0.400000
//       // id52                          1.600000
//       // id2                           4.000000
//       // id4                          13.800000
//       // All other variables in the range 1-105 are 0.
      
//       // The constraints for nonzero misclassification error are:
//       // id182: id42 - 3.3 id2 + 3.3 id1 + id4 - id3 >= 1
//       // id202: id52 + 3.3 id2 - 3.3 id1 - id4 + id3 >= 1
//       // The first two vars are misclassification errors, the other two are W_1 and G_1 (W_2 and G_2 are zero due to k=2).
//       // For those W_1 and G_1 we get id42 >= 0.4 from id183 and id52 >= 1.6 from id202. Their sum is 2: nonzero.
//       //Putting a ">=" solves the problem. With ">", a point is reassigned only if "max > value_of_current_cluster". Therefore, the two points with the same first component are NOT reassigned to the same cluster. With a ">=" they are: when there are ties, the cluster with largest index is always chosen.
//     }
//   }

//   // cout << "CPLEX-END vector<pair<int,int> > pamf::perform_RLP(Problem_instance & pi)" << endl;

//   env.end();

//   return reassignments;
// }


// bool
// pamf::reassign_with_RLP(Problem_instance & pi) {
  
//   bool solution_changed = false;
  
//   vector<pair<int,int> > reassignments = pamf::perform_RLP(pi);
//   // cout << pi.assignment << endl;
//   // cout << "reassignments " << endl;
//   // for (int i = 0; i < reassignments.size(); i++) {
//   //   cout << reassignments[i].first << " " << reassignments[i].second << endl;
//   // }

//   if (reassignments.size() != 0) {
//     //! reassign missclassified points to their 'closest' cluster
//     for (int i = 0; i < int(reassignments.size()); i++) {
//       int pt = reassignments[i].first;
//       int old_cl = pi.J(pt);
      
//       int cl = reassignments[i].second;
//       pi.assignment[pt][old_cl] = false;
//       pi.assignment[pt][cl] = true;
//     }
//     solution_changed = true;
//   }
//   // cout << pi.assignment << endl;
  
//   return solution_changed;
// }

void
pamf::RLP_based_combinatorial_reassignment( Problem_instance & pi, double time, vector<bool> * globalIllAssigned ) {

  // cout << "CPLEX-BEGIN void pamf::RLP_based_combinatorial_reassignment( Problem_instance & pi, double time ) " << endl;

  //! general variables
  vector<pair<int,int> >	reassignments;

  IloInt m = IloInt( pi.M() );
  IloInt k = IloInt( pi.K() );
  IloInt n = IloInt (pi.N() );		//! we are working on the independent dimensions ONLY

  IloEnv env;

  IloModel reassignerModel = IloModel( env );
  IloCplex reassigner( reassignerModel );
  // reassigner.setOut(env.getNullStream());
  // reassigner.setError(env.getNullStream());
  // reassigner.setWarning(env.getNullStream());

  // cout << "TIME " << time << endl;

  if ( time < 1e300 )
    reassigner.setParam ( IloCplex::TiLim, (time > 0.01 ? time : 0) );
  reassigner.setParam ( IloCplex::Threads, 1 );
  reassigner.setParam ( IloCplex::ClockType, 1 );
  reassigner.setParam ( IloCplex::MIPEmphasis, 4 );
  reassigner.setParam ( IloCplex::IntSolLim, pi.ops.RLP_based_combinatorial_reassignment_sol_limit );


  typedef IloArray<IloNumVarArray> IloNumVarArrayArray;

  IloNumVarArrayArray W( env, k ); 
  for ( IloInt j = 0; j < k; j++) {
    W[j] = IloNumVarArray ( env, n-1, -IloInfinity, IloInfinity, IloNumVar::Float );
  }
  
  IloNumVarArray G( env, k, -IloInfinity, IloInfinity, IloNumVar::Float );
  // IloNumVarArray e( env, m, 0, IloInfinity, IloNumVar::Float );

  typedef IloArray<IloNumVarArray> IloNumVarArrayArray;
  IloNumVarArrayArray x( env, m ); 
  for ( IloInt i = 0; i < m; i++) {
    x[i] = IloNumVarArray ( env, k, 0, 1, IloNumVar::Bool );
  }

  IloExpr obj( env );

  for ( IloInt i = 0; i < m; i++) {
    // if ( globalIllAssigned != NULL )
    //   if ( !reassignIllAssigned && (*globalIllAssigned)[i] == true )
    // 	continue;
    //   else if ( reassignIllAssigned && (*globalIllAssigned)[i] == false )
    // 	continue;
    for ( IloInt j = 0; j < k; j++) {
         obj += pi.distance[i][j] * x[i][j];
    }
  }
  reassignerModel.add( IloMinimize( env, obj ) );

  //reassigner = IloCplex ( reassignerModel );


  for ( IloInt i = 0; i < m; i++ ) {
    // if ( globalIllAssigned != NULL )
    //   if ( !reassignIllAssigned && (*globalIllAssigned)[i] == true )
    // 	continue;
    //   else if ( reassignIllAssigned && (*globalIllAssigned)[i] == false )
    // 	continue;

    reassignerModel.add ( IloSum(x[i]) == 1 );	

    //ill-assigned points are NOT reassigned to their current cluster
    if ( globalIllAssigned != NULL ) {
      // if ( reassignIllAssigned && (*globalIllAssigned)[i] == true ) {
      if ( (*globalIllAssigned)[i] == true ) {
	IloInt jOfI = IloInt( pi.J( i ) );
	reassignerModel.add ( x[i][jOfI]  == 0 );	
      }
    }

    IloNumArray a( env, n-1 );
    for ( int l = 0; l < n-1; l++)
      a[l] = pi.a[i](l);
    
    for ( IloInt j = 0; j < k; j++ ) {
      // if ( pi.assignment[i][j] == false )
      // 	continue; //! no constraint if the ith point is not assigned to the jth cluster
      for ( IloInt jj = 0; jj < k; jj++ ) {
	if ( j == jj )
	  continue;
	//ampl_model << "sum{l in N0} ( - a[" << i+1 << ",l] * (W[" << j+1 << ",l]-W[" << jj+1 << ",l]) ) + (G[" << j+1 << "]-G[" << jj+1 << "]) + 1  <= e[" << i+1<< "];" << endl;
	// reassignerModel.add ( e[i] + IloScalProd(a, W[j])  - IloScalProd(a, W[jj])  - G[j] + G[jj] + 9999 * ( 1 - x[i][j]) >= 1);	
	reassignerModel.add ( IloScalProd(a, W[j])  - IloScalProd(a, W[jj])  - G[j] + G[jj] + 9999 * ( 1 - x[i][j]) >= 1);	
	// reassignerModel.add ( IloScalProd(a, W[j])  - IloScalProd(a, W[jj])  - G[j] + G[jj] + 999 * ( 1 - x[i][j]) >= 1);	
      }
    }
  }
  
  reassigner.extract( reassignerModel );
  //cout << reassignerModel << endl;

  try {
    cout << "Before issuing the solve command, with timelimit = " << time << endl;
    reassigner.solve();
    cout << "Solver status = " << reassigner.getStatus() << endl;
    cout << "Solver status = " << reassigner.getCplexStatus() << endl;
  }
  catch (  IloCplex::Exception e) {
    cout << e << endl;
  }

  //if ( reassigner.getStatus() == IloCplex::AbortTimeLim ) { //at timelimit, it exits with status = Unknown...
  //if ( reassigner.getStatus() != IloCplex::Optimal ) {
  if ( reassigner.getCplexStatus() != IloCplex::Optimal && reassigner.getCplexStatus() != IloCplex::SolLim) {
    cout << " Cplex returned without either status IloCplex::Optimal or IloCplex::SolLim (possibly due to timeout): terminating" << endl;
    env.end();
    return;
  } else {


    for ( IloInt j = 0; j < k; j++ ) {
      for ( IloInt l = 0; l < n-1; l++ ) {
	if ( reassigner.isExtracted( W[j][l] ) )
	  pi.W[j](l) = reassigner.getValue( W[j][l] );
	else
  	pi.W[j](l) = 0;
      }
      pi.G[j] = reassigner.getValue( G[j] );
    }

    for ( IloInt i = 0; i < m; i++) {
      // if ( globalIllAssigned != NULL  )
      //   if ( !reassignIllAssigned && (*globalIllAssigned)[i] == true )
      // 	continue;
      //   else if ( reassignIllAssigned && (*globalIllAssigned)[i] == false )
      // 	continue;
      
      for ( IloInt j = 0; j < k; j++) {
	if ( reassigner.isExtracted( x[i][j] ) )
	  pi.assignment[i][j] = reassigner.getValue( x[i][j] );
	else
	  pi.assignment[i][j] = 0;
      }
    }

#ifdef __DEBUG__
    double misclass_error = pamf::check_RLP ( pi );
    cout << "Assertion of separability: " << misclass_error << endl;
    assert ( misclass_error == 0 );
#endif

    // cout << "CPLEX-END void pamf::RLP_based_combinatorial_reassignment( Problem_instance & pi, double time ) " << endl;
    
  }

  env.end();

  return;
}



int
Problem_instance::index_of_RLP_induced_cluster( const int i ) {
 
  int k = this->K();
  int n = this->N();

  double	  max;
  double	  curr; 
  int		  max_cluster = -1;
  
  double value_of_current_cluster = this->W[this->J(i)] * this->a[i] + this->a[i](n-1) - this->G[this->J(i)];
  max = -1e300;
  for ( int j = 0; j < k; j++) {
    curr = this->W[j] * this->a[i] + this->a[i](n-1) - this->G[j]; // the last component of W is set to -1 (the last component of 'a' is the dependent variable y, which is not involved in the partitioning). In the scalar product, it is summed up multiplied by (-1); adding the last component again we remove it.
    if (curr >= max - 1e-06) { //small tolerance
      max = curr;
      max_cluster = j;
    }
  }
  // cout << "max j " << max_cluster << endl;
  assert (max > -1e300 && max_cluster != -1);
  
  return max_cluster;
  
}

void
algorithm::papavero_combinatorial_reassignment( Problem_instance & pi, base_generator_type & generator ) {

  // pamf::perform_RLP( pi, true );
  pamf::plain_RLP_parameter_update( pi );

  int k = pi.K();
  
  for (int i = 0; i < pi.M(); i++) {

    int old_j = pi.J(i);

    vector<double> probs( k, 0 );

 
    // cout << "old_j " << old_j << endl;

    int num_of_zero_prob_clusters = 0;
    for ( int j = 0; j < k; j++ ) {
      //we assume that point i belongs to the portion of the domain currently assignd to cluster j
      double RLP_error =  0;
      for ( int jj = 0; jj < k; jj++ ) {
	double curr_RPL_error = pi.a[i] * ( pi.W[j] - pi.W[jj] ) - (pi.G[j] - pi.G[jj]) + 1;
	if ( curr_RPL_error > RLP_error ) 
	  RLP_error = curr_RPL_error;
      }

      // cout << " j " << j << " RLP error " << RLP_error << endl;

      // if ( j == old_j ||
      if  ( pi.ops.TS.flag_use_TS && pi.tabu_list.is_tabu(i,j)
	    && pi.tabu_list.does_aspire(i,j,pi.distance[i][j]) == false ) {
        probs[j] = 0;
	num_of_zero_prob_clusters++;
      }
      else
	probs[j] = 1.0/(1 + pi.distance[i][j]) * (1 + RLP_error );
    }

    // cout << "probs " << probs << endl;

    int new_j;
    if ( num_of_zero_prob_clusters == 0 ) {
      new_j = random_generator::integer( k, generator ); //a sort of implicit shaker is applied if there is no possible reassignment
      // cout << "All-zero probs: uniformly sampling" << endl;
    }
    else {
      utility::normalize_1norm( probs );
      // cout << "probs " << probs << endl;
      new_j = random_generator::sample_discrete_distribution( probs, generator );
    }

    // cout << new_j << endl;

    pi.assignment[i][old_j] = false;
    pi.assignment[i][new_j] = true;
    //cout << "[" << i << "," << old_j << "]-->[" << i << "," << new_j << "]" << endl;

    if ( pi.ops.TS.flag_use_TS && pi.distance[i][old_j] < pi.distance[i][new_j] )
      pi.tabu_list.add_old_and_remove_new( i, old_j, pi.distance[i][old_j], new_j);   
  }
  return;
}


void
pamf::plain_RLP_induced_combinatorial_reassignment( Problem_instance & pi ) {

  int m = pi.M();

  double	max;
  double	curr; 
  int		max_cluster = -1;
  pair<int,int> tmp_pair;
  
  for ( int i = 0; i < m; i++) {
      
    int max_cluster = pi.index_of_RLP_induced_cluster( i );

    pi.assignment[i][pi.J( i )] = false;
    pi.assignment[i][max_cluster] = true;
  }
}


//circular buffer
#include <boost/circular_buffer.hpp>

Return_value
algorithm::TSpamf(Problem_instance & pi, bool warmStar, double maxTime, base_generator_type & generator ) {

  pi.ops.TS.flag_use_TS = false; //we do not use the regular tabu search here.

  boost::circular_buffer<pair<vector<int>, int> > tabu_list( pi.ops.TS.length );

  const int m = pi.M();
  const int n = pi.N();
  const int k = pi.K();
  assert(k > 1);

  Chronometer timer;
  timer.start();

  if ( warmStar == false)
    pi.create_random_solution_and_clean_TS ( generator );

  cout << "LogInitTS_pamf\t"
      << " initObj "      << pi.solution_measure()
      << endl;

  Problem_instance best_pi = pi;

  cout << "TS_pamf:"
      << " conservative = " << pi.ops.APR.conservative
      << " samples = "<< pi.ops.TS.neighborhood_samples
      << " swaps = " << pi.ops.TS.swaps
      << endl;


  int iteration = 0;
  while ( timer.elapsed_time() < maxTime - TIMETOL ) {

    // if ( tabu_list.empty() == false )
    //   tabu_list.pop_front();

    // cout << "TS PRINT" << endl;
    // for ( boost::circular_buffer<pair< vector<int>, int> >::iterator it = tabu_list.begin(); it != tabu_list.end(); it++ ) {
    //   cout <<  it->first << "; " << it->second << endl;
    // }

    vector<double> probs = ill_assigned_probabilities( pi, generator );

    Problem_instance best_neighboring_pi = pi;
    Problem_instance prev_pi = pi;

    cout << "LogIntTS_pamf\t"
        << " origNeighSol "      << pi.solution_measure()
        << endl;
    
    bool tabu = false;
    int number_of_tabu = 0;

    double best_neighboring_pi_obj = 1e300;

    //generate neighborhood
    for ( int s = 0; ( timer.elapsed_time() < maxTime - TIMETOL ) && 
	    (s < pi.ops.TS.neighborhood_samples || tabu == true ); s++) {  //build a new sol until a nontabu one is found
      
      if ( number_of_tabu >= 10 ) {
	cout << "10 tabu solutions in a row found; applying shaker and going to the next iteration" << endl;
	TSpamf_shaker( pi, generator );
	break;
      }	

      cout << "Neighboring solution # " << s << endl;
      cout << "Previous was tabu? " << tabu << endl;

      cout << "Creating neighboring solution (new_pi) without Intertwining RLP " << endl;

      Problem_instance new_pi = pts_neighboring_solution( pi, probs, generator );

      if ( pi.ops.TS.cooper == false ) {
	cout << "new_pi obj value before Intertwining RLP = " << new_pi.solution_measure() << endl;
	
	// cout << "Intertwining RLP" << endl;
	pamf::plain_RLP_parameter_update ( new_pi );
	pamf::plain_RLP_induced_combinatorial_reassignment( new_pi );
	plain_parameter_update( new_pi );
	new_pi.update_distances();
	
	cout << "new_pi obj value after Intertwining RLP = " << new_pi.solution_measure() << endl;
      }
      else
	cout << "No Intertwining RLP needed, due to using \"Cooper\"" << endl;

     
      double new_pi_obj = new_pi.solution_measure();

      if ( new_pi_obj < best_neighboring_pi_obj ) { //store best neighboring sol
	best_neighboring_pi = new_pi;
	best_neighboring_pi_obj = new_pi_obj;
      }
      
      cout << "LogIntTSpamf-interNeighborhood\t"
	   << " newNeig "      << new_pi_obj
	   << " bestNeig "     << best_neighboring_pi_obj
	   << endl;

      vector< int > solution( m, -1 );
      for ( int i = 0; i < m; i++ )
	solution[i] = new_pi.J ( i );

      
      // if ( tabu_list.find( solution ) == tabu_list.end() )
      // 	tabu = false;

      tabu = false;
      for ( boost::circular_buffer<pair<vector<int>, int> >::iterator it = tabu_list.begin(); it != tabu_list.end(); it++ ) {
      	// if ( it->first == solution && it->second >= iteration-pi.ops.TS.length) {
      	if ( it->first == solution && it->second >= iteration-pi.ops.TS.length ) {
	  //NOTE: I checked the correctness of vector::operator== on the c++ distro installed on MY MACHINE, and it is ok
	  // if ( pi.tabu_list.use_aspiration_criterion() == true ) {
	  //   if ( best_neighboring_pi_obj < best_pi.solution_measure () ) {
	  //     cout << "Aspiration kicked in: tabu status overridded" << endl;
	  //     continue;
	  //   }
	  // }
      	  tabu = true;
	  number_of_tabu++;
	  cout << "The current solution is tabu (" number_of_tabu << "-th in a row)" << endl;
	  // cout << "Current solution " << solution << endl;
	  // cout << "Tabu solution " << it->first << endl;
      	  break;
      	}
	// tabu = false;
      }
      if ( tabu == false )
	cout << "The current solution is NOT tabu " << endl;       
      
    } //neighborhood explored
    

    cout << "LogIntTSpamf-exitNeighborhood\t"
	 << " bestPi "      << best_pi.solution_measure() 
	 << " bestNeig "    << best_neighboring_pi.solution_measure()
	 << endl;
    
    if ( best_neighboring_pi.solution_measure() < best_pi.solution_measure() ) { //store best
      best_pi = best_neighboring_pi;
      cout << "New best sol found" << endl;
    }

    bool worsening = false;
    if ( best_neighboring_pi.solution_measure() > prev_pi.solution_measure() ) { //worsening solution found
      cout << "Best neigh sol is not better than the previous one: adding solution to tabu list " << endl;
      worsening = true;
      vector< int > solution( m, -1 );
      for ( int i = 0; i < m; i++ )
	solution[i] = best_neighboring_pi.J ( i );
      
      pair<vector<int>, int> new_pair( solution, iteration );
      tabu_list.push_back( new_pair );

      // cout << "TS PRINT" << endl;
      // for ( boost::circular_buffer<pair<vector<int>, int> >::iterator it = tabu_list.begin(); it != tabu_list.end(); it++ )
      // 	cout <<  it->first << endl;
 
    }

    pi = best_neighboring_pi; //next solution is build AROUND the current one...
    prev_pi = pi;

    iteration++;

    //cout << "It: " << iteration << "\t" << obj << endl;
    cout << "LogCurrTS_pamf\t"
        << " currObj "      << pi.solution_measure()
        << " worsening "    << worsening
        << " bestObj "      << best_pi.solution_measure()
        << " iterations "   << iteration
        << " timeLeft "     << maxTime - timer.elapsed_time()
        << endl;

    // pi.tabu_list.full_print_tabu_list();

  } // end of the iterations

  cout << "LogFinalTS_pamf\t"
      // << " currObj "      << obj
      << " currObj "      << pi.solution_measure()
      << " bestObj "      << best_pi.solution_measure()
      << " iterations "   << iteration
      << endl;


  Return_value retval;
  retval.obj = best_pi.solution_measure();
  retval.time_to_best = timer.elapsed_time();
  retval.total_time = retval.time_to_best; //the same (last sol = best sol)
  retval.iterations_to_best = iteration;
  retval.total_iterations = retval.iterations_to_best;  //the same (last sol = best sol)

  return retval;
}

void
algorithm::TSpamf_random_shaker(Problem_instance & pi, base_generator_type & generator ) {

  int m = pi.M();
  int k = pi.K();

  // vector<double> probs = ill_assigned_probabilities( pi ); //12-09-2011
  vector<double> probs = ill_assigned_probabilities( pi, generator );
  cout << "ill-assignment probs " << probs << endl;

  for ( int i = 0; i < m; i++) {

    if ( random_generator::random( generator ) <= probs[i] ) {

      cout << "Shaking point " << i << endl;

      int new_j = -1;

      if ( pi.ops.APR.reassignment_target == CLOSEST ) {
        double best_dist = 1e300;

        for (int jj = 0; jj < k; jj++) {
          if (jj == pi.J( i ) )
            continue;
          if (pi.distance[i][jj] < best_dist) {
            best_dist = pi.distance[i][jj];
            new_j = jj;
          }
        }
      }
      else if ( pi.ops.APR.reassignment_target == RAND_DIST ) {
        vector<double> cluster_probs( k, -1);

        for (int jj = 0; jj < k; jj++) {
          if (jj == pi.J( i ) )
            cluster_probs[jj] = 0;
          else
            cluster_probs[jj] = 1.0/pi.distance[i][jj]; //1 was amiss
        }

        utility::normalize_1norm( cluster_probs );

        // do {
        //   new_j = random_generator::sample_discrete_distribution( cluster_probs, generator );
        // } while ( new_j == pi.J( i ) );
        //find SEPTEMBER_ISSUE in this file to get a comment on the error that was fixed
        new_j = random_generator::sample_discrete_distribution( cluster_probs, generator );

      }
      else if ( pi.ops.APR.reassignment_target == RAND ) {
        new_j = random_generator::integer( k, generator );
      }

      assert ( new_j != -1 );

      // cout << i << "from " << pi.J(i) << " to " << new_j << endl;

      int old_j = pi.J(i);
      pi.assignment[i][pi.J(i)] = false;
      pi.assignment[i][new_j] = true;

      //! the GOOD assignment becomes tabu ....
      // cout << "TS before" << endl;
      // pi.tabu_list.full_print_tabu_list(); cout << endl;

      if ( old_j != new_j )
        pi.tabu_list.add( i, pi.J(i), pi.distance[i][pi.J(i)] );

      // cout << "TS after" << endl;
      // pi.tabu_list.full_print_tabu_list(); cout << endl;

    }
  }

  plain_parameter_update( pi );

  pi.update_distances( );
}
