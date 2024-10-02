#include "instance_generator.h"

void
Instance_generator::populate(Problem_instance & pi, Instance_type & type, Instance_generation_printing_data & igpd)
{
  // adding the name of the instance to the igpd structure (used for printing purposes only)
  switch (type.data_type)
    {
    case RANDOM:
      igpd.name = "Random";
      break;
    case WAVE:
      igpd.name = "Wave";
      break;
    case RWAVE:
      igpd.name = "Random_wave";
      break;
    case SEMI_RANDOM:
      igpd.name = "Semi_random";
      break;
    case SEMI_RANDOM_CONTINUOS_REGRESSION:
      igpd.name = "Semi_random_continuous_regression";
      break;
    case SEMI_RANDOM_NON_CONTINUOS_REGRESSION:
      igpd.name = "Semi_random_non_continuous_regression";
      break;
    case SEMI_RANDOM_NON_CONTINUOS_REGRESSION_PANDA:
      igpd.name = "Semi_random_non_continuous_regression_panda";
      break;
    case PCLUSTERING_INSTANCE:
      igpd.name = "Semi_random_pclustering";
      pi.problem_type.problem = PCLUSTERING;
      break;
    default:
      cout << "Unrecognized parameter in populate" << endl;
      abort();
    }

  if (type.data_type == PCLUSTERING_INSTANCE)
    return populate_pclustering(pi, type, igpd);

  const int m = pi.M();
  const int n = pi.N();
  const int k = pi.K();

  LaVectorDouble population_probabilities(k);
  random_generator::randomize01(population_probabilities);
  utility::normalize_1norm(population_probabilities);


  if (type.data_type == RANDOM)
    {
      for (int i = 0; i < int(pi.a.size()); i++)
	random_generator::randomize01(pi.a[i]);
      return;
    }
  if (type.data_type == WAVE)
    {
      wave(pi);
      return;
    }
    
  if (type.data_type == RWAVE)
    {
      rwave(pi);
      return;
    }  
    
   if (type.data_type == SEMI_RANDOM_NON_CONTINUOS_REGRESSION)
    {
      populate_srnc(pi,type, igpd.file_path);
      return;
    }
    
   if (type.data_type == SEMI_RANDOM_CONTINUOS_REGRESSION)
    {
      populate_src(pi,type, igpd.file_path);
      return;
    }

  int rand_j = -1;
	
  //! madness preventer
  int max_good_point_shooting_trials = 1000000;

  //! guassian error variances (different per hyperplane)
  vector<double> variances;
  for (int kk = 0; kk < k; kk++)
    {
      double discount = random_generator::random();
      discount *= 0.3;
      discount += 0.7;
      variances.push_back(discount*type.variance);
    }

  //! generate K simplices
  LaGenMatDouble current_simplex(n,n);
  vector<LaGenMatDouble> simplices;
  for (int j = 0; j < k; j++)
    {
      random_generator::randomize01(current_simplex);
      simplices.push_back(current_simplex);
    }

  //! generate K hyperplanes from the respective simplex

  // the vector e ia used in the linear system Ax = b; in doing this we only lose the possibility of expressing the equation of a hyperplane passing throught the origin, i.e. we are fixing gamma = 1.
  LaVectorDouble e = LaGenMatDouble::ones(n,1);
  LaVectorDouble current_hyperplane_parameters(n,1);
  vector<LaVectorDouble> hyperplane_parameters;
  vector<double> hyperplane_gammas;

  for (int j = 0; j < k; j++)
    {
      LaLinearSolve(simplices[j], current_hyperplane_parameters, e);
      double gamma = 1;
      utility::normalize_2norm(current_hyperplane_parameters, gamma);
      hyperplane_parameters.push_back(current_hyperplane_parameters);
      hyperplane_gammas.push_back(gamma);
    }

  //! generating the noisy points
	
  //! at start the rand_j cluster is not fixed yet
  bool old_point_has_been_shot = true;
  for (int i = 0, shooting_trials = 0; i < m; i++, shooting_trials++)
    {
      if (shooting_trials > max_good_point_shooting_trials)
	{
	  cout << "ERROR: it seems that the method is unable to generate a good point for the current, randomly chosen, cluster; a FORCE exit has been imposed" << endl;
	  abort();
	}
		
      // changing the cluster we are shooting onto if the last shooting trial finished with a success; otherwise we don't change rand_j, keeping on shooting points on the same hyperplane
      if (old_point_has_been_shot)
	{
	  //rand_j = random_generator::integer(k);
	  rand_j = random_generator::sample_integer(population_probabilities);
	  assert(rand_j != -1);
	  old_point_has_been_shot = false;
	  shooting_trials = 0;
	}

      LaVectorDouble current_point(n);
      random_generator::randomize01(current_point);

      //! setting the last coordinate of current_point so that it lies on the rand_j-th hyperplane
      current_point(n-1) = ( hyperplane_gammas[rand_j] - current_point(LaIndex(0,n-2))  * hyperplane_parameters[rand_j] (LaIndex(0,n-2)) ) / hyperplane_parameters[rand_j] (n-1);

      //! assure that the last non-noisy coordinate is INSIDE the box
      if ( current_point(n-1) < 0 || current_point(n-1) > 1 )		//! out of the box!
	{
	  i--;
	  continue;
	}

      if (type.data_type == SEMI_RANDOM_CONTINUOS_REGRESSION)
	{
	  //! assuring that the randomly choosen cluster is maximum for that point
	  bool flag_wrong_point = false;
	  for (int j = 0; j < k; j++)
	    {
	      if (j == rand_j)
		continue;
	
	      //! assure that the last coord of the newly created point IS maximum on rand_j cluster
	      //! otherwise, the point is to be trown away and rebuilt...
	      if (current_point(n-1) < ( 1 - current_point(LaIndex(0,n-2)) * hyperplane_parameters[j](LaIndex(0,n-2)) ) / hyperplane_parameters[j](n-1) )
		{
		  flag_wrong_point = true;
		  break;
		}
	    }

	  //! enforcing misclassification, by sampling and choosing whether the current point has to be misclassified or not...
	  if (type.fraction_of_misclassified == 0)
	    {
	      if (flag_wrong_point)			// we don't want misclassified points
		{
		  i--;
		  continue;
		}
	    }
	  else
	    {
	      bool flag_misclassify_current_point = (random_generator::random() < type.fraction_of_misclassified) ? true : false;

	      if (flag_misclassify_current_point == false)	//! we were shooting a good (non misclassified) point
		{
		  if (flag_wrong_point)	//! bad point, and we wanted a GOOD one
		    {
		      i--;			//! another points has to be created
		      continue;
		    }
		}
	      else
		{
		  if (!flag_wrong_point)	//! good point, and we wanted a BAD one
		    {
		      i--;			//! another points has to be created
		      continue;
		    }
		}
	    }
	}

      //------------------
      //! adding noise
      //------------------
      double error;
      if (type.error == GUASSIAN){
	// guassian distribution error
	error = random_generator::nrandom() * sqrt (variances[rand_j]);
      } else if (type.error == UNIFORM){
	// KANIKA uniform distrubution error
	error = (2*random_generator::random()-1) * 2*sqrt(variances[rand_j]);
      } else abort();
		
      //! assuring that the random_distance is no larger than sigma_level_threshold*sqrt(variance)
      if (type.sigma_level_threshold != 0 && fabs(error) > type.sigma_level_threshold*sqrt(variances[rand_j]) )
	{
	  i--;
	  continue;
	}

      // ADDING THE NOISE      

      if (type.data_type == SEMI_RANDOM)   // normal error on every component
	{
	  for (int l = 0; l < n; l++)
	    {
	      current_point(l) = current_point(l) + hyperplane_parameters[rand_j](l) / 
		Blas_Norm2 (hyperplane_parameters[rand_j]) * error;
	    }
	  //! asserting exactness of the distance in the newly generated point
	  assert ( fabs ( Formulas::l2_norm_distance(current_point, hyperplane_parameters[rand_j], hyperplane_gammas[rand_j]) - fabs (error) ) <= 1e-6 );
	}
      else if (type.data_type == SEMI_RANDOM_CONTINUOS_REGRESSION) // vertical error on the last component only
	{	
	  current_point(n-1) = current_point(n-1) + error;
	}

      //! storing the newly created point
      pi.a[i] = current_point;
      if (pi.J(i) != -1)
	pi.assignment[i][pi.J(i)] = false;
      pi.assignment[i][rand_j] = true;
      old_point_has_been_shot = true;

    }

  if (type.fraction_of_random_points != 0)
    randomize(pi, type.fraction_of_random_points);

  //! storing the parameters of the generators
  pi.w = hyperplane_parameters;
  pi.gamma = hyperplane_gammas;

  //! storing the variances
  igpd.variances = variances;
  igpd.variance = type.variance;

  if (type.minimum_population_guarantee){
    // KANIKA
    //asserting that each cluster is at least populated to 'some' degree
    bool not_good = false;
    for (int j = 0; j < k; j++)
      {
	if (pi.M(j) < double(m)/k*0.3)
	  {
	    not_good = true;
	    break;
	  }
      }
    if (not_good == true)
      {
	cout << "The instance has some clusters with too few points. Regenerating it" << endl;
	Instance_generator::populate(pi, type, igpd);
      }
  }
}

void
Instance_generator::populate_pclustering(Problem_instance & pi, Instance_type & type, Instance_generation_printing_data & igpd)
{
  const int m = pi.M();
  const int n = pi.N();
  const int k = pi.K();

  LaVectorDouble population_probabilities(k);
  random_generator::randomize01(population_probabilities);
  utility::normalize_1norm(population_probabilities);

  int rand_j = -1;
	
  //! madness preventer
  int max_good_point_shooting_trials = 1000000;

  //! guassian error variances (different per hyperplane)
  vector<double> variances;
  for (int kk = 0; kk < k; kk++)
    {
      double discount = random_generator::random();
      discount *= 0.3;
      discount += 0.7;
      variances.push_back(discount*type.variance);
    }

  //! generate K centroids
  LaVectorDouble current_centroid(n);
  vector<LaVectorDouble> centroids;
  for (int j = 0; j < k; j++)
    {
      random_generator::randomize01(current_centroid);
      centroids.push_back(current_centroid);
    }

  //! generating the noisy points
	
  //! at start the rand_j cluster is not fixed yet
  bool old_point_has_been_shot = true;
  for (int i = 0, shooting_trials = 0; i < m; i++, shooting_trials++)
    {
      if (shooting_trials > max_good_point_shooting_trials)
	{
	  cout << "ERROR: it seems that the method is unable to generate a good point for the current, randomly chosen, cluster; a FORCE exit has been imposed" << endl;
	  abort();
	}
		
      // changing the cluster we are shooting onto if the last shooting trial finished with a success; otherwise we don't change rand_j, keeping on shooting points on the same hyperplane
      if (old_point_has_been_shot)
	{
	  //rand_j = random_generator::integer(k);
	  rand_j = random_generator::sample_integer(population_probabilities);
	  assert(rand_j != -1);
	  old_point_has_been_shot = false;
	  shooting_trials = 0;
	}

      LaVectorDouble current_point = centroids[rand_j];

      //------------------
      //! adding noise
      //------------------

      double error = random_generator::nrandom() * sqrt (variances[rand_j]);
      LaVectorDouble error_direction(n);
      random_generator::randomize(error_direction);
      utility::normalize_2norm(error_direction);
		
      //! assuring that the random_distance is no larger than sigma_level_threshold*sqrt(variance)
      if (type.sigma_level_threshold != 0 && fabs(error) > type.sigma_level_threshold*sqrt(variances[rand_j]) )
	{
	  i--;
	  continue;
	}

      // ADDING THE NOISE      

      for (int l = 0; l < n; l++)
	{
	  current_point(l) = current_point(l) + error_direction(l)* error;
	}
      //! asserting exactness of the distance in the newly generated point
      assert ( fabs ( Formulas::l2_norm_prototypal_distance(current_point, centroids[rand_j]) - fabs(error) ) <= 1e-6 );
      
      //! storing the newly created point
      pi.a[i] = current_point;
      if (pi.J(i) != -1)
	pi.assignment[i][pi.J(i)] = false;
      pi.assignment[i][rand_j] = true;
      old_point_has_been_shot = true;
    }

  if (type.fraction_of_random_points != 0)
    randomize(pi, type.fraction_of_random_points);

  //! storing the parameters of the generators
  pi.w = centroids;

  //! storing the variances
  igpd.variances = variances;
  igpd.variance = type.variance;
}


void 
Instance_generator::wave(Problem_instance & pi)
{
  //! initializations
  int k = 7;
  assert (k == pi.K());
  int m = pi.M();
  int n = 2;
  assert (n == pi.N());					
  int n_line = m / k;			//! number of points in a single line
  double h_space = 0.25;			//! step on x axes
  double v_space = 0.9;			//! step on y axes
  double h_step = h_space / n_line;
  double v_step = v_space / n_line;

  int added_points = 0;

  //! h_low_left
  for (double i = 0; i < h_space; i += h_step, added_points++)
    {
      LaVectorDouble v(2);
      v(0) = i;
      v(1) = 0.05;
      pi.a[added_points] = v;
    }

  //! v_left
  for (double i = 0.05 + v_step; i < v_space; i += v_step, added_points++)
    {
      LaVectorDouble v(2);
      v(0) = 0.25 + 0.1 * i;
      v(1) = i;
      pi.a[added_points] = v;
    }

  //! h_high_left
  for (double i = 0.25 + 0.1 * (0.1+v_space+v_step); i < h_space + 0.25 + 0.1 * (v_space+v_step); i += h_step, added_points++)
    {
      LaVectorDouble v(2);
      v(0) = i;
      v(1) = 0.95;
      pi.a[added_points] = v;
    }

  //! v_center
  for (double i = v_space + 0.05 - v_step; i > 0.05; i -= v_step, added_points++)
    {
      LaVectorDouble v(2);
      v(0) = 0.5 + 0.1 * (v_space+v_step) + 0.1 * (1-i);
      v(1) = i;
      pi.a[added_points] = v;
    }

  //! h_low_right
  for (double i = 0.5 + 0.1 * (v_space+v_step) + 0.1; i < h_space + 0.5 + 0.1 * (v_space+v_step) + 0.1; i += h_step, added_points++)
    {
      LaVectorDouble v(2);
      v(0) = i;
      v(1) = 0.05;
      pi.a[added_points] = v;
    }

  //! v_right
  for (double i = 0.75; i < v_space; i += v_step, added_points++)
    {
      LaVectorDouble v(2);
      v(0) = 0.75 + 0.05 * i;
      v(1) = i;
      pi.a[added_points] = v;
    }

  //! h_high_right
  for (double i = 0.75; i < h_space + 0.75; i += h_step, added_points++)
    {
      LaVectorDouble v(2);
      v(0) = i;
      v(1) = 1;
      pi.a[added_points] = v;
    }

  if (added_points != m)
    {
      cout << "Error! Francesco's awful way to create a wave has been unable to generate the right amount of points!" << endl;
      cout << "Actually, it generated " << added_points << " out of " << m << "!. Halting!" << endl;
      abort();
    }
}


void 
Instance_generator::populate_srnc(Problem_instance & pi, Instance_type& type, string true_param)
{
  
  
  
   int k=pi.K();
   
   int n=pi.N();
	int m = pi.M();
  
  cout<<"K: "<<k<<endl;
  cout<<"N: "<<n<<endl;
  cout<<"M: "<<m<<endl;
  
   vector<double> variances;
  for (int kk = 0; kk < k; kk++)
    {
      double discount = random_generator::random();
      discount *= 0.3;
      discount += 0.7;
      variances.push_back(discount*type.variance);
    }
  
  
   LaVectorDouble e = LaGenMatDouble::ones(n-1,1);
  LaVectorDouble current_centroid(n-1);
  vector<LaVectorDouble> centroids;
  for (int j = 0; j < k; j++)
    {
      
      random_generator::randomize01(current_centroid);
     
     bool ok=false;
     
     if(j>=1)
		while(fabs(current_centroid(n-2) - centroids[j-1](n-2))< (0.5/k) || !ok)//at least a bit 'far' on a dimension (not all centroids collapsed in one point)
				{
					random_generator::randomize01(current_centroid);
					
					ok=true; //avoid having identical centroids
					for (int jj = 0; jj < j; jj++)
							{
								bool uguali=true;
								for(int i=0;i<n-1;i++)
									if(fabs(current_centroid(i)-centroids[jj](i))>=0.0005)
										{uguali=false;}
									
								if(uguali==true)
									{ok=false; break;}
							}
					}
     
      centroids.push_back(current_centroid);
      
    }
  
  vector <double> Gs;
  
  vector<LaVectorDouble> Ws;
  
  max_margin_RLP(centroids, Ws, Gs);

//per appoggiarli da qualche parte.
for (int j = 0; j < k; j++)
    {
		for(int i=0;i<n-1;i++)
			pi.w[j](i) = centroids[j](i);
			
	pi.w[j](n-1) = 0;
	}


cout<<"Centroids:"<<endl;
cout<<centroids<<endl;

 /* LaGenMatDouble current_simplex(n-1,n-1);
  vector<LaGenMatDouble> simplices;
  for (int j = 0; j < k; j++)
    {
      random_generator::randomize01(current_simplex);
      simplices.push_back(current_simplex);
    }

  //! generate K hyperplanes from the respective simplex

  // the vector e ia used in the linear system Ax = b; in doing this we only lose the possibility of expressing the equation of a hyperplane passing throught the origin, i.e. we are fixing gamma = 1.
  LaVectorDouble e = LaGenMatDouble::ones(n-1,1);
  LaVectorDouble current_hyperplane_parameters(n-1,1);
  vector<LaVectorDouble> hyperplane_parameters;
  vector<double> hyperplane_gammas;

  for (int j = 0; j < k; j++)
    {
      LaLinearSolve(simplices[j], current_hyperplane_parameters, e);
      double gamma = 1;
      utility::normalize_2norm(current_hyperplane_parameters, gamma);
      hyperplane_parameters.push_back(current_hyperplane_parameters);
      hyperplane_gammas.push_back(gamma);
    }
*/  
  
  string true_param_file = true_param + "_true";
  cerr<<"bl"<<true_param_file;
  ofstream ampl_data(true_param_file.c_str());
  cerr<<"bl";
  LaVectorDouble current_hyperplane_parameters(n-1,1);
  vector<LaVectorDouble> hyperplane_parameters;
  vector<double> hyperplane_gammas;
   for (int j = 0; j < k; j++)
    {
		random_generator::randomize01(current_hyperplane_parameters);
		current_hyperplane_parameters = current_hyperplane_parameters - (0.5)*e;
		
		
		
		cout<<current_hyperplane_parameters<<endl;
		double gamma=0.5;
		//utility::normalize_2norm(current_hyperplane_parameters, gamma);
		hyperplane_parameters.push_back(current_hyperplane_parameters);
		cout<<current_hyperplane_parameters<<endl;
		hyperplane_gammas.push_back(gamma);
	
	}
  
   cerr<<"bl";
  ampl_data << "# True hyperplanes parameters \n";
  ampl_data << "w := \n";
  for(int j=0;j<k;j++)
	for (int i = 0; i < n-1; i++)
		ampl_data << j+1 <<" "<<i+1<<"\t"<<hyperplane_parameters[j](i)<<"\n";
  ampl_data << ";" << endl;
  
  
  ampl_data << "gamma [*] := \n" << endl;
  for(int j=0;j<k;j++)
    {
      ampl_data << j+1 << " "<<hyperplane_gammas[j]<<"\n";
      
    }
  ampl_data << ";" << endl;
  
   
  
  
  
  
  
  
  
  int misclass=0;
    bool old_point_has_been_shot = true;
  for (int i = 0, shooting_trials = 0; i < m; i++, shooting_trials++)
    {
      if (shooting_trials > 100000)
	{
	  cout << "ERROR: it seems that the method is unable to generate a good point for the current, randomly chosen, cluster; a FORCE exit has been imposed" << endl;
	  abort();
	}
		
      LaVectorDouble current_point(n-1);
      random_generator::randomize01(current_point);
      
		double	max = -1e300;
		double		curr=-1; 
		int			max_cluster = -1;
		   
      
      for (int j = 0; j < k; j++)
	  {
		  if(0)
		  {cout<<"vettore W: "<<Ws[j]<<endl;
		  cout<<"punto corrente: "<<current_point<<endl;
		  cout<<"ultimo componente punto: "<<current_point(n)<<endl;
		 // cout<<"w*a+pan "<<pi.W[j] * pi.a[i] + pi.a[i](n)<<endl;
			}
		
		// check N/n/n-1...
			assert(Ws[j].size()==n-1);
			curr = Ws[j] * current_point - Gs[j]; // the last component of W is set to -1 (the last component of 'a' is the dependent variable y, which is not involved in the partitioning). In the scalar product, it is summed up multiplied by (-1); adding the last component again we remove it.
																		
		
															
	    if (curr > max)
	      {
		max = curr;
		max_cluster = j;
	      }
	  }
      
      
      
      LaVectorDouble point_and_y(n);
      
      point_and_y = -1;
      
      for(int ii=0;ii<n-1;ii++)
		{point_and_y(ii)=current_point(ii); }
      
      if (type.fraction_of_misclassified == 0 || (misclass*100/m)>=type.fraction_of_misclassified)
			point_and_y(n-1) = hyperplane_parameters[max_cluster]*current_point + hyperplane_gammas[max_cluster];
		else
			{  int rand_cluster= rand()%k;
				
				point_and_y(n-1) = hyperplane_parameters[rand_cluster]*current_point + hyperplane_gammas[rand_cluster];
				if(rand_cluster!=max_cluster)
					{misclass++;
					cout<<"miscla "<<misclass<<endl;
					}
				
			}
      

      //------------------
      //! adding noise
      //------------------
      double error;
      
      
     	// guassian distribution error
	error = random_generator::nrandom() * sqrt (variances[max_cluster]);
    		
   

      // ADDING THE NOISE      

     point_and_y(n-1) = point_and_y(n-1)  + error;
	

      //! storing the newly created point
      pi.a[i] = point_and_y;
      //cout<<i<<" "<<pi.a[i]<<endl;
      if (pi.J(i) != -1)
	pi.assignment[i][pi.J(i)] = false;
      pi.assignment[i][max_cluster] = true;
      

    }

  
  
  
  
  
}  
  
  void 
Instance_generator::populate_src(Problem_instance & pi, Instance_type& type, string true_param)
{
  
  
  
   int k=pi.K();
   
   int n=pi.N();
	int m = pi.M();
  
  cout<<"K: "<<k<<endl;
  cout<<"N: "<<n<<endl;
  cout<<"M: "<<m<<endl;
  
   vector<double> variances;
  for (int kk = 0; kk < k; kk++)
    {
      double discount = random_generator::random();
      discount *= 0.3;
      discount += 0.7;
      variances.push_back(discount*type.variance);
    }
  
  
   LaVectorDouble e = LaGenMatDouble::ones(n-1,1);
  
  

 /* LaGenMatDouble current_simplex(n-1,n-1);
  vector<LaGenMatDouble> simplices;
  for (int j = 0; j < k; j++)
    {
      random_generator::randomize01(current_simplex);
      simplices.push_back(current_simplex);
    }

  //! generate K hyperplanes from the respective simplex

  // the vector e ia used in the linear system Ax = b; in doing this we only lose the possibility of expressing the equation of a hyperplane passing throught the origin, i.e. we are fixing gamma = 1.
  LaVectorDouble e = LaGenMatDouble::ones(n-1,1);
  LaVectorDouble current_hyperplane_parameters(n-1,1);
  vector<LaVectorDouble> hyperplane_parameters;
  vector<double> hyperplane_gammas;

  for (int j = 0; j < k; j++)
    {
      LaLinearSolve(simplices[j], current_hyperplane_parameters, e);
      double gamma = 1;
      utility::normalize_2norm(current_hyperplane_parameters, gamma);
      hyperplane_parameters.push_back(current_hyperplane_parameters);
      hyperplane_gammas.push_back(gamma);
    }
*/  
  
  string true_param_file = true_param + "_true";
  cerr<<"bl"<<true_param_file;
  ofstream ampl_data(true_param_file.c_str());
  cerr<<"bl";
  LaVectorDouble current_hyperplane_parameters(n-1,1);
  vector<LaVectorDouble> hyperplane_parameters;
  vector<double> hyperplane_gammas;
  vector<int> points;
   for (int j = 0; j < k; j++)
    {
		points.push_back(0);
		random_generator::randomize01(current_hyperplane_parameters);
		current_hyperplane_parameters = 2*current_hyperplane_parameters - (0.5)*e;
		
		
		
		cout<<current_hyperplane_parameters<<endl;
		double gamma=random_generator::random();
		//utility::normalize_2norm(current_hyperplane_parameters, gamma);
		hyperplane_parameters.push_back(current_hyperplane_parameters);
		cout<<current_hyperplane_parameters<<endl;
		hyperplane_gammas.push_back(gamma);
	
	}
  
   cerr<<"bl";
  ampl_data << "# True hyperplanes parameters \n";
  ampl_data << "w := \n";
  for(int j=0;j<k;j++)
	for (int i = 0; i < n-1; i++)
		ampl_data << j+1 <<" "<<i+1<<"\t"<<hyperplane_parameters[j](i)<<"\n";
  ampl_data << ";" << endl;
  
  
  ampl_data << "gamma [*] := \n" << endl;
  for(int j=0;j<k;j++)
    {
      ampl_data << j+1 << " "<<hyperplane_gammas[j]<<"\n";
      
    }
  ampl_data << ";" << endl;
  
   
  
  
   
  
  
  
  int misclass=0;
    bool old_point_has_been_shot = true;
  for (int i = 0, shooting_trials = 0; i < m; i++, shooting_trials++)
    {
      if (shooting_trials > 100000)
	{
	  cout << "ERROR: it seems that the method is unable to generate a good point for the current, randomly chosen, cluster; a FORCE exit has been imposed" << endl;
	  abort();
	}
		
      LaVectorDouble current_point(n-1);
      random_generator::randomize01(current_point);
      
		double	max = -1e300;
		double		curr=-1e300; 
		int			max_cluster = -1;
		   
      
      for (int j = 0; j < k; j++)
	  {
		  if(0)
		  {		  cout<<"punto corrente: "<<current_point<<endl;
		  cout<<"ultimo componente punto: "<<current_point(n)<<endl;
		 
			}
		
																		
		curr=hyperplane_parameters[j]*current_point + hyperplane_gammas[j];
			
			
															
	    if (curr > max)
	      {
		max = curr;
		max_cluster = j;
	      }
	  }
      
      points[max_cluster]++;
      
      LaVectorDouble point_and_y(n);
      
      point_and_y = -1;
      
      for(int ii=0;ii<n-1;ii++)
		{point_and_y(ii)=current_point(ii); }
      
      if (type.fraction_of_misclassified == 0 || (misclass*100/m)>=type.fraction_of_misclassified)
			point_and_y(n-1) = hyperplane_parameters[max_cluster]*current_point + hyperplane_gammas[max_cluster];
		else
			{  int rand_cluster= rand()%k;
				
				point_and_y(n-1) = hyperplane_parameters[rand_cluster]*current_point + hyperplane_gammas[rand_cluster];
				if(rand_cluster!=max_cluster)
					{misclass++;
					cout<<"miscla "<<misclass<<endl;
					}
				
			}
      

      //------------------
      //! adding noise
      //------------------
      double error;
      
      
     	// guassian distribution error
	error = random_generator::nrandom() * sqrt (variances[max_cluster]);
    		
   

      // ADDING THE NOISE      

     point_and_y(n-1) = point_and_y(n-1)  + error;
	

      //! storing the newly created point
      pi.a[i] = point_and_y;
      //cout<<i<<" "<<pi.a[i]<<endl;
      if (pi.J(i) != -1)
	pi.assignment[i][pi.J(i)] = false;
      pi.assignment[i][max_cluster] = true;
      

    }

  
  cout<<points;
  int nonempty=0;
      for (int j = 0; j < k; j++)
			if(points[j]>1)
					nonempty++;
  
  cout<<"Real k: "<<nonempty<<endl;
  
  
}  
  
 void 
Instance_generator::rwave(Problem_instance & pi)
{
  
  
  //! initializations
  int k = 0;
  int m = pi.M();
  int n = 2;
  assert (n == pi.N());					
  int n_line = m / pi.K();			//! number of points in a single line
  double h_space = 0.25;			//! step on x axes
  double v_space = 0.9;			//! step on y axes
  //double h_step = h_space / n_line;
  //double v_step = v_space / n_line;
  double x=0;
  double y=0;

  int added_points = 0;
double prob=0.8;
cout<<"seed: "<<random_generator::get_seed()<<endl;
	random_generator::set_seed(random_generator::get_seed()); //leo: the seed is USELESS?
while(added_points<m)
{
	LaVectorDouble v(2);
      v(0) = x;
      v(1) = y;
      pi.a[added_points] = v;
      added_points++;
	
	
	if(added_points%n_line==0 && k<pi.K()-1)
	//if(random_generator::random()>prob)
	{
		//h_space=random_generator::random()/2;
		v_space=random_generator::random()*pow((double)-1,k+1);
		k++;	
		if(k==pi.K())
			prob=1;
		}
	
	x+=h_space;
	y+=v_space;
		
	
	}

  
  if (pi.K()!=k+1)
    {
      cout << "Error! Leonardo's awful way to create a wave has been unable to generate the right amount of clusters!" << endl;
      cout << "Actually, it generated " << k+1 << " out of " << pi.K() << "!." << endl;
      
    }

		


}


void
Instance_generator::write_points_to_file(Problem_instance & pi, string & file_path, Instance_generation_printing_data & igpd)
{
  ofstream instance_file(file_path.c_str());
  instance_file << igpd.name << endl;
  instance_file << "m: " << pi.M() << endl;
  instance_file << "n: " << pi.N() << endl;
  instance_file << "k: " << pi.K() << endl;
  instance_file << "seed: " << igpd.seed << endl;
  instance_file << "Variance: " << igpd.variance << endl;
  for (int k = 0; k < pi.K(); k++)
    {
      instance_file << "Var[" << k << "] = " << igpd.variances[k] << endl;
    }
  instance_file << "Data_set:" << endl;
  for (int i = 0; i < pi.M(); i++)
    {
      instance_file << "p:";
      for (int l = 0; l < pi.N(); l++)
	instance_file << "  " << pi.a[i](l);
      instance_file << endl;
    }
}

void
Instance_generator::randomize(Problem_instance & pi, double ratio)
{	
  int m = pi.M();
  int n = pi.N();
  vector<double>	min_coords(n, -1e300);
  vector<double>	max_coords(n, 1e300);
	
  //! find min and max coordinates in the point set
  pi.get_min_max_dimensions(min_coords, max_coords);

  for (int i = 0; i < m; i++)
    if (ratio < random_generator::random())
      for (int l = 0; l < n; l++)
	pi.a[i](l) = random_generator::random_in_box(min_coords[l], max_coords[l]);
}

void
Instance_generator::max_margin_RLP(const vector<LaVectorDouble> & setof_domain_representative_points, vector<LaVectorDouble> & Ws, vector<double>  & Gs)
{
  int k = setof_domain_representative_points.size();
  cerr<<"k: "<<k;
  int n = setof_domain_representative_points[0].size();	// D is already D-1, i.e. it is the number of dimensions of the DOMAIN

  //! deleting the vectors
  Ws.clear();
  Gs.clear();		
  
  LaVectorDouble dummy_LaVectorDouble(n);
  dummy_LaVectorDouble = -1;
  Ws = vector<LaVectorDouble>(k, dummy_LaVectorDouble);
  Gs = vector<double>(k, -1);

  //--------------------------
  //! generate model
  //--------------------------

  ofstream ampl_model("instance_generation_RLP.mod");

  //! variables and parameters
  ampl_model << "\
	set N;	#space dimensions						\n\
	set K;	#cluster (and point) space					\n\
	var G {K};								\n\
	var W {K, N};								\n\
	var inverse_margin {K, K} >= 0;						\n\
	param a {K, N};								\n\
	                                                                        \n\
	minimize obj:								\n\
	sum{i in K, j in K: i <> j} (inverse_margin[i,j]);			\n\
	                                                                        \n\
	subject to								\n\
	missclassification_error {i in K, j in K: i <> j}:			\n\
	sum{l in N} (a[i,l] * ( W[i,l] - W[j,l])  ) - (G[i] - G[j]) >= 1;	\n\
	                                                                        \n\
	inverse_margin_constraint_1 {i in K, j in K, l in N: i < j}:	        \n\
	inverse_margin[i,j] >= W[i,l] - W[j,l];					\n\
	                                                                        \n\
	inverse_margin_constraint_2 {i in K, j in K, l in N: i < j}:	        \n\
	inverse_margin[i,j] >= - W[i,l] + W[j,l];				\n\
	" << endl;

  //--------------------------
  //! generate data
  //--------------------------

  ofstream ampl_data("instance_generation_RLP.dat");
  
  ampl_data << "set N := "; for (int i = 1; i <= n; i++) ampl_data << i << " "; ampl_data << ";" << endl;
  ampl_data << "set K := "; for (int i = 1; i <= k; i++) ampl_data << i << " "; ampl_data << ";" << endl;

   ampl_data << "param a : \n";
  for (int i = 1; i <= n; i++)
    ampl_data << " " << i;
  ampl_data << ":=" << endl;
  for (int i = 0; i < k; i++)
    {
      ampl_data << i+1 << " ";
      for (int l = 0; l < n; l++)
	ampl_data << " " << setof_domain_representative_points[i](l);
      ampl_data << "\n";
    } 
    
  ampl_data << ";" << endl;
 
  //--------------------
  //! generate run
  //--------------------
  ofstream ampl_run("instance_generation_RLP.run");
	
  ampl_run << "\
	option solver cplexamp;													\n\
	model instance_generation_RLP.mod;									\n\
	data instance_generation_RLP.dat;									\n\
	solve;													\n\
	display obj > instance_generation_RLP.output;							\n\
	display W > instance_generation_RLP.output;							\n\
	display G > instance_generation_RLP.output;							\n\
	display inverse_margin > instance_generation_RLP.output;					\n\
	" << endl;

  //--------------------
  //! execute the linear program
  //--------------------

  system("ampl instance_generation_RLP.run");

  //--------------------
  //! fetch the results
  //--------------------


  ampl_output_reader ampl_output("instance_generation_RLP.output");
  string name_obj("obj");
  double obj;
  ampl_output.read_value(name_obj, obj);

  string name_W("W");
  ampl_output.read_matrix(name_W, Ws, k, n);

  string name_G("G");
  ampl_output.read_vector(name_G, Gs, k);
}


