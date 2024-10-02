#include "formulas.h"

double
Formulas::l2_norm_prototypal_distance(LaVectorDouble & a, LaVectorDouble & w)
{
  return Blas_Norm2(a - w);
}

double
Formulas::l2_norm_distance(LaVectorDouble & a, LaVectorDouble & w, double gamma) {

  // cout << "Formulas::l2_norm_distance" << endl;
  // cout << "w " << w << endl;
  // cout << "gamma " << gamma << endl;
  // cout << fabs 	( a * w - gamma ) << endl;
  // cout << Blas_Norm2(w) << endl;
  
  return 	( fabs 	( a * w - gamma )
		  /
		  Blas_Norm2(w) 
		  );
}

double
Formulas::l1_norm_distance(LaVectorDouble & a, LaVectorDouble & w, double gamma)
{
	return 	( fabs 	( a * w - gamma )
// 					/
// 					Blas_Norm1(*w) 
// ?!?? if it's two norm, we should divide it by inf norm....
			);
}

//! NOTE given ax + by = gamma, in fixing b = +1 or -1 we get
//! 1) b = -1 ==> ax -y = gamma ==> y = +ax - gamma ==> r = y - ax + gamma
//! 2) b = +1 ==> ax +y = gamma ==> y = -ax + gamma ==> r = y + ax - gamma
//! ------> they are CLEARLY equivalent, put that r is taken in abs
double
Formulas::residual(LaVectorDouble & a, LaVectorDouble & w, double gamma)
{
	int n = a.size();
	if (w(n-1) != -1)
		cout << "diverso! " << w(n-1) << endl;
	assert( w(n-1) == -1 );
	double to_return;
	to_return = fabs 	( a * w - gamma );
	double sure_residual = 0;
	sure_residual = -a(n-1);
	sure_residual -= gamma;
	for (int l = 0; l < n-1; l++)
		sure_residual += a(l) * w(l);
	sure_residual = fabs(sure_residual);
	if ( fabs(to_return - sure_residual) >= 1e-06 )
	{
		cout << "ERROR in double formulas::residual(LaVectorDouble * a, LaVectorDouble * w, double gamma)" << endl;
		cout << to_return << endl;
		cout << sure_residual << endl;
		abort();
	}
	return to_return;
}
