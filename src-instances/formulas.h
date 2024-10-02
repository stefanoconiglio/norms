#include"inclusions.h"

#ifndef __FORMULAS__
#define __FORMULAS__

class Formulas
{
 public:
  // 'a la' k-means distance
  static double l2_norm_prototypal_distance(LaVectorDouble & a, LaVectorDouble & w);

  //! l2 norm orthogonal distance	--returned in abs
  static double l2_norm_distance(LaVectorDouble & a, LaVectorDouble & w, double gamma);
	
  //! l1 norm orthogonal distance	--returned in abs
  static double l1_norm_distance(LaVectorDouble & a, LaVectorDouble & w, double gamma);
	
  //! residual (it is norm independent) --returned in abs
  static double residual(LaVectorDouble & a, LaVectorDouble & w, double gamma);
};

#endif


