#include "inclusions.h"

#ifndef __INSTANCE_GENERATOR__
#define __INSTANCE_GENERATOR__

class Instance_generator {
 public:
  // PUBLIC INTERFACE: random dataset generators
  //! generates all kind of instances
  static void populate(Problem_instance & pi, Instance_type & type, Instance_generation_printing_data & igpd, base_generator_type & generator );
  
  static void populate_pclustering(Problem_instance & pi, Instance_type & type, Instance_generation_printing_data & igpd, base_generator_type & generator );

  static void populate_srnc(Problem_instance & pi, Instance_type& type, base_generator_type & generator );
  
  
  // PUBLIC INSTANCE: data set text file generators
  
  //! outputs to file
  static void write_points_to_file(Problem_instance & pi, string & file_path, Instance_generation_printing_data & igpd);
  
  //! solves a max_margin RLP problem on a set of linearly separable points; used when generating instances by the instance generator
  static void max_margin_RLP(const vector<LaVectorDouble> & setof_domain_representative_points, vector<LaVectorDouble> & Ws, vector<double> & Gs);
  
 private:
  //! generates a 2-d wave with a fixed n = 2 and k = 5;
  static void wave(Problem_instance & pi);
  
  //! randomizes ratio*m of the points in the instance, substituting to their coordinates new random ones
  static void randomize(Problem_instance & pi, double ratio, base_generator_type & generator );
  

};

#endif
