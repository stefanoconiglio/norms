#include "inclusions.h"

#ifndef __EXACT__
#define __EXACT__


namespace exact {

    //! outputs a 0-1 vector with 0 if the point is sampled out, 1 if kept, on the basis of the probabilities in row_probabilities and the policy specified in sampling_threshold; if == 0, samples (out) m times according to the cumulative distribution built over row_probabilities; if > 0, it samples out the points with a probability higher than the AVERAGE probability multiplied by sampling_threshold; if < 0, it samples out the m*sampling_threshold points with higher probability (sorting them)
    std::vector<bool> sample_by_policy(double sampling_threshold, vector<double> & row_probabilities, base_generator_type & generator );

    //! computes the exact gamma according to Mangasarian
    double gamma(const LaGenMatDouble & A, const LaVectorDouble & w);

    //! computes the gamma for a weigthed problem with weight p
    double gamma(const LaGenMatDouble & A, const LaVectorDouble & w, const LaVectorDouble & p);

    //! ricalualtes optimal w and gamma according to Mangasarian's algorithm; it peforms the weighted update in case a weight vector p is given AND it performs sampling based on svd or exact probabilities
    // void recalculate_l2_norm_hyperplane_parameters(Problem_instance & pi, int j, LaVectorDouble * p, bool flag_do_sampling, bool flag_worsening, base_generator_type & generator );

    void plain_recalculate_l2_norm_hyperplane_parameters(Problem_instance & pi, int j);

    //! centroids update
    //void recalculate_l2_norm_centroids(Problem_instance & pi, int j);

    //! ricalualtes optimal w and gamma according to LS; it peforms the weighted update in case a weight vector p is given AND it performs sampling based on svd or exact probabilities
    // void recalculate_l2_norm_affine_submodel_parameters(Problem_instance & pi, int j, LaVectorDouble * p, bool flag_do_sampling, base_generator_type & generator );

    //! ricalualtes optimal w and gamma according to LS; it peforms the weighted update in case a weight vector p is given AND it performs sampling based on svd or exact probabilities; it enforces gamma = 0 and w >= 0
    // void recalculate_l2_norm_linear_submodel_parameters(Problem_instance & pi, int j, LaVectorDouble * p, bool flag_do_sampling, base_generator_type & generator );

    //! updates the cluster parameters by means of the Inverse Power Method
    void update_l2_norm_hyperplane_parameters(Problem_instance & pi, int j);
 

    //! returns the exact sampling probabilities (obtained via the power method) --note: they must be rescaled s.t. they sum up to 1
    //! it can even produce a comparing print of the sampling probabilities per row with the decrease in the objective obtained by removing that row
    //! the probabilities are simply the cluster measures after removing each point and reupdating the clusters (that's why they MUST be rescaled)
    vector<double> exact_l2_norm_hyperplane_sampling_probabilities(const LaGenMatDouble & A);

    //! NOTE: possibly unused/useless
    //! exact sampling (to be rescaled) probabilities for a 2-norm min ||Ax-b|| problem
    vector<double> exact_l2_norm_regression_sampling_probabilities(const LaGenMatDouble & A, const LaVectorDouble & b);


    //! standard B = A' * (I - ee'/m) * A.
    LaGenMatDouble B(const LaGenMatDouble & A);

    //! computes weighted B matrix: B = A^T P^T * (I - p p^T / p^Tp) * A *P, accordingly to the weights vector p;
    //! P is the diagonal matrix containing the weights; it rescales the rows in A by those values
    //! in case some points are sampled OUT of A, it is sufficient to pass a vector of ones as p, setting it to 0 per each sampled out point
    LaGenMatDouble B(const LaGenMatDouble & A, const LaVectorDouble & p);

    //! computes the norm-independent hyperplane parameters for the case in which m_j == n (even m_j <. n in reality); NOT WORKING FOR HYPERPLANES CONTAINING THE ORIGIN!
    //void calculate_exact_fitting_parameters(Problem_instance & pi, int j);


}

#endif


