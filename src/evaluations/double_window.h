#ifndef EVAL_DOUBLE_WINDOW_H
#define EVAL_DOUBLE_WINDOW_H

#include <vector>

/* What to assess for double window
 * - How fast is it?
 * - How accurate is it?
 *
 * Speed:
 *  - is it faster than d_p
 *  - is it faster than pla
 *  How does it change when you vary:
 *  - number of parameters
 *  - size of windows
 *  - size of data
 *  How to assess speed:
 *  - Wall Clock Times
 *  - CPU Times
 *  - Page Accesses?
 *
 *  Accuracy:
 *  - is it more accurate than d_p min mse
 *  - is it more accurate than d_p min maxdev
 *  - is it more accurate than pla
 *  How does it change when you vary:
 *  - number of parameters
 *  - size of windows
 *  - size of data
 *  How to assess accuracy
 *  - MEAN Squared Error
 *  - Maximum Deviation 
 *
 */

typedef std::vector<double> Seq;

namespace double_window_eval {

double mse_of_method(const Seq& s, unsigned int num_params, const Seq& ldist, const Seq& rdist);
double maxdev_of_method(const Seq& s, unsigned int num_params, const Seq& ldist, const Seq& rdist);
double cputime_ms_of_method(const Seq& s, unsigned int num_params, const Seq& ldist, const Seq& rdist);

Seq get_mse_over_sizes(const Seq& s, std::vector<unsigned int> sizes, unsigned int num_params, const Seq& ldist, const Seq& rdist);
Seq get_maxdev_over_sizes(const Seq& s, std::vector<unsigned int> sizes, unsigned int num_params, const Seq& ldist, const Seq& rdist);
Seq get_cputime_over_sizes(const Seq& s, std::vector<unsigned int> sizes, unsigned int num_params, const Seq& ldist, const Seq& rdist);

Seq get_mse_over_num_params(const Seq& s, std::vector<unsigned int> vec_params, const Seq& ldist, const Seq& rdist);
Seq get_maxdev_over_num_params(const Seq& s, std::vector<unsigned int> vec_params, const Seq& ldist, const Seq& rdist);
Seq get_cputime_over_num_params(const Seq& s, std::vector<unsigned int> vec_params, const Seq& ldist, const Seq& rdist);

Seq get_mse_over_window_window_sizes(const Seq& s, std::vector<unsigned int> window_sizes, unsigned int num_params);
Seq get_maxdev_over_window_window_sizes(const Seq& s, std::vector<unsigned int> window_sizes, unsigned int num_params);
Seq get_cputime_over_window_window_sizes(const Seq& s, std::vector<unsigned int> window_sizes, unsigned int num_params);

};

#endif
