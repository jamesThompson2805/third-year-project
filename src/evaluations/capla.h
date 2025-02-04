#ifndef EVAL_CAPLA_H
#define EVAL_CAPLA_H

#include <vector>
#include "general.h"

/* What to assess for CAPLA
 * - How fast is it?
 * - How accurate is it?
 *
 * Speed:
 *  - is it faster than d_w
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
 *  - is it more accurate than d_w
 *  - is it more accurate than d_p min mse
 *  - is it more accurate than d_p min maxdev
 *  - is it more accurate than pla
 *  How does it change when you vary:
 *  - number of parameters
 *  - size of windows
 *  - shape of windows
 *  - size of data
 *  How to assess accuracy
 *  - MEAN Squared Error
 *  - Maximum Deviation 
 *
 */

namespace capla_eval {
double mse_of_method(const Seqd& s, unsigned int num_params, const Seqd& ldist, const Seqd& rdist);
double maxdev_of_method(const Seqd& s, unsigned int num_params, const Seqd& ldist, const Seqd& rdist);
double cputime_ms_of_method(const Seqd& s, unsigned int num_params, const Seqd& ldist, const Seqd& rdist);

Seqd get_mse_over_sizes(const Seqd& s, Sequi sizes, unsigned int num_params, const Seqd& ldist, const Seqd& rdist);
Seqd get_maxdev_over_sizes(const Seqd& s, Sequi sizes, unsigned int num_params, const Seqd& ldist, const Seqd& rdist);
Seqd get_cputime_over_sizes(const Seqd& s, Sequi sizes, unsigned int num_params, const Seqd& ldist, const Seqd& rdist);

Seqd get_mse_over_num_params(const Seqd& s, Sequi vec_params, const Seqd& ldist, const Seqd& rdist);
Seqd get_maxdev_over_num_params(const Seqd& s, Sequi vec_params, const Seqd& ldist, const Seqd& rdist);
Seqd get_cputime_over_num_params(const Seqd& s, Sequi vec_params, const Seqd& ldist, const Seqd& rdist);

Seqd get_mse_over_window_window_sizes_method_mean(const Seqd& s, Sequi window_sizes, unsigned int num_params);
Seqd get_maxdev_over_window_window_sizes_method_mean(const Seqd& s, Sequi window_sizes, unsigned int num_params);
Seqd get_cputime_over_window_window_sizes_method_mean(const Seqd& s, Sequi window_sizes, unsigned int num_params);

Seqd get_mse_over_window_window_sizes_method_miss_start(const Seqd& s, Sequi window_sizes, unsigned int num_params);
Seqd get_maxdev_over_window_window_sizes_method_miss_start(const Seqd& s, Sequi window_sizes, unsigned int num_params);
Seqd get_cputime_over_window_window_sizes_method_miss_start(const Seqd& s, Sequi window_sizes, unsigned int num_params);

Seqd get_mse_over_window_window_sizes_method_tri(const Seqd& s, Sequi window_sizes, unsigned int num_params);
Seqd get_maxdev_over_window_window_sizes_method_tri(const Seqd& s, Sequi window_sizes, unsigned int num_params);
Seqd get_cputime_over_window_window_sizes_method_tri(const Seqd& s, Sequi window_sizes, unsigned int num_params);

Seqd get_mse_over_window_window_sizes_method_tri_miss_start(const Seqd& s, Sequi window_sizes, unsigned int num_params);
Seqd get_maxdev_over_window_window_sizes_method_tri_miss_start(const Seqd& s, Sequi window_sizes, unsigned int num_params);
Seqd get_cputime_over_window_window_sizes_method_tri_miss_start(const Seqd& s, Sequi window_sizes, unsigned int num_params);
};


#endif
