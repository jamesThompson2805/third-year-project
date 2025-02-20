#ifndef EVAL_CAPLA_H
#define EVAL_CAPLA_H


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

#include "pla.h"
using DRT_COMPR = std::function<std::vector<std::tuple<DoublePair, unsigned int>> (const Seqd&,unsigned int)>;

namespace capla_eval {
  DRT generate_mean_DRT(unsigned int win_size);
  DRT_COMPR generate_mean_DRT_COMPR(unsigned int win_size);
  DRT generate_mean_skip_one_DRT(unsigned int win_size);
  DRT generate_tri_DRT(unsigned int win_size);
  DRT generate_tri_skip_one_DRT(unsigned int win_size);
};


#endif
