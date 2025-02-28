#ifndef BOTTOM_UP_H
#define BOTTOM_UP_H



#include "pla.h"
#include "apla_segment_and_merge.h"


namespace bottom_up {
  using ERROR_F = double (*)(const double *const, const DoublePair&, unsigned int);
  double se( const double *const, const DoublePair&, unsigned int);
  Seqddt bottom_up(const Seqd&, double, ERROR_F);
};




#endif
