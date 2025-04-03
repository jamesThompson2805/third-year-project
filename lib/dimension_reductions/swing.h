#ifndef SWING_H
#define SWING_H
#include "pla.h"

namespace swing {
  Seqddt swing(const Seqd& s, double epsilon);
  std::vector<std::tuple<double, unsigned int>> swing_compr(const Seqd& s, double epsilon);

  Seqddt slide(const Seqd& s, double epsilon);

};
#endif

