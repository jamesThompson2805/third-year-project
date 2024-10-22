#include <iostream>

#include "random_walk.h"

int main()
{
  UnifFunctor uni;
  RandomWalk walk(uni);

  walk.gen_steps(100);

  std::cout << walk.get_walk().front() << std::endl;

  return 0;
}
