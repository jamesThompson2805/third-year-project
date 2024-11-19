#include "random_walk.h"

RandomWalk::RandomWalk( std::function<double()> next_incr_generator) : m_gen_next_incr(next_incr_generator)
{
  m_walk.push_back(0.0);
}

void RandomWalk::gen_steps(unsigned int n)
{
  for (int i=0; i<n; ++i) {
    m_walk.push_back(m_walk.back() + m_gen_next_incr());
  }
}

const deque<double>& RandomWalk::get_walk() {
  return m_walk;
}
