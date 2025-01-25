#include "random_walk.h"

#include <fstream>
#include <iterator>

RandomWalk::RandomWalk( std::function<double()> next_incr_generator) : m_gen_next_incr(next_incr_generator)
{
  m_walk.push_back(0.0);
  num_steps = 0;
}

void RandomWalk::gen_steps(unsigned int n)
{
  for (int i=0; i<n; ++i) {
    m_walk.push_back(m_walk.back() + m_gen_next_incr());
    num_steps++;
  }
}

const deque<double>& RandomWalk::get_walk() {
  return m_walk;
}

void RandomWalk::save_walk(std::string filename)
{
  std::ofstream output(filename);
  std::ostream_iterator<double> output_it(output, "\t");
  std::copy(m_walk.begin(), m_walk.end(), output_it);
}
