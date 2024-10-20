#include <random>
#include <functional>
#include <deque>
#include <string>
using std::deque;


class RandomWalk {
private:
  deque<float> m_walk;
  std::function< float() > m_gen_next_incr;

public:
  RandomWalk(std::function<float()> next_incr_generator);
  void gen_steps(unsigned int n);
  const deque<float>& get_walk();
};

RandomWalk::RandomWalk( std::function<float()> next_incr_generator) : m_gen_next_incr(next_incr_generator)
{
  m_walk.push_back(0.0);
}

void RandomWalk::gen_steps(unsigned int n)
{
  for (int i=0; i<n; ++i) {
    m_walk.push_back(m_walk.back() + m_gen_next_incr());
  }
}

const deque<float>& RandomWalk::get_walk() {
  return m_walk;
}

class UnifFunctor {
private:
  std::mt19937 m_gen;
  std::uniform_real_distribution<> m_dis;

public:
  UnifFunctor() 
  {
    std::random_device rd;
    m_gen = std::mt19937(rd());
    m_dis = std::uniform_real_distribution<>(-1.0, 1.0); 
  }
  
  float operator () ()
  {
    return m_dis(m_gen);
  }
};
class NormalFunctor {
private:
  std::mt19937 m_gen;
  std::normal_distribution<> m_dis;

public:
  NormalFunctor() 
  {
    std::random_device rd;
    m_gen = std::mt19937(rd());
    m_dis = std::normal_distribution<>(0.0, 1.0); 
  }
  
  float operator () ()
  {
    return m_dis(m_gen);
  }
};
class CauchyFunctor {
private:
  std::mt19937 m_gen;
  std::cauchy_distribution<> m_dis;

public:
  CauchyFunctor() 
  {
    std::random_device rd;
    m_gen = std::mt19937(rd());
    m_dis = std::cauchy_distribution<>(0.0, 1.0); 
  }
  
  float operator () ()
  {
    return m_dis(m_gen);
  }
};


