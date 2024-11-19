#ifndef RANDOM_WALK_H
#define RANDOM_WALK_H

#include <random>
#include <functional>
#include <deque>
#include <string>
using std::deque;


class RandomWalk {
private:
  deque<double> m_walk;
  std::function< double() > m_gen_next_incr;

public:
  RandomWalk(std::function<double()> next_incr_generator);
  void gen_steps(unsigned int n);
  const deque<double>& get_walk();
};

class BernFunctor{
private:
  std::mt19937 m_gen;
  std::bernoulli_distribution m_dis;

public:
  BernFunctor(unsigned int seed=0, double p=0.5) 
  {
    std::random_device rd;
    m_gen = std::mt19937(rd());
    m_gen.seed(seed);
    m_dis = std::bernoulli_distribution(p); 
  }
  
  double operator () ()
  {
    return 2*m_dis(m_gen) - 1.0;
  }
};
class UnifFunctor {
private:
  std::mt19937 m_gen;
  std::uniform_real_distribution<> m_dis;

public:
  UnifFunctor(unsigned int seed=0, double min=-1.0, double max=1.0) 
  {
    std::random_device rd;
    m_gen = std::mt19937(rd());
    m_gen.seed(seed);
    m_dis = std::uniform_real_distribution<>(min, max); 
  }
  
  double operator () ()
  {
    return m_dis(m_gen);
  }
};
class NormalFunctor {
private:
  std::mt19937 m_gen;
  std::normal_distribution<> m_dis;

public:
  NormalFunctor(unsigned int seed=0, double mean=0.0, double var=1.0) 
  {
    std::random_device rd;
    m_gen = std::mt19937(rd());
    m_gen.seed(seed);
    m_dis = std::normal_distribution<>(mean, var); 
  }
  
  double operator () ()
  {
    return m_dis(m_gen);
  }
};
class CauchyFunctor {
private:
  std::mt19937 m_gen;
  std::cauchy_distribution<> m_dis;

public:
  CauchyFunctor(unsigned int seed=0, double mean=0.0, double var=1.0) 
  {
    std::random_device rd;
    m_gen = std::mt19937(rd());
    m_gen.seed(seed);
    m_dis = std::cauchy_distribution<>(mean, var); 
  }
  
  double operator () ()
  {
    return m_dis(m_gen);
  }
};

#endif
