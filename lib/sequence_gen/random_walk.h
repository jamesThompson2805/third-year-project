#ifndef RANDOM_WALK_H
#define RANDOM_WALK_H

#include <random>
#include <functional>
#include <deque>
#include <string>
using std::deque;

/**
 * @file random_walk.h
 * @brief Header file the RandomWalk class and various functors to use inside RandomWalk.
 */


/**
 * Models a random walk of independent, identically distributed random variables.
 */
class RandomWalk {
private:
  deque<double> m_walk;
  std::function< double() > m_gen_next_incr;
  unsigned int num_steps;

public:
  /**
   * @brief Construct a random walk using a function that generates the next number when called.
   * @param series is a mutable reference to a sequence.
   */
  RandomWalk(std::function<double()> next_incr_generator);
  /**
   * @brief gen_steps simulates n steps of the random walk.
   * @param n The number of steps to simulate.
   */
  void gen_steps(unsigned int n);
  /**
   * @brief get_walk returns a view into the random walk
   * @return A const reference to a double ended queue representing the walk so far.
   */
  const deque<double>& get_walk();
  /**
   * @brief save_walk saves the current walk so far in TSV format.
   * @param filename The string to be opened as a file to save the walk in, should be of format "file_path/name.format".
   */
  void save_walk(std::string filename);
};

/**
 * @brief Functor to simulate independent trials of Radermacher variable.
 */
class BernFunctor{
private:
  std::mt19937 m_gen;
  std::bernoulli_distribution m_dis;

public:
  /**
   * @brief Constructs the Functor based off the parameters.
   * @param seed The seed for the pseudo random generator.
   * @param p The probability of the variable being +1.
   */
  BernFunctor(unsigned int seed=0, double p=0.5) 
  {
    std::random_device rd;
    m_gen = std::mt19937(rd());
    m_gen.seed(seed);
    m_dis = std::bernoulli_distribution(p); 
  }
  
  /**
   * @brief Overloaded () operator allowing Functor to act like function.
   * @return result of independent trial of variable.
   */
  double operator () ()
  {
    return 2*m_dis(m_gen) - 1.0;
  }
};
/**
 * @brief Functor to simulate independent trials of Uniform variable.
 */
class UnifFunctor {
private:
  std::mt19937 m_gen;
  std::uniform_real_distribution<> m_dis;

public:
  /**
   * @brief Constructs the Functor based off the parameters.
   * @param seed The seed for the pseudo random generator.
   * @param min The minimum value of the interval.
   * @param max The maximum value of the interval.
   */
  UnifFunctor(unsigned int seed=0, double min=-1.0, double max=1.0) 
  {
    std::random_device rd;
    m_gen = std::mt19937(rd());
    m_gen.seed(seed);
    m_dis = std::uniform_real_distribution<>(min, max); 
  }
  
  /**
   * @brief Overloaded () operator allowing Functor to act like function.
   * @return result of independent trial of variable.
   */
  double operator () ()
  {
    return m_dis(m_gen);
  }
};
/**
 * @brief Functor to simulate independent trials of Normal variable.
 */
class NormalFunctor {
private:
  std::mt19937 m_gen;
  std::normal_distribution<> m_dis;

public:
  /**
   * @brief Constructs the Functor based off the parameters.
   * @param seed The seed for the pseudo random generator.
   * @param mean The average value of the trial.
   * @param var The variance of the normal random variable.
   */
  NormalFunctor(unsigned int seed=0, double mean=0.0, double var=1.0) 
  {
    std::random_device rd;
    m_gen = std::mt19937(rd());
    m_gen.seed(seed);
    m_dis = std::normal_distribution<>(mean, var); 
  }
  
  /**
   * @brief Overloaded () operator allowing Functor to act like function.
   * @return result of independent trial of variable.
   */
  double operator () ()
  {
    return m_dis(m_gen);
  }
};
/**
 * @brief Functor to simulate independent trials of Cauchy variable.
 */
class CauchyFunctor {
private:
  std::mt19937 m_gen;
  std::cauchy_distribution<> m_dis;

public:
  /**
   * @brief Constructs the Functor based off the parameters.
   * @param seed The seed for the pseudo random generator.
   * @param a Parameter a of Cauchy distribution.
   * @param b Parameter b of Cauchy distribution.
   */
  CauchyFunctor(unsigned int seed=0, double a=0.0, double b=1.0) 
  {
    std::random_device rd;
    m_gen = std::mt19937(rd());
    m_gen.seed(seed);
    m_dis = std::cauchy_distribution<>(a, b); 
  }
  
  /**
   * @brief Overloaded () operator allowing Functor to act like function.
   * @return result of independent trial of variable.
   */
  double operator () ()
  {
    return m_dis(m_gen);
  }
};

#endif
