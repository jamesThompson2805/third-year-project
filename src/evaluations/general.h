#ifndef EVAL_GENERAL_H
#define EVAL_GENERAL_H

#include <functional>
#include <vector>

typedef std::vector<double> Seqd;
typedef std::vector<unsigned int> Sequi;

typedef std::function<Seqd(const Seqd&, unsigned int)> DRT;
typedef std::function<double(const Seqd&, const Seqd&)> COMP;


namespace general_eval {
  double comp_of_method(const Seqd& s, unsigned int num_params, DRT f, COMP comp);
  double mse_of_method(const Seqd& s, unsigned int num_params, DRT f);
  double maxdev_of_method(const Seqd& s, unsigned int num_params, DRT f);
  double cputime_ms_of_method(const Seqd& s, unsigned int num_params, DRT f);

  Seqd get_comp_over_sizes(const Seqd& s, Sequi sizes, unsigned int num_params, DRT f, COMP comp);
  Seqd get_mse_over_sizes(const Seqd& s, Sequi sizes, unsigned int num_params, DRT f);
  Seqd get_maxdev_over_sizes(const Seqd& s, Sequi sizes, unsigned int num_params, DRT f);
  Seqd get_cputime_over_sizes(const Seqd& s, Sequi sizes, unsigned int num_params, DRT f);

  Seqd get_comp_over_num_params(const Seqd& s, Sequi vec_params, DRT f, COMP comp);
  Seqd get_mse_over_num_params(const Seqd& s, Sequi vec_params, DRT f);
  Seqd get_maxdev_over_num_params(const Seqd& s, Sequi vec_params, DRT f);
  Seqd get_cputime_over_num_params(const Seqd& s, Sequi vec_params, DRT f);
}


#endif
