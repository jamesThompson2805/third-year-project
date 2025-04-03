#include "capla.h"

/**
 * @file capla.cpp
 * @brief File implementing functions of capla.h
 */


#include "pla.h"
#include "conv_double_window.h"

DRT capla_eval::generate_mean_DRT(unsigned int win_size)
{
  Seqd l(win_size, 1/(double)win_size);
  Seqd r(win_size, 1/(double)win_size);
  return [l, r](const Seqd& s, unsigned int num_params){ return pla::apla_to_seq( c_d_w::conv_pla(s, num_params, l, r) ); };
}
DRT_COMPR capla_eval::generate_mean_DRT_COMPR(unsigned int win_size)
{
  Seqd l(win_size, 1/(double)win_size);
  Seqd r(win_size, 1/(double)win_size);
  return [l, r](const Seqd& s, unsigned int num_params){ return c_d_w::conv_pla(s, num_params, l, r); };
}

DRT capla_eval::generate_mean_skip_one_DRT(unsigned int win_size)
{
  Seqd l(win_size, 1/(double)win_size);
  Seqd r(win_size, 1/(double) (win_size-1) );
  r[0] = 0.0;
  return [l, r](const Seqd& s, unsigned int num_params){ return pla::apla_to_seq( c_d_w::conv_pla(s, num_params, l, r) ); };
}

DRT capla_eval::generate_tri_DRT(unsigned int win_size)
{
  Seqd l(win_size);
  Seqd r(win_size);
  for (int i=0; i<win_size; i++) {
    l[i] = r[win_size -1 -i] = (double) 2*(i+1) / (double) (win_size * (win_size + 1));
  }
  return [l, r](const Seqd& s, unsigned int num_params){ return pla::apla_to_seq( c_d_w::conv_pla(s, num_params, l, r) ); };
}

DRT capla_eval::generate_tri_skip_one_DRT(unsigned int win_size)
{
  Seqd l(win_size);
  Seqd r(win_size);
  for (int i=0; i<win_size; i++) {
    l[i] = (double) 2*(i+1) / (double) (win_size * (win_size + 1));
    r[win_size -1 -i] = (double) 2*(i+1) / (double) (win_size * (win_size - 1));
  }
  r[0] = 0.0;
  return [l, r](const Seqd& s, unsigned int num_params){ return pla::apla_to_seq( c_d_w::conv_pla(s, num_params, l, r) ); };
}
