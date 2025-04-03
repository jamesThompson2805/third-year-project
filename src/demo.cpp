/**
 * @file demo.cpp
 *
 * @brief This file offers a demo of some of the capabilities of the project designed for the presentation.
 */
#include "random_walk.h"
#include "ucr_parsing.h"
#include "z_norm.h"

#include "bottom_up.h"

#include "r_tree.h"
#include "lower_bounds_apla.h"

#include "plotting/series_plotting.h"
#include "plotting/plot_dimreduct_pla.h"
#include "evaluations/capla.h"

#include <vector>
#include <string>
using std::vector, std::string;

/**
 * @brief The demo function automatically runs all sections of the demo in order, displaying a full end to end demonstration. It requires firefox to run.
 */
void demo()
{
  /**** GENERATE RANDOM WALK **********/
  RandomWalk walk( NormalFunctor(11) ); 
  walk.gen_steps(100'000);
  walk.save_walk("./tsv/testwalk1.tsv");
  vector<double> synth_data = ucr_parsing::parse_tsv("./tsv/testwalk1.tsv",-1);
  z_norm::z_normalise(synth_data);
  vector<double> short_synth( synth_data.begin(), synth_data.begin()+2000);

  Series s_walk = { synth_data, "Normal Walk" };
  Series s_walk_short = { short_synth, "Normal Walk" };
  PlotDetails p_walk = { "Plot of generated random walk", "Time", "Value", "./img/", PDF };
  PlotDetails p_walk_short = { "Plot of subsequence of generated random walk", "Time", "Value", "./img/", PDF };
  plot::plot_series(s_walk, p_walk);
  plot::plot_series(s_walk_short, p_walk_short);

  /**** LOAD REAL DATA ****************/
  string ucr_datasets_loc = "external/data/UCRArchive_2018/";
  vector<string> datasets = ucr_parsing::parse_folder_names(ucr_datasets_loc);

  unsigned int di = 30;
  vector<double> real_data = parse_ucr_dataset(datasets[di], ucr_datasets_loc,  ucr_parsing::DatasetType::TEST);
  z_norm::z_normalise(real_data);
  vector<double> short_real( real_data.begin(), real_data.begin()+2000);

  Series s_real = { real_data, "ECG Scan" };
  Series s_real_short = { short_real, "ECG Scan" };
  PlotDetails p_real = { "Plot of ECG Data", "Time", "Voltage", "./img/", PDF };
  PlotDetails p_real_short = { "Plot of subsequence of ECG Data", "Time", "Voltage", "./img/", PDF };
  plot::plot_series(s_real, p_real);
  plot::plot_series(s_real_short, p_real_short);

  /***** TEST SOME DRT'S ************/
  auto capla_tri_f = capla_eval::generate_tri_DRT(5);
  auto capla_drt = [](const Seqd& s){ return capla_eval::generate_mean_DRT_COMPR(5)(s, 33); };
  auto bottom_up_f_uncompr = [&](const Seqd& s, unsigned int parameter){ 
    auto apla = bottom_up::bottom_up(s, 0.1, bottom_up::se);
    if (apla.size() < parameter/3) segmerge::segment_to_dim(s,apla,parameter);
    if (apla.size() > parameter/3) segmerge::merge_to_dim(s,apla,parameter);
    return apla;
  };
  auto b_u_drt = [&bottom_up_f_uncompr](const Seqd& s){ return bottom_up_f_uncompr(s, 63); };

  PlotDetails p_capla_drt = { "Plot of Sliding Window upon ECG", "Time", "Voltage", "./img/", PDF };
  //plot_pla::plot_any_apla(short_real, "ECG Scan", capla_drt, "Sliding Window", p_capla_drt);
  PlotDetails p_b_u_drt = { "Plot of Bottom Up Approach upon ECG", "Time", "Voltage", "./img/", PDF };
  //plot_pla::plot_any_apla(short_real, "ECG Scan", b_u_drt, "Bottom Up", p_capla_drt);

  /****** FIND SOME HEART-BEATS *****/
  const unsigned int seq_size = 50;
  const unsigned int NS = 5;
  auto retrieval_f = [](const unsigned int& i, const vector<double>& q) {
    return std::vector<std::array<const double*,2>>( {{ q.data()+i, q.data()+i+seq_size-1 }} );
  };
  RTree<apla_bounds::AplaMBR<NS>, unsigned int> r_tree(40,10 , apla_bounds::mbr_area<NS> , apla_bounds::mbr_merge<NS> , apla_bounds::dist_to_mbr_sqr<NS>);
  auto vec_of_mbrs = apla_bounds::vec_to_subseq_mbrs<NS>(short_real,seq_size,bottom_up_f_uncompr);
  for (int i=0; i<vec_of_mbrs.size(); i++) {
    r_tree.insert(vec_of_mbrs[i], i);
  }
  Seqd query(short_real.begin() + 51, short_real.begin() + 101);
  Seqd short_real_x; for (int i=0; i<short_real.size(); i++) short_real_x.push_back(i);
  std::vector<Line> vl;
  Line short_real_line = { short_real_x, short_real, "ECG Data" };
  vl.push_back(short_real_line);
  auto sim_results = r_tree.sim_search_exact(query, 7.0, retrieval_f, short_real);
  for (const auto& [l_ptr,r_ptr] : sim_results) {
    std::cout << l_ptr - short_real.data() << std::endl;
    Seqd res(l_ptr, r_ptr+1);
    Seqd res_x; for (int i=0; i<res.size(); i++) res_x.push_back( l_ptr - short_real.data() + i );
    vl.push_back({res_x, res, "distance 2 away"});
  }
  PlotDetails p_sim = { "Plot of closest to first heart beat upon ECG", "Time", "Voltage", "./img/", PDF };
  plot::plot_lines(vl,p_sim);

  /************ Find Some Similar Subsequences ***********/
  RTree<apla_bounds::AplaMBR<NS>, unsigned int> r_tree2(40,10 , apla_bounds::mbr_area<NS> , apla_bounds::mbr_merge<NS> , apla_bounds::dist_to_mbr_sqr<NS>);
  auto vec_of_mbrs_synth = apla_bounds::vec_to_subseq_mbrs<NS>(short_synth,seq_size,bottom_up_f_uncompr);
  for (int i=0; i<vec_of_mbrs_synth.size(); i++) {
    r_tree2.insert(vec_of_mbrs_synth[i], i);
  }
  Seqd query_synth(short_synth.begin(), short_synth.begin() + 50);
  Seqd short_synth_x; for (int i=0; i<short_synth.size(); i++) short_synth_x.push_back(i);
  vl.clear();
  Line short_synth_line = { short_synth_x, short_synth, "Synthetic Data" };
  vl.push_back(short_synth_line);
  auto knn_results = r_tree2.knn_search(query_synth, 8, retrieval_f, short_synth);
  std::cout << std::endl;
  for (const auto& [l_ptr,r_ptr] : knn_results) {
    std::cout << l_ptr - short_synth.data() << std::endl;
    Seqd res(l_ptr, r_ptr+1);
    Seqd res_x; for (int i=0; i<res.size(); i++) res_x.push_back( l_ptr - short_synth.data() + i );
    vl.push_back({res_x, res, "k-nn close"});
  }
  PlotDetails p_knn = { "Plot of 8 closest to beginning", "Time", "Value", "./img/", PDF };
  plot::plot_lines(vl,p_knn);

}
