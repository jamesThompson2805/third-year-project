#include "swing.h"
using std::vector, std::tuple;
using std::min, std::max;

Seqddt swing::swing(const Seqd& s, double epsilon)
{
  if (epsilon < 0) return {};
  if (s.size() < 1) return {};
  if (s.size() == 1) return {{{s[0],0.0},0}};
  if (s.size() == 2) return {{{s[0],s[1]-s[0]},1}};

  Seqddt breakpoints;
  unsigned int seg_start = 0;
  unsigned int c_ind = 2;
  DoublePair u = {s[0], s[1]-s[0]+epsilon};
  DoublePair l = {s[0], s[1]-s[0]-epsilon};
  while (c_ind < s.size()) {
    u = {u[0], min( (s[c_ind]-u[0]+epsilon)/(c_ind-seg_start), u[1])};
    l = {l[0], max( (s[c_ind]-l[0]-epsilon)/(c_ind-seg_start), l[1])};


    if (u[1] >= l[1]) {
      // good case, still within epsilon of all points
      c_ind++;
    } else {
      // not still within epsilon of all points => last in segment was index before this one
      // => find the line of best fit for previous
      // 	and as lines connected we start next segment from previous index
      double best_grad = pla::regression_thru_point(u[0] ,s.data() + seg_start+1, s.data() + c_ind - 1);
      if (breakpoints.size() == 0) {
	// only  the first line segment actually 'starts' at seg-start
	breakpoints.push_back( {{u[0],best_grad}, c_ind-1} );
      } else {
	breakpoints.push_back( {{u[0]+best_grad,best_grad}, c_ind-1} );
      }

      // calculate the next upper and lower lines via projected start of prev line
      double l_start = u[0] + best_grad*(c_ind-1-seg_start);
      u = {l_start, (s[c_ind]+epsilon-l_start)/(1.0) };
      l = {l_start, (s[c_ind]-epsilon-l_start)/(1.0) };
      // adjust the segment start and final increment c_ind
      seg_start = c_ind - 1;
      c_ind++;
    }
  }
  // algorithm doesn't record final line as terminates beforehand so manually add it
  double best_grad = pla::regression_thru_point(u[0], s.data() + seg_start + 1, s.data() + c_ind - 1);
  breakpoints.push_back( {{u[0]+best_grad,best_grad}, c_ind-1} );

  return breakpoints;
}

vector<tuple<double, unsigned int>> swing::swing_compr(const Seqd& s, double epsilon)
{
  if (epsilon < 0) return {};
  if (s.size() < 1) return {};
  if (s.size() == 1) return {{0.0,0}};
  if (s.size() == 2) return {{s[1]-s[0],1}};

  vector<tuple<double, unsigned int>> breakpoints;
  unsigned int seg_start = 0;
  unsigned int c_ind = 2;
  DoublePair u = {s[0], s[1]-s[0]+epsilon};
  DoublePair l = {s[0], s[1]-s[0]-epsilon};
  while (c_ind < s.size()) {
    u = {u[0], min( (s[c_ind]-u[0]+epsilon)/(c_ind-seg_start), u[1])};
    l = {l[0], max( (s[c_ind]-l[0]-epsilon)/(c_ind-seg_start), l[1])};


    if (u[1] >= l[1]) {
      // good case, still within epsilon of all points
      c_ind++;
    } else {
      // not still within epsilon of all points => last in segment was index before this one
      // => find the line of best fit for previous
      // 	and as lines connected we start next segment from previous index
      double best_grad = pla::regression_thru_point(u[0] ,s.data() + seg_start+1, s.data() + c_ind - 1);
      if (breakpoints.size() == 0) {
	// this is the only line segment that actually 'starts' at seg-start
	breakpoints.push_back( {best_grad, c_ind-1} );
      } else {
	breakpoints.push_back( {best_grad, c_ind-1} );
      }

      // calculate the next upper and lower lines via projected start of prev line
      double l_start = u[0] + best_grad*(c_ind-1-seg_start);
      u = {l_start, (s[c_ind]-l_start+epsilon)/(c_ind-seg_start) };
      l = {l_start, (s[c_ind]-l_start-epsilon)/(c_ind-seg_start) };
      // adjust the segment start and final increment c_ind
      seg_start = c_ind - 1;
      c_ind++;
    }
  }
  // algorithm doesn't record final line as terminates beforehand so manually add it
  double best_grad = pla::regression_thru_point(u[0], s.data() + seg_start + 1, s.data() + c_ind - 1);
  breakpoints.push_back( {best_grad, c_ind-1} );

  return breakpoints;
}

Seqddt swing::slide(const Seqd &s, double epsilon)
{
  if (epsilon < 0) return {};
  if (s.size() < 1) return {};
  if (s.size() == 1) return {{{s[0],0.0},0}};
  if (s.size() == 2) return {{{s[0],s[1]-s[0]},1}};

  Seqddt breakpoints;
  unsigned int seg_start = 0;
  unsigned int c_ind = 2;
  DoublePair u = {s[0]-epsilon, s[1]-s[0]+2*epsilon};
  DoublePair l = {s[0]+epsilon, s[1]-s[0]-2*epsilon};
  while (c_ind < s.size()) {
    u = {u[0], min( (s[c_ind]-u[0]+2*epsilon)/(c_ind-seg_start), u[1])};
    l = {l[0], max( (s[c_ind]-l[0]-2*epsilon)/(c_ind-seg_start), l[1])};


    if (u[1] >= l[1]) {
      // good case, still within epsilon of all points
      c_ind++;
    } else {
      // not still within epsilon of all points => last in segment was index before this one
      // => find the line of best fit for previous
      // 	and as lines connected we start next segment from previous index
      double best_grad = pla::regression_thru_point(u[0] ,s.data() + seg_start+1, s.data() + c_ind - 1);
      if (breakpoints.size() == 0) {
	// this is the only line segment that actually 'starts' at seg-start
	breakpoints.push_back( {{u[0],best_grad}, c_ind-1} );
      } else {
	breakpoints.push_back( {{u[0]+best_grad,best_grad}, c_ind-1} );
      }

      // calculate the next upper and lower lines via projected start of prev line
      double l_start = u[0] + best_grad*(c_ind-1-seg_start);
      u = {l_start, (s[c_ind]-l_start+epsilon)/(c_ind-seg_start) };
      l = {l_start, (s[c_ind]-l_start-epsilon)/(c_ind-seg_start) };
      // adjust the segment start and final increment c_ind
      seg_start = c_ind - 1;
      c_ind++;
    }
  }
  // algorithm doesn't record final line as terminates beforehand so manually add it
  double best_grad = pla::regression_thru_point(u[0], s.data() + seg_start + 1, s.data() + c_ind - 1);
  breakpoints.push_back( {{u[0]+best_grad,best_grad}, c_ind-1} );

  return breakpoints;
}

