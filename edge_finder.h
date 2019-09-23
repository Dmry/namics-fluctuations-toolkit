#ifndef edge_finder_h
#define edge_finder_h

#include <vector>
#include "lattice_accessor.h"
#include <algorithm>
#include <iostream>

using namespace std;

class Edge_finder {
public:
  Edge_finder(Lattice_accessor&, int threshold);
  ~Edge_finder();
  int detect_edges(const vector<double>&, size_t);

  vector<double> edges;
  const int threshold;

private:
  Lattice_accessor geometry;
  vector<double> sobel_edge_detector(const vector<double>&, size_t);
  vector<double> gaussian_blur(const vector<double>&);
  double convolution(const vector<double>&, const vector<double>&);
  vector<double> get_xy_plane(const vector<double>&, int, int, int, int = 3);
  vector<double> get_xz_plane(const vector<double>&, int, int, int, int = 3);
};

#endif
