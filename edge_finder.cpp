# include "edge_finder.h"

constexpr uint8_t HALO = 1;

/******* Edge_finder ********/
// TODO: expects rho including boundaries!

Edge_finder::Edge_finder(Lattice_accessor& geometry, int threshold)
  : threshold{threshold}, geometry{geometry}
{
}

Edge_finder::~Edge_finder() {
}

int Edge_finder::detect_edges(const vector<double>& rho, size_t threshold) {
  edges.clear();
  //vector<double> temp((MX-3)*(MY-3)*(MZ-3));
  //temp = gaussian_blur(component[0]->rho);
  edges = sobel_edge_detector(rho, threshold);
  return 0;
}

vector<double> Edge_finder::sobel_edge_detector(const vector<double>& rho, size_t tolerance) {
  //TODO: Generalize to 2D and 1D or add warning?
  vector<double> result( rho.size() );
  size_t threshold = tolerance;

  vector<double> Gx_minus = {-1, -3, -1, -3, -6, -3, -1, -3, -1};
  vector<double> Gx_mid = {0, 0, 0, 0, 0, 0, 0, 0, 0};
  vector<double> Gx_plus = {1, 3, 1, 3, 6, 3, 1, 3, 1};

  vector<double> Gy_minus = {-1, 0, 1, -3, 0, 3, -1, 0, 1};
  vector<double> Gy_mid = {-3, 0, 3, -6, 0, 6, -3, 0, 3};
  vector<double> Gy_plus = {-1, 0, 1, -3, 0, 3, -1, 0, 1};

  vector<double> Gz_minus = {1, 3, 1, 0, 0, 0, -1, -3, -1};
  vector<double> Gz_mid = {3, 6, 3, 0, 0, 0, -3, -6, -3};
  vector<double> Gz_plus = {1, 3, 1, 0, 0, 0, -1, -3, -1};

  int i{0};

  for (size_t x = 0; x < geometry.MX; ++x)
    for (size_t y = 0; y < geometry.MY; ++y)
      for (size_t z = 0 ; z < geometry.MZ; ++z ) {

        double conv_x = convolution(Gx_minus, get_xy_plane(rho, x, y, z));
        conv_x += convolution(Gx_mid, get_xy_plane(rho, x, y, z+1));
        conv_x += convolution(Gx_plus, get_xy_plane(rho, x, y, z+2));

        double conv_y = convolution(Gy_minus, get_xy_plane(rho, x, y, z));
        conv_y += convolution(Gy_mid, get_xy_plane(rho, x, y, z+1));
        conv_y += convolution(Gy_plus, get_xy_plane(rho, x, y, z+2));

        double conv_z = convolution(Gz_minus, get_xz_plane(rho, x, y, z));
        conv_z += convolution(Gz_mid, get_xz_plane(rho, x, y+1, z));
        conv_z += convolution(Gz_plus, get_xz_plane(rho, x, y+2, z));

        *geometry.val_ptr(result, x+1, y+1, z+1) = abs(conv_x) + abs(conv_y) + abs(conv_z);

        ++i;
    }

  //normalize between 0 and 255
  double min = *min_element(result.begin(), result.end());
  double max = *max_element(result.begin(), result.end());

  for (double& all_elements : result)
    all_elements = (255 - 0) * ((all_elements - min) / (max - min)) + 0;

  //cut-off at threshold
  for (double& all_elements : result)
    if (all_elements < threshold) {
      all_elements = 0;
    }

  return result;
}

vector<double> Edge_finder::gaussian_blur(const vector<double>& rho) {
  vector<double> result( (geometry.MX-HALO)*(geometry.MY-HALO)*(geometry.MY-HALO) );
  vector<double> G1 = {1, 2, 1, 2, 4, 2, 1, 2, 1};
  vector<double> G2 = {1, 1, 1, 1, 2, 1, 1, 1, 1};
  vector<double> G3 = {1, 1, 1, 1, 2, 1, 1, 1, 1};

  for (double& all_values : G1)
    all_values = all_values * 1/16;

  for (double& all_values : G2)
    all_values = all_values * 1/16;

  for (double& all_values : G3)
    all_values = all_values * 1/16;

  int i{0};

  for (size_t x = 0; x < geometry.MX - HALO; ++x)
    for (size_t y = 0; y < geometry.MY - HALO; ++y)
      for (size_t z = 0 ; z < geometry.MZ - HALO ; ++z ) {
        double conv_1 = convolution(G1, get_xy_plane(rho, x, y, z));
        double conv_2 = convolution(G2, get_xy_plane(rho, x, y, z+1));
        double conv_3 = convolution(G3, get_xy_plane(rho, x, y, z+2));
        result[i] = (conv_1 + conv_2 + conv_3);
        ++i;
    }
  return result;
}

double Edge_finder::convolution(const vector<double>& kernel, const vector<double>& pixel) {
  if (kernel.size() != pixel.size()) {
    std::cerr << "Convolution: pixel and kernel not of equal size!" << std::endl;
  }

  double accumulator = 0;

  for (size_t i = 0 ; i < kernel.size() ; ++i) {
    accumulator += kernel[i] * pixel[i];
  }

  return accumulator;
}

vector<double> Edge_finder::get_xy_plane(const vector<double>& rho, int x, int y, int z, int size) {
  vector<double> pixel(size*size);

  int i = 0;

  for (int horizontal = 0 ; horizontal < size ; ++horizontal)
    for (int vertical = 0 ; vertical < size ; ++vertical) {
        pixel[i] = geometry.val(const_cast<vector<double>&>(rho), x+horizontal, y+vertical, z);
        ++i;
    }

  return pixel;
}

vector<double> Edge_finder::get_xz_plane(const vector<double>& rho, int x, int y, int z, int size) {
  vector<double> pixel(size*size);

  int i = 0;

  for (int horizontal = 0 ; horizontal < size ; ++horizontal)
    for (int depth = 0 ; depth < size ; ++depth) {
        pixel[i] = geometry.val(const_cast<vector<double>&>(rho), x+horizontal, y, z+depth);
        ++i;
    }

  return pixel;
}
