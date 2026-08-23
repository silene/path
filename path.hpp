#ifndef PATH_HPP
#define PATH_HPP

#include "image.hpp"

struct integrator {
  Image::buffer img;
  integrator(int w, int h);
  void pixel(int x, int y);
  void tiled();
};

namespace Solver {

void prepare();

}

#endif
