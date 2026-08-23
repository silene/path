#ifndef CAMERA_HPP
#define CAMERA_HPP

#include "base.hpp"

namespace Camera {

namespace {
using Vector::vec;
using Matrix::mat;
}

struct simple: base {
  vec pos;
  mat rot;
  double zoom;
  simple(vec const &p, mat const &r, double z);
  simple(vec const &p, vec const &t, double z, vec const &y);
  simple(vec const &p, vec const &t, double z);
  std::pair<vec, vec> get(double fx, double fy) const;
};

struct simple_lens: simple {
  double lens, focal;
  simple_lens(vec const &p, mat const &r, double z, double l, double f);
  simple_lens(vec const &p, vec const &t, double z, double l, vec const &y);
  simple_lens(vec const &p, vec const &t, double z, double l);
  std::pair<vec, vec> get(double fx, double fy) const;
};

}

#endif
