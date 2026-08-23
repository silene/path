#ifndef SAMPLER_HPP
#define SAMPLER_HPP

#include <random>
#include <vector>

#include "linalg.hpp"

extern thread_local std::mt19937_64 *rng;

template<class T>
using biased = std::pair<T, double>;

// Computations should be invariant wrt this arbitrarily large value.
// Taken as 1 to avoid numerical issues.
double constexpr Dirac = 1.;

namespace Sampler {

using Vector::vec;

struct discrete {
  std::vector<float> probas;
  std::vector<std::pair<float, int>> aliases;
  double sum;
  discrete() = default;
  discrete(int nb, float const *);
  biased<int> sample() const;
  double pdf(int i) const { return sum ? probas[i] : 0.; }
};

struct discrete2D {
  discrete data;
  std::vector<discrete> rows;
  discrete2D() = default;
  discrete2D(int w, int h, float const *);
  biased<std::pair<int,int>> sample() const;
  double pdf(int x, int y) const;
};

point2 disk_uniform();

struct sphere_uniform {
  biased<vec> sample() const;
  double pdf(vec const &) const;
};

struct hemisphere_uniform {
  vec u;
  hemisphere_uniform(vec const &u_): u(u_) {}
  biased<vec> sample() const;
  double pdf(vec const &) const;
};

struct hemisphere_linear {
  vec u;
  hemisphere_linear(vec const &u_): u(u_) {}
  biased<vec> sample() const;
  double pdf(vec const &) const;
};

struct hemisphere_power {
  vec u;
  double p;
  hemisphere_power(vec const &u_, double p_): u(u_), p(p_) {}
  biased<vec> sample() const;
  double pdf(vec const &) const;
};

struct cone_uniform {
  vec u;
  double cmax;
  cone_uniform(vec const &u_, double c_): u(u_), cmax(c_) {}
  biased<vec> sample() const;
  double pdf(vec const &) const;
};

}

#endif
