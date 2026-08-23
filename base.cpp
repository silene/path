#include <cassert>

#include "base.hpp"

namespace Spectrum {

sampled_wl::sampled_wl(double s) {
  for (int i = 0; i < nb_s; ++i) {
    double l = s + i * (double)(max_wl - min_wl) / nb_s;
    if (l >= max_wl) l = l - max_wl + min_wl;
    assert(min_wl <= l && l < max_wl);
    lambda[i] = l;
    pdf[i] = 1. / (max_wl - min_wl);
  }
}

}

namespace Material {

Spectrum::sampled_spectrum base::bxdf(interaction const &, vec const &inc, sampled_wl const &) const {
  assert(false);
}

double base::pdf(interaction const &, vec const &inc) const {
  assert(false);
}

}

namespace Light {

double base::pdf(vec const &pos, vec const &n, vec const &dir) const {
  assert(false);
}

}
