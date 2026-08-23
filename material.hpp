#ifndef MATERIAL_HPP
#define MATERIAL_HPP

#include "base.hpp"
#include "linalg.hpp"
#include "color.hpp"

namespace Material {

using Vector::vec;
using Spectrum::sampled_wl;

struct emissive: base {
  Light::ptr light;

  emissive(Light::ptr l)
    : base(Emissive, false, false), light(l) {}

  biased<ray> sample(interaction const &, sampled_wl const &wl) const;
};

struct lambertian: base {
  Spectrum::ptr sp;
  lambertian(Spectrum::ptr s)
    : base(Solid, true, true), sp(s) {}

  sampled_spectrum bxdf(interaction const &pt, vec const &inc, sampled_wl const &wl) const;
  double pdf(interaction const &pt, vec const &inc) const;
  biased<ray> sample(interaction const &pt, sampled_wl const &wl) const;
};

struct reflective: base {
  Spectrum::ptr eta, extinct;
  reflective(Spectrum::ptr n, Spectrum::ptr k)
    : base(Solid, false, false), eta(n), extinct(k) {}

  biased<ray> sample(interaction const &pt, sampled_wl const &wl) const;

#if 0
  base const *regularize() const {
    return new rough(sp, 0.4);
  }
#endif
};

struct refractive: base {
  double eta;
  refractive(double e): base(Transmitive, false, false), eta(e) {}

  biased<ray> sample(interaction const &pt, sampled_wl const &wl) const;
};

struct thin_refractive: base {
  double eta;
  base const *mat;
  thin_refractive(double e, base const *m);
  sampled_spectrum bxdf(interaction const &pt, vec const &inc, sampled_wl const &wl) const;
  double pdf(interaction const &pt, vec const &inc) const;
  biased<ray> sample(interaction const &pt, sampled_wl const &wl) const;
};

using ptr = Material::base const *;

}

#endif
