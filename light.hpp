#ifndef LIGHT_HPP
#define LIGHT_HPP

#include "base.hpp"
#include "image.hpp"
#include "solid.hpp"

namespace Light {

struct point: base {
  Spectrum::ptr sp;
  vec pos;

  point(Spectrum::ptr s, vec const &p)
    : base(false, false), sp(s), pos(p) {}

  biased<ray> sample(vec const &p, vec const &) const;
  sampled_spectrum get_sp(vec const &p, vec const &, sampled_wl const &wl) const;
};

struct spot: base {
  Spectrum::ptr sp;
  vec pos, dir;
  double angle1, angle2;

  spot(Spectrum::ptr s, vec const &p, vec const &d, double a1, double a2)
    : base(false, false), sp(s), pos(p), dir(d), angle1(a1), angle2(a2) {}

  biased<ray> sample(vec const &p, vec const &) const;
  sampled_spectrum get_sp(vec const &p, vec const &, sampled_wl const &wl) const;
};

struct directional: base {
  Spectrum::ptr sp;
  vec dir;

  directional(Spectrum::ptr s, vec const &d)
    : base(false, false), sp(s), dir(d) {}

  biased<ray> sample(vec const &, vec const &) const;
  sampled_spectrum get_sp(vec const &, vec const &, sampled_wl const &wl) const;
};

struct multidirectional: base {
  Spectrum::ptr sp;
  vec dir;
  double cmax, inv_area;

  multidirectional(Spectrum::ptr s, vec const &d, double a);
  biased<ray> sample(vec const &, vec const &) const;
  double pdf(vec const &, vec const &, vec const &d) const;
  sampled_spectrum get_sp(vec const &, vec const &d, sampled_wl const &wl) const;
};

struct uniform: base {
  Spectrum::ptr sp;

  uniform(Spectrum::ptr s): base(true, true), sp(s) {}

  biased<ray> sample(vec const &, vec const &n) const;
  double pdf(vec const &, vec const &n, vec const &d) const;
  sampled_spectrum get_sp(vec const &, vec const &, sampled_wl const &wl) const;
};

struct environment: base {
  Image::base const *img;
  Sampler::discrete2D samp;
  double strength, area;
  Matrix::mat rot;
  environment(Image::base const *i, double s, double a = 0.);
  biased<ray> sample(vec const &, vec const &n) const;
  double pdf(vec const &, vec const &n, vec const &d) const;
  sampled_spectrum get_sp(vec const &, vec const &dir, sampled_wl const &wl) const;
};

struct spherical: base {
  Spectrum::ptr sp;
  Solid::sphere *sph;
  double strength;

  spherical(Spectrum::ptr s, vec const &c, double r);
  biased<ray> sample(vec const &pos, vec const &) const;
  double pdf(vec const &pos, vec const &, vec const &d) const;
  sampled_spectrum get_sp(vec const &pos, vec const &dir, sampled_wl const &wl) const;
};

}

#endif
