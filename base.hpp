#ifndef BASE_HPP
#define BASE_HPP

#include <array>

#include "geometry.hpp"
#include "linalg.hpp"
#include "sampler.hpp"

namespace Settings {

enum skind { None, Plain, Weighted };

extern double max_depth;
extern skind shadows;
extern int min_steps, max_steps;
extern int min_samples, max_samples;
extern double variance;
extern bool regularize;

}

namespace Spectrum {

int constexpr min_wl = 360, max_wl = 750;
int constexpr nb_wl = 390;
int constexpr nb_s = 4;

struct sampled_spectrum: std::array<float, nb_s> {
  sampled_spectrum() = default;
  sampled_spectrum(double d) { fill(d); }
  bool zero() const {
    for (int i = 0; i < nb_wl; ++i) {
      if ((*this)[i]) return false;
    }
    return true;
  }
};

struct sampled_wl {
  std::array<double, nb_s> lambda, pdf;
  sampled_wl(double s);
};

struct base {
  virtual sampled_spectrum sample(point2 const &, sampled_wl const &) const = 0;
  virtual ~base() = default;
};

using ptr = base const *;

}

namespace Solid {

namespace {
using Vector::vec;
}

struct contact {
  vec pos, normal; // in the local basis
  point2 uv;
  int data;
};

struct base {
  virtual double distance(vec const &, vec const &, contact &, int) const = 0;
  virtual bool complete(contact &, int) const = 0;
  virtual vec snormal(contact const &co) const { return co.normal; }
  virtual point2 uv(contact const &) const { return { 0., 0. }; }
  virtual int subparts() const { return 1; }
  virtual Box::box bounds(int, Transform::ptr) const = 0;
  virtual Ball::ball sbounds(int, Transform::ptr) const = 0;
  virtual bool inside(vec const &) const = 0;
  virtual ~base() = default;
};

using ptr = base const *;

}

namespace Material {

namespace {
using Vector::vec;
using Spectrum::sampled_spectrum;
using Spectrum::sampled_wl;
}

struct interaction {
  vec pos, normal, out;
  point2 uv;
};

enum skind { Specular, SemiSpecular, Diffuse };

struct ray {
  sampled_spectrum sp;
  vec dir;
  skind specular;
};

enum mkind { Solid, Emissive, Transmitive };

struct base {
  mkind kind;
  bool has_bxdf, has_pdf;
  base(mkind k, bool b, bool p)
    : kind(k), has_bxdf(b), has_pdf(p) {}
  virtual sampled_spectrum bxdf(interaction const &, vec const &inc, sampled_wl const &) const;
  virtual double pdf(interaction const &, vec const &inc) const;
  virtual biased<ray> sample(interaction const &, sampled_wl const &) const = 0;
  virtual base const *regularize() const { return this; }
  virtual ~base() = default;
};

using ptr = base const *;

}

namespace Light {

namespace {
using Vector::vec;
using Spectrum::sampled_spectrum;
using Spectrum::sampled_wl;
}

struct ray {
  vec dir;
  double dist;
};

using Vector::vec;

struct base {
  bool has_pdf, surrounding;
  base(bool b, bool s): has_pdf(b), surrounding(s) {}
  virtual biased<ray> sample(vec const &pos, vec const &n) const = 0;
  virtual double pdf(vec const &pos, vec const &n, vec const &dir) const;
  virtual sampled_spectrum get_sp(vec const &pos, vec const &dir, sampled_wl const &wl) const = 0;
  virtual ~base() = default;
};

using ptr = base const *;

}

struct object {
  Solid::ptr solid;
  Material::ptr material;
  Transform::ptr transf;
};

namespace Scene {

extern std::vector<Light::ptr> lights;
extern std::vector<object> objects;

}

namespace Camera {

namespace {
  using Vector::vec;
}

struct base {
  virtual std::pair<vec, vec> get(double fx, double fy) const = 0;
  virtual ~base() = default;
};

using ptr = base const *;

extern ptr camera;

}

#endif
