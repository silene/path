#ifndef SDF_HPP
#define SDF_HPP

#include <vector>

#include "geometry.hpp"
#include "linalg.hpp"
#include "solid.hpp"

namespace SDF {

using Vector::vec;
using Matrix::mat;

struct base {
  virtual double distance(vec const &) const = 0;
  virtual Box::box bounds() const = 0;
  virtual vec normal(vec const &pos, double d) const;
  virtual ~base() = default;
};

using ptr = base const *;

struct union_: base {
  std::vector<ptr> sdfs;
  union_(std::initializer_list<ptr> s);
  double distance(vec const &pos) const;
  vec normal(vec const &pos, double d_) const;
  Box::box bounds() const;
};

struct smooth_union: base {
  double k, ik;
  ptr s1, s2;
  smooth_union(double k_, ptr s1_, ptr s2_)
    : k(k_), ik(1 / k_), s1(s1_), s2(s2_) {}

  double distance(vec const &pos) const;
  vec normal(vec const &pos, double d_) const;
  Box::box bounds() const;
};

struct intersection: base {
  ptr s1, s2;
  intersection(ptr s1_, ptr s2_): s1(s1_), s2(s2_) {}

  double distance(vec const &pos) const;
  vec normal(vec const &pos, double d_) const;
  Box::box bounds() const;
};

struct difference: base {
  ptr s1, s2;
  difference(ptr s1_, ptr s2_): s1(s1_), s2(s2_) {}

  double distance(vec const &pos) const;
  vec normal(vec const &pos, double d_) const;

  Box::box bounds() const {
    return s1->bounds();
  }
};

struct sphere: base {
  vec center;
  double radius;
  sphere(vec const &c, double r): center(c), radius(r) {}

  double distance(vec const &pos) const;
  vec normal(vec const &pos, double) const;
  Box::box bounds() const;
};

struct cylinder: base {
  vec center, dir;
  double radius;
  cylinder(vec const &c, vec const &d, double r)
    : center(c), dir(d), radius(r) {}

  double distance(vec const &pos) const;
  vec normal(vec const &pos, double) const;

  Box::box bounds() const {
    return Box::whole;
  }
};

struct plane: base {
  vec normal_;
  double dist;
  plane(vec const &n, double d): normal_(n), dist(d) {}

  double distance(vec const &pos) const;

  vec normal(vec const &, double) const {
    return normal_;
  }

  Box::box bounds() const {
    return Box::whole;
  }
};

struct box: base {
  vec center, dim;
  box(vec const &c, vec const &d): center(c), dim(d) {}

  double distance(vec const &pos) const;
  Box::box bounds() const;
};

#if 0
struct iso: base {
  ptr s;
  vec center;
  mat rotate;

  iso(ptr s_, vec const &c, vec const &a, double r):
    s(s_), center(c), rotate(Matrix::rotation(a, r)) {}

  double distance(vec const &pos) const {
    vec p = rotate * (pos - center);
    return s->distance(p);
  }

  vec normal(vec const &pos, double d_) const {
    vec p = rotate * (pos - center);
    return transpose(rotate) * s->normal(p, d_);
  }
};
#endif

struct inflated: base {
  ptr s;
  double radius;
  inflated(double r, ptr s_): s(s_), radius(r) {}

  double distance(vec const &pos) const {
    return s->distance(pos) - radius;
  }

  vec normal(vec const &pos, double d_) const {
    return s->normal(pos, d_);
  }

  Box::box bounds() const;
};

}

namespace Solid {

struct sdf: base {
  SDF::ptr s;
  sdf(SDF::ptr s_): s(s_) {}
  double distance(vec const &pos, vec const &dir, contact &, int) const;

  bool complete(contact &co, int) const {
    co.normal = s->normal(co.pos, 0.);
    return true;
  }

  Box::box bounds(int, Transform::ptr t) const {
    Box::box b = s->bounds();
    return t ? t->bounds(b) : b;
  }

  Ball::ball sbounds(int, Transform::ptr t) const;

  bool inside(vec const &pos) const {
    return s->distance(pos) <= 1e-5;
  }
};

}

#endif
