#include <algorithm>
#include <cassert>

#include "sdf.hpp"

namespace SDF {

vec base::normal(vec const &pos, double d) const {
  d = std::max(d * 1e-6, 1e-6);
  double md = -d;
  double
    d0 = distance(pos + vec { d, d, d }),
    d1 = distance(pos + vec { md, md, d }),
    d2 = distance(pos + vec { md, d, md }),
    d3 = distance(pos + vec { d, md, md });
  vec v = {
    d0 - d1 - d2 + d3,
    d0 - d1 + d2 - d3,
    d0 + d1 - d2 - d3 };
  return normalize(v);
}

union_::union_(std::initializer_list<ptr> s): sdfs(s) {}

double union_::distance(vec const &pos) const {
  double d = INFINITY;
  for (ptr s: sdfs) { d = std::min(d, s->distance(pos)); }
  return d;
}

vec union_::normal(vec const &pos, double d_) const {
  double bd = INFINITY;
  ptr bs = NULL;
  for (ptr s : sdfs) {
    double d = s->distance(pos);
    if (d < bd) { bd = d; bs = s; }
  }
  return bs->normal(pos, d_);
}

Box::box union_::bounds() const {
  Box::box b = Box::empty;
  for (ptr s: sdfs) { b = Box::merge(b, s->bounds()); }
  return b;
}

double smooth_union::distance(vec const &pos) const {
  double d1 = s1->distance(pos);
  double d2 = s2->distance(pos);
  double d = std::abs(d1 - d2);
  if (d > k) return std::min(d1, d2);
  double h = 1 - d * ik;
  return std::min(d1, d2) - h * h * k * 0.25;
}

vec smooth_union::normal(vec const &pos, double d_) const {
  double d1 = s1->distance(pos);
  double d2 = s2->distance(pos);
  double d = std::abs(d1 - d2);
  if (d > k)
    return d1 < d2 ? s1->normal(pos, d_) : s2->normal(pos, d_);
  double h = (1 - d * ik) * 0.5;
  vec n1 = s1->normal(pos, d_);
  vec n2 = s2->normal(pos, d_);
  vec n = d1 < d2 ? mix(n1, n2, h) : mix(n2, n1, h);
  return normalize(n);
}

Box::box smooth_union::bounds() const {
  return Box::merge(s1->bounds(), s2->bounds());
}

double intersection::distance(vec const &pos) const {
  double d1 = s1->distance(pos);
  double d2 = s2->distance(pos);
  return std::max(d1, d2);
}

vec intersection::normal(vec const &pos, double d_) const {
  double d1 = s1->distance(pos);
  double d2 = s2->distance(pos);
  return d1 > d2 ? s1->normal(pos, d_) : s2->normal(pos, d_);
}

Box::box intersection::bounds() const {
  return Box::intersect(s1->bounds(), s2->bounds());
}

double difference::distance(vec const &pos) const {
  double d1 = s1->distance(pos);
  double d2 = s2->distance(pos);
  return std::max(d1, -d2);
}

vec difference::normal(vec const &pos, double d_) const {
  double d1 = s1->distance(pos);
  double d2 = s2->distance(pos);
  return d1 > -d2 ? s1->normal(pos, d_) : - s2->normal(pos, d_);
}

double sphere::distance(vec const &pos) const {
  return norm(pos - center) - radius;
}

vec sphere::normal(vec const &pos, double) const {
  return normalize(pos - center);
}

Box::box sphere::bounds() const {
  vec r { radius, radius, radius };
  return { center - r, center + r };
}

double cylinder::distance(vec const &pos) const {
  vec v = pos - center;
  return norm(v - (v | dir) * dir) - radius;
}

vec cylinder::normal(vec const &pos, double) const {
  vec v = pos - center;
  return normalize(v - (v | dir) * dir);
}

double plane::distance(vec const &pos) const {
  return (pos | normal_) + dist;
}

double box::distance(vec const &pos) const {
  vec p = pos - center;
  for (int i = 0; i < 3; ++i) { p[i] = std::abs(p[i]); }
  p = p - dim;
  double d = std::max({ p[0], p[1], p[2] });
  if (d <= 0) return d;
  vec q;
  for (int i = 0; i < 3; ++i) { q[i] = std::max(p[i], 0.); }
  return norm(q);
}

Box::box box::bounds() const {
  return { center - dim, center + dim };
}

Box::box inflated::bounds() const {
  Box::box b = s->bounds();
  vec p1 { b.p1[0] - radius, b.p1[1] - radius, b.p1[2] - radius };
  vec p2 { b.p2[0] + radius, b.p2[1] + radius, b.p2[2] + radius };
  return { p1, p2 };
}

}

namespace Solid {

double sdf::distance(vec const &pos, vec const &dir, contact &, int) const {
  double prev = INFINITY;
  bool optimistic = false;
  for (double l = 0; l < Settings::max_depth; ) {
    vec npos = pos + l * dir;
    double d = s->distance(npos);
    if (optimistic && d <= prev) {
      // The bet was wrong, cancel the extra jump.
      optimistic = false;
      l -= prev;
      continue;
    }
    if (d < 1e-6) return l;
    // Optimistically jump twice as far, if the distance is growing.
    optimistic = (prev < d);
    prev = d;
    l += d;
    if (!optimistic) continue;
    l += d;
  }
  return INFINITY;
}

Ball::ball sdf::sbounds(int, Transform::ptr t) const {
  assert(false);
}

}
