#include <cassert>

#include "solid.hpp"

namespace Solid {

double sphere::distance(vec const &pos, vec const &dir) const {
  vec p = pos - center;
  double b = dir | p;
  double c = (p | p) - radius * radius;
  if (b >= 0 && c >= 0) return INFINITY;
  double d = b * b - c;
  if (d < 0) return INFINITY;
  d = sqrt(d);
  double t1 = -b + (b > 0 ? -d : d);
  double t2 = c / t1;
  if (t1 < 0 || (t2 >= 0 && t2 < t1)) t1 = t2;
  return t1 < 0 ? INFINITY : t1;
}

bool sphere::complete(contact &co, int) const {
  co.normal = normalize(co.pos - center);
  return true;
}

Box::box sphere::bounds(int, Transform::ptr t) const {
  assert(t == NULL);
  vec r { radius, radius, radius };
  return { center - r, center + r };
}

Ball::ball sphere::sbounds(int, Transform::ptr t) const {
  assert(t == NULL);
  return { center, radius };
}

bool sphere::inside(vec const &pos) const {
  vec p = pos - center;
  return (p | p) - radius * radius <= 1e-10;
}

double cylinder::distance(vec const &pos, vec const &dir, contact &, int) const {
  vec p = pos - basis;
  double ia2 = 1 / (axis | axis);
  vec q = p - (p | axis) * ia2 * axis;
  vec r = dir - (dir | axis) * ia2 * axis;
  double a = (r | r);
  double b = (q | r);
  double c = (q | q) - radius * radius;
  if (b >= 0 && c >= 0) return INFINITY;
  double d = b * b - a * c;
  if (d < 0) return INFINITY;
  d = sqrt(d);
  double t1 = (-b + (b > 0 ? -d : d)) / a;
  double t2 = a * c / t1;
  if (t1 < 0 || (t2 >= 0 && t2 < t1)) t1 = t2;
  if (t1 < 0) return INFINITY;
  double f = ((p + t1 * dir) | axis) * ia2;
  if (f >= 0 && f <= 1) return t1;
  if (t2 < 0) return INFINITY;
  f = ((p + t2 * dir) | axis) * ia2;
  if (f >= 0 && f <= 1) return t2;
  return INFINITY;
}

bool cylinder::complete(contact &co, int) const {
  vec p = co.pos - basis;
  double f = (p | axis) / (axis | axis);
  co.normal = normalize(p - f * axis);
  return true;
}

Box::box cylinder::bounds(int, Transform::ptr t) const {
  assert(t == NULL);
  vec r { radius, radius, radius };
  Box::box b = Box::merge(basis - r, basis + r);
  b = Box::merge(b, basis + axis - r);
  b = Box::merge(b, basis + axis + r);
  return b;
}

Ball::ball cylinder::sbounds(int, Transform::ptr t) const {
  assert(false);
}

bool cylinder::inside(vec const &pos) const {
  vec p = pos - basis;
  double f = (p | axis) / (axis | axis);
  if (f < 0 || f > 1) return false;
  p = p - f * axis;
  return (p | p) - radius * radius <= 1e-10;
}

double plane::distance(vec const &pos, vec const &dir, contact &, int) const {
  double d = (pos | normal_) + dist;
  double s = (dir | normal_);
  if (s == 0) return d == 0 ? 0 : INFINITY;
  double r = - d / s;
  return r < 0 ? INFINITY : r;
}

Box::box plane::bounds(int, Transform::ptr t) const {
  assert(t == NULL);
  Box::box b = Box::empty;
  double d = Settings::max_depth;
  for (int i = 0; i < 3; ++i) {
    if (normal_[i] == 0) continue;
    for (int j = 0; j < 4; ++j) {
      vec v { d, d, d };
      v[i] = 0.;
      if (j & 1) { v[(i + 1) % 3] = -d; }
      if (j & 2) { v[(i + 2) % 3] = -d; }
      v[i] = - ((v | normal_) + dist) / normal_[i];
      if (v[i] < -d) continue;
      if (v[i] > d) continue;
      b = Box::merge(b, v);
    }
  }
  return b;
}

Ball::ball plane::sbounds(int, Transform::ptr t) const {
  assert(false);
}

bool plane::inside(vec const &pos) const {
  return (pos | normal_) + dist <= 1e-10;
}

double union_::distance(vec const &pos, vec const &dir, contact &co, int p) const {
  if (p < parts1) return obj1->distance(pos, dir, co, p);
  return obj2->distance(pos, dir, co, p - parts1);
}

bool union_::complete(contact &co, int p) const {
  if (p < parts1)
    return obj1->complete(co, p) && !obj2->inside(co.pos);
  return obj2->complete(co, p - parts1) && !obj1->inside(co.pos);
}

Box::box union_::bounds(int p, Transform::ptr t) const {
  assert(t == NULL);
  if (p < parts1) return obj1->bounds(p, t);
  return obj2->bounds(p - parts1, t);
}

Ball::ball union_::sbounds(int, Transform::ptr) const {
  assert(false);
}

intersection::intersection(ptr o1, ptr o2)
  : obj1(o1), obj2(o2), parts1(obj1->subparts()),
    parts(parts1 + obj2->subparts())
{
  Box::box b = Box::empty;
  for (int n = 0; n < parts1; ++n) {
    b = Box::merge(b, obj1->bounds(n, NULL));
  }
  bnds = b;
  b = Box::empty;
  for (int n = parts1; n < parts; ++n) {
    b = Box::merge(b, obj2->bounds(n - parts1, NULL));
  }
  bnds = Box::intersect(bnds, b);
}

double intersection::distance(vec const &pos, vec const &dir, contact &co, int p) const {
  if (p < parts1) return obj1->distance(pos, dir, co, p);
  return obj2->distance(pos, dir, co, p - parts1);
}

bool intersection::complete(contact &co, int p) const {
  if (p < parts1)
    return obj1->complete(co, p) && obj2->inside(co.pos);
  return obj2->complete(co, p - parts1) && obj1->inside(co.pos);
}

Box::box intersection::bounds(int p, Transform::ptr t) const {
  assert(t == NULL);
  if (p < parts1) return Box::intersect(bnds, obj1->bounds(p, t));
  return Box::intersect(bnds, obj2->bounds(p - parts1, t));
}

Ball::ball intersection::sbounds(int, Transform::ptr) const {
  assert(false);
}

difference::difference(ptr o1, ptr o2)
  : obj1(o1), obj2(o2), parts1(obj1->subparts()),
    parts(parts1 + obj2->subparts())
{
  Box::box b = Box::empty;
  for (int n = 0; n < parts1; ++n) {
    b = Box::merge(b, obj1->bounds(n, NULL));
  }
  bnds1 = b;
}

double difference::distance(vec const &pos, vec const &dir, contact &co, int p) const {
  if (p < parts1) return obj1->distance(pos, dir, co, p);
  return obj2->distance(pos, dir, co, p - parts1);
}

bool difference::complete(contact &co, int p) const {
  if (p < parts1)
    return obj1->complete(co, p) && !obj2->inside(co.pos);
  if (!obj2->complete(co, p - parts1)) return false;
  co.normal = -co.normal;
  return obj1->inside(co.pos);
}

Box::box difference::bounds(int p, Transform::ptr t) const {
  assert(t == NULL);
  if (p < parts1) return obj1->bounds(p, t);
  return Box::intersect(bnds1, obj2->bounds(p - parts1, t));
}

Ball::ball difference::sbounds(int, Transform::ptr) const {
  assert(false);
}

}
