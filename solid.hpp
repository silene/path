#ifndef SOLID_HPP
#define SOLID_HPP

#include "base.hpp"
#include "geometry.hpp"
#include "linalg.hpp"

namespace Solid {

namespace {
using Vector::vec;
using Geometry::contact;
}

struct sphere: base {
  vec center;
  double radius;
  sphere(vec const &c, double r)
    : center(c), radius(r) {}

  double distance(vec const &pos, vec const &dir) const;

  double distance(vec const &pos, vec const &dir, contact &, int) const {
    return distance(pos, dir);
  }

  bool complete(contact &co, int) const;
  Box::box bounds(int, Transform::ptr t) const;
  Ball::ball sbounds(int, Transform::ptr t) const;
  bool inside(vec const &pos) const;
};

struct cylinder: base {
  vec basis, axis;
  double radius;
  cylinder(vec const &b, vec const &a, double r)
    : basis(b), axis(a), radius(r) {}

  double distance(vec const &pos, vec const &dir, contact &, int) const;
  bool complete(contact &co, int) const;
  Box::box bounds(int, Transform::ptr t) const;
  Ball::ball sbounds(int, Transform::ptr t) const;
  bool inside(vec const &pos) const;
};

struct plane: base {
  vec normal_;
  double dist;
  plane(vec const &n, double d): normal_(n), dist(d) {}

  double distance(vec const &pos, vec const &dir, contact &, int) const;

  bool complete(contact &co, int) const {
    co.normal = normal_;
    return true;
  }

  Box::box bounds(int, Transform::ptr t) const;
  Ball::ball sbounds(int, Transform::ptr t) const;
  bool inside(vec const &pos) const;
};

struct union_: base {
  ptr obj1, obj2;
  int parts1, parts;

  union_(ptr o1, ptr o2)
    : obj1(o1), obj2(o2), parts1(obj1->subparts()),
      parts(parts1 + obj2->subparts()) {}

  double distance(vec const &pos, vec const &dir, contact &co, int p) const;
  bool complete(contact &co, int p) const;
  int subparts() const { return parts; }
  Box::box bounds(int p, Transform::ptr t) const;
  Ball::ball sbounds(int, Transform::ptr) const;

  bool inside(vec const &pos) const {
    return obj1->inside(pos) || obj2->inside(pos);
  }
};

struct intersection: base {
  ptr obj1, obj2;
  int parts1, parts;
  Box::box bnds;

  intersection(ptr o1, ptr o2);

  double distance(vec const &pos, vec const &dir, contact &co, int p) const;
  bool complete(contact &co, int p) const;
  int subparts() const { return parts; }
  Box::box bounds(int p, Transform::ptr t) const;
  Ball::ball sbounds(int, Transform::ptr) const;

  bool inside(vec const &pos) const {
    return obj1->inside(pos) && obj2->inside(pos);
  }
};

struct difference: base {
  ptr obj1, obj2;
  int parts1, parts;
  Box::box bnds1;

  difference(ptr o1, ptr o2);
  double distance(vec const &pos, vec const &dir, contact &co, int p) const;
  bool complete(contact &co, int p) const;
  int subparts() const { return parts; }
  Box::box bounds(int p, Transform::ptr t) const;
  Ball::ball sbounds(int, Transform::ptr) const;

  bool inside(vec const &pos) const {
    return obj1->inside(pos) && !obj2->inside(pos);
  }
};

}

#endif
