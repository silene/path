#ifndef GEOMETRY_HPP
#define GEOMETRY_HPP

#include <cmath>

#include "debug.hpp"
#include "linalg.hpp"

namespace Box {

struct box {
  Vector::small p1, p2;
};

box constexpr empty { { INFINITY, INFINITY, INFINITY }, { -INFINITY, -INFINITY, -INFINITY } };

box constexpr whole { { -INFINITY, -INFINITY, -INFINITY }, { INFINITY, INFINITY, INFINITY } };

using Vector::vec;

struct checker {
  Vector::small pos, invdir;
  checker(vec const &p, vec const &d)
    : pos(p)
  { for (int i = 0; i < 3; ++i) { invdir[i] = 1.f / d[i]; } }
  bool operator()(box const &b, double dmax) const;
};

box merge(vec const &v1, vec const &v2);
box merge(box const &b1, vec const &v);
box merge(box const &b1, box const &b2);
box intersect(box const &b1, box const &b2);
vec center(box const &b);

}

namespace Ball {

using Vector::vec;

struct ball {
  vec center;
  double radius;
};

ball constexpr empty { vec(), -INFINITY };

void inflate(ball &b, vec const &v);

}

namespace Transform {

using Vector::vec;
using Matrix::mat;
using Box::box;

struct iso {
  vec center;
  double scale, iscale;
  mat rotate;
  iso(vec const &c, double sc, vec const &a, double r);

  vec to_relative(vec const &pos) const;
  vec of_relative(vec const &pos) const;
  box bounds(box const &b) const;

  vec to_relative_d(vec const &dir) const {
    return rotate * dir;
  }

  vec of_relative_d(vec const &dir) const {
    return transpose(rotate) * dir;
  }
};

using ptr = iso const *;

}

point2 toUV(Vector::vec const &p, Vector::vec const &q, Vector::vec const &r);

struct object;

namespace Geometry {

namespace {
using Vector::vec;
}

struct contact {
  vec pos, normal; // in the local basis
  point2 uv;
  int data;
};

struct full_contact {
  contact co;
  double dist = INFINITY;
  object const *obj = NULL;
  int data;
};

struct contact_finder {
  vec pos, dir;
  Box::checker checker;
  contact_finder(vec const &p, vec const &d);
  bool check_range(int ib, int ie, full_contact &) const;
  bool check_split(int sp, full_contact &) const;
  full_contact operator()(double dmax) const;
};

void prepare_bounds();

}

#endif
