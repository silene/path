#include "geometry.hpp"

namespace Box {

bool checker::operator()(box const &b, double dmax) const {
  ++dbg->boxes;
  float v1 = 0, v2 = dmax;
  for (int i = 0; i < 3; ++i) {
    float t1 = (b.p1[i] - pos[i]) * invdir[i];
    float t2 = (b.p2[i] - pos[i]) * invdir[i];
    if (t2 < t1) std::swap(t1, t2);
    if (t1 > v1) v1 = t1;
    if (t2 < v2) v2 = t2;
    if (v2 < v1) return false;
  }
  return true;
}

box merge(vec const &v1, vec const &v2) {
  box b;
  for (int i = 0; i < 3; ++i) { b.p1[i] = std::min<float>(v1[i], v2[i]); }
  for (int i = 0; i < 3; ++i) { b.p2[i] = std::max<float>(v1[i], v2[i]); }
  return b;
}

box merge(box const &b1, vec const &v) {
  box b;
  for (int i = 0; i < 3; ++i) { b.p1[i] = std::min<float>(b1.p1[i], v[i]); }
  for (int i = 0; i < 3; ++i) { b.p2[i] = std::max<float>(b1.p2[i], v[i]); }
  return b;
}

box merge(box const &b1, box const &b2) {
  box b;
  for (int i = 0; i < 3; ++i) { b.p1[i] = std::min<float>(b1.p1[i], b2.p1[i]); }
  for (int i = 0; i < 3; ++i) { b.p2[i] = std::max<float>(b1.p2[i], b2.p2[i]); }
  return b;
}

box intersect(box const &b1, box const &b2) {
  box b;
  for (int i = 0; i < 3; ++i) { b.p1[i] = std::max<float>(b1.p1[i], b2.p1[i]); }
  for (int i = 0; i < 3; ++i) { b.p2[i] = std::min<float>(b1.p2[i], b2.p2[i]); }
  return b;
}

vec center(box const &b) {
  vec v;
  for (int i = 0; i < 3; ++i) { v[i] = 0.5 * (b.p1[i] + b.p2[i]); }
  return v;
}

}


namespace Ball {

void inflate(ball &b, vec const &v) {
  if (b.radius < 0.) {
    b.center = v;
    b.radius = 0.;
    return;
  }
  double d = norm(v - b.center);
  if (d <= b.radius) return;
  b.center = mix(b.center, v, 0.5 - 0.5 * b.radius / d);
  b.radius = 0.5 * (b.radius + d);
}

}

namespace Transform {

iso::iso(vec const &c, double sc, vec const &a, double r):
  center(c), scale(sc), iscale(1 / sc), rotate(Matrix::rotation(a, r))
{}

vec iso::to_relative(vec const &pos) const {
  return iscale * (rotate * (pos - center));
}

vec iso::of_relative(vec const &pos) const {
  return scale * (transpose(rotate) * pos) + center;
}

box iso::bounds(box const &b) const {
  box r = Box::empty;
  for (int i = 0; i < 8; ++i) {
    vec v;
    for (int j = 0; j < 3; ++j) {
      v[j] = (i & (1 << j)) ? b.p2[j] : b.p1[j];
    }
    r = Box::merge(r, of_relative(v));
  }
  return r;
}

}

point2 toUV(Vector::vec const &p, Vector::vec const &q, Vector::vec const &r) {
  // the barycentric coordinates u, v, 1-u-v are proportional
  // to the areas of the three subtriangles opr, oqr, pqr
  Vector::vec n = cross(p, q);
  double d = (n | n);
  double u = (n | cross(r, q)) / d;
  double v = (n | cross(p, r)) / d;
  return { u, v };
}
