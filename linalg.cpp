#include "linalg.hpp"

namespace Vector {

point2 from_sphere(vec const &dir) {
  vec d = (1. / (std::abs(dir[0]) + std::abs(dir[1]) + std::abs(dir[2]))) * dir;
  double u = d[0], v = -d[2];
  if (d[1] < 0.) {
    double t = u;
    u = (1. - std::abs(v)) * std::copysign(1., u);
    v = (1. - std::abs(t)) * std::copysign(1., v);
  }
  return { (u + 1.) * 0.5, (v + 1.) * 0.5 };
}

vec to_sphere(point2 const &uv) {
  double x = uv[0] * 2. - 1., z = 1. - uv[1] * 2.;
  double y = 1. - (std::abs(x) + std::abs(z));
  if (y < 0) {
    double t = x;
    x = (1. - std::abs(z)) * std::copysign(1., x);
    z = (1. - std::abs(t)) * std::copysign(1., z);
  }
  return normalize(vec { x, y, z });
}

}

namespace Matrix {

vec operator*(mat const &m, vec const &v) {
  return {
    m[0] * v[0] + m[1] * v[1] + m[2] * v[2],
    m[3] * v[0] + m[4] * v[1] + m[5] * v[2],
    m[6] * v[0] + m[7] * v[1] + m[8] * v[2] };
}

mat operator+(mat const &m, mat const &n) {
  mat r;
  for (int i = 0; i < 9; ++i) { r[i] = m[i] + n[i]; }
  return r;
}

mat operator*(double a, mat const &m) {
  mat r;
  for (int i = 0; i < 9; ++i) { r[i] = a * m[i]; }
  return r;
}

mat diag(double a, double b, double c) {
  return { a, 0., 0., 0., b, 0., 0., 0., c };
}

mat transpose(mat const &m) {
  return {
    m[0], m[3], m[6],
    m[1], m[4], m[7],
    m[2], m[5], m[8] };
}

mat rotation(vec const &v, double a) {
  double c = cos(a), s = sin(a), cc = 1 - c;
  mat res = { 0., -v[2], v[1], v[2], 0., -v[0], -v[1], v[0], 0. };
  res = diag(c, c, c) + s * res;
  for (int i = 0; i < 9; ++i) {
    res[i] += cc * v[i / 3] * v[i % 3];
  }
  return res;
}

mat rotation_zy(vec const &z, vec const &y) {
  vec p = normalize(z);
  vec q = normalize(cross(y, p));
  vec r = normalize(cross(p, q));
  return { q[0], r[0], p[0], q[1], r[1], p[1], q[2], r[2], p[2] };
}

}
