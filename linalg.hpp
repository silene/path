#ifndef LINALG_HPP
#define LINALG_HPP

#include <array>
#include <cmath>

using point2 = std::array<double, 2>;

namespace Vector {

struct small: std::array<float, 3> {};

struct vec: std::array<double, 3> {
  operator small() const {
    return small { (float)(*this)[0], (float)(*this)[1], (float)(*this)[2] };
  }
};

inline vec &operator+=(vec &u, vec const &v) {
  for (int i = 0; i < 3; ++i) { u[i] += v[i]; }
  return u;
}

inline vec operator+(vec const &u, vec const &v) {
  vec w(u);
  return w += v;
}

inline vec operator-(vec const &u, vec const &v) {
  vec w;
  for (int i = 0; i < 3; ++i) { w[i] = u[i] - v[i]; }
  return w;
}

inline vec operator-(vec const &u) {
  vec w;
  for (int i = 0; i < 3; ++i) { w[i] = -u[i]; }
  return w;
}

inline vec &operator*=(vec &u, double a) {
  for (int i = 0; i < 3; ++i) { u[i] *= a; }
  return u;
}

inline vec operator*(double a, vec const &u) {
  vec w;
  for (int i = 0; i < 3; ++i) { w[i] = a * u[i]; }
  return w;
}

inline double operator|(vec const &u, vec const &v) {
  double r = 0;
  for (int i = 0; i < 3; ++i) { r += u[i] * v[i]; }
  return r;
}

inline double norm(vec const &u) {
  return sqrt(u|u);
}

inline vec normalize(vec const &u) {
  return (1 / norm(u)) * u;
}

inline vec mix(vec const &u, vec const &v, double k) {
  return (1 - k) * u + k * v;
}

inline vec operator*(vec const &u, vec const &v) {
  vec w;
  for (int i = 0; i < 3; ++i) { w[i] = u[i] * v[i]; }
  return w;
}

inline vec cross(vec const &u, vec const &v) {
  vec w;
  w[0] = u[1] * v[2] - u[2] * v[1];
  w[1] = u[2] * v[0] - u[0] * v[2];
  w[2] = u[0] * v[1] - u[1] * v[0];
  return w;
}

point2 from_sphere(vec const &dir);
vec to_sphere(point2 const &uv);

}


namespace Matrix {

using Vector::vec;

struct mat: std::array<double, 9> {};

vec operator*(mat const &m, vec const &v);
mat operator+(mat const &m, mat const &n);
mat operator*(double a, mat const &m);

mat diag(double a, double b, double c);
mat transpose(mat const &m);
mat rotation(vec const &v, double a);
mat rotation_zy(vec const &z, vec const &y);

}

#endif
