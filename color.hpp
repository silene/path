#ifndef COLOR_HPP
#define COLOR_HPP

#include <array>
#include <vector>

#include "base.hpp"
#include "linalg.hpp"
#include "image.hpp"

namespace Color {

extern Matrix::mat RGBtoXYZ;
extern Vector::vec RGBtoY;
extern Matrix::mat XYZtoRGB;

Vector::vec toXYZ(double wl);
double toY(double wl);

}

namespace Spectrum {

inline sampled_spectrum &operator+=(sampled_spectrum &u, sampled_spectrum const &v) {
  for (int i = 0; i < nb_s; ++i) { u[i] += v[i]; }
  return u;
}

inline sampled_spectrum operator+(sampled_spectrum const &u, sampled_spectrum const &v) {
  sampled_spectrum w(u);
  return w += v;
}

inline sampled_spectrum &operator*=(sampled_spectrum &u, float a) {
  for (int i = 0; i < nb_s; ++i) { u[i] *= a; }
  return u;
}

inline sampled_spectrum operator*(float a, sampled_spectrum const &u) {
  sampled_spectrum w;
  for (int i = 0; i < nb_s; ++i) { w[i] = a * u[i]; }
  return w;
}

inline sampled_spectrum operator*(sampled_spectrum const &u, sampled_spectrum const &v) {
  sampled_spectrum w;
  for (int i = 0; i < nb_s; ++i) { w[i] = u[i] * v[i]; }
  return w;
}

struct full_spectrum: std::array<float, nb_wl> {
  full_spectrum();
};

void add(full_spectrum &fs, sampled_wl const &wl, sampled_spectrum const &s);

Vector::vec toXYZ(full_spectrum const &s);

struct blackbody: base {
  double strength;
  int temp;
  //full_spectrum fsp;

  blackbody(double f, int t);
  double aux(double wl) const;
  sampled_spectrum sample(point2 const &, sampled_wl const &l) const;
};

struct uniform: base {
  double strength;
  uniform(double s): strength(s) {}
  sampled_spectrum sample(point2 const &, sampled_wl const &) const {
    return sampled_spectrum(strength);
  }
};

struct sampled: base {
  typedef std::vector<std::pair<double,double>> samples_t;
  samples_t samples;
  sampled(samples_t const &s);
  sampled_spectrum sample(point2 const &, sampled_wl const &wl) const;
};

struct from_wl: base {
  double strength, fact;
  double lambda, width;
  from_wl(double s, double l, double w)
    : strength(s), lambda(l), width(w) {
    fact = 1.; // / (width * sqrt(2 * M_PI));
  }

  sampled_spectrum sample(point2 const &, sampled_wl const &wl) const;
};

struct xyY: base {
  double strength;
  double x, y;
  xyY(double x_, double y_, double s)
    : strength(s), x(x_), y(y_) {}
  sampled_spectrum sample(point2 const &, sampled_wl const &wl) const;
};

sampled_spectrum fromXYZ(Vector::vec const &c, sampled_wl const &wl);

struct from_texture: base {
  Image::ptr img;
  from_texture(Image::ptr i): img(i) {}
  sampled_spectrum sample(point2 const &uv, sampled_wl const &wl) const;
};

}

#endif
