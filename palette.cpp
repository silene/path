#include <algorithm>
#include <atomic>
#include <cassert>
#include <fstream>
#include <functional>
#include <iostream>
#include <random>
#include <sstream>
#include <string>
#include <thread>
#include <vector>

struct image {
  int width, height;
  std::string data;
  image(int w, int h): width(w), height(h), data(w * h * 3, '\0') {}
  void write(int x, int y, int r, int g, int b) {
    if (!(0 <= x && x < width && 0 <= y && y < height)) {
      std::cerr << x << ' ' << y << '\n';
    }
    assert(0 <= x && x < width && 0 <= y && y < height);
    int i = (y * width + x) * 3;
    data[i + 0] = r;
    data[i + 1] = g;
    data[i + 2] = b;
  }
  void save(char const *name) {
    std::ofstream out(name, std::ios::binary);
    out << "P6 " << width << ' ' << height << " 255\n";
    out << data << '\n';
  }
};

namespace Vector {

struct small: std::array<float, 3> {};

struct vec: std::array<double, 3> {
  operator small() const {
    return small { (float)(*this)[0], (float)(*this)[1], (float)(*this)[2] };
  }
};

vec &operator+=(vec &u, vec const &v) {
  for (int i = 0; i < 3; ++i) { u[i] += v[i]; }
  return u;
}

vec operator+(vec const &u, vec const &v) {
  vec w(u);
  return w += v;
}

vec operator-(vec const &u, vec const &v) {
  vec w;
  for (int i = 0; i < 3; ++i) { w[i] = u[i] - v[i]; }
  return w;
}

vec operator-(vec const &u) {
  vec w;
  for (int i = 0; i < 3; ++i) { w[i] = -u[i]; }
  return w;
}

vec operator*(double a, vec const &u) {
  vec w;
  for (int i = 0; i < 3; ++i) { w[i] = a * u[i]; }
  return w;
}

double operator|(vec const &u, vec const &v) {
  double r = 0;
  for (int i = 0; i < 3; ++i) { r += u[i] * v[i]; }
  return r;
}

double norm(vec const &u) {
  return sqrt(u|u);
}

vec normalize(vec const &u) {
  return (1 / norm(u)) * u;
}

vec mix(vec const &u, vec const &v, double k) {
  return (1 - k) * u + k * v;
}

vec operator*(vec const &u, vec const &v) {
  vec w;
  for (int i = 0; i < 3; ++i) { w[i] = u[i] * v[i]; }
  return w;
}

vec cross(vec const &u, vec const &v) {
  vec w;
  w[0] = u[1] * v[2] - u[2] * v[1];
  w[1] = u[2] * v[0] - u[0] * v[2];
  w[2] = u[0] * v[1] - u[1] * v[0];
  return w;
}

}

using vec = Vector::vec;

namespace Matrix {

struct mat: std::array<double, 9> {};

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

}

using mat = Matrix::mat;

namespace Spectrum {

int const min_wl = 360, max_wl = 750;
int const nb_wl = 390;

struct full_spectrum: std::array<double, nb_wl> {
  full_spectrum() { fill(0.); }
};

double bbaux(int temp, int wl) {
  double c = 299792458.;
  double h = 6.62606957e-34;
  double kb = 1.3806488e-23;
  double l = wl * 1e-9;
  double l2 = l * l, l5 = l2 * l2 * l;
  return 2 * h * c * c / (l5 * (exp((h * c) / (l * kb * temp)) - 1));
}

full_spectrum blackbody(int temp) {
  full_spectrum sp;
  int ml = 2.8977721e6 / temp;
  if (ml > max_wl) ml = max_wl;
  double strength = 1 / bbaux(temp, ml);
  for (int i = 0; i < nb_wl; ++i) {
    double t = (double)i / nb_wl;
    sp[i] = bbaux(temp, (1 - t) * min_wl + t * max_wl) * strength;
  }
  return sp;
}

double const sumY = 106.946;

// Wyman, Sloane, Shirley, JCGT 2013
double wss(double w, double c1, double c2) {
  double v = w * (w < 0 ? c1 : c2);
  return exp(-0.5 * v * v);
}

vec toXYZ(double wl) {
  return {
    0.362  * wss(wl - 442.0, 0.0624, 0.0374) +
    1.056  * wss(wl - 599.8, 0.0264, 0.0323) +
    -0.065 * wss(wl - 501.1, 0.0490, 0.0382),
    0.821  * wss(wl - 568.8, 0.0213, 0.0247) +
    0.286  * wss(wl - 530.9, 0.0613, 0.0322),
    1.217  * wss(wl - 437.0, 0.0845, 0.0278) +
    0.681  * wss(wl - 459.0, 0.0385, 0.0725)
  };
}

vec toXYZ(full_spectrum const &s) {
  vec c { 0., 0., 0. };
  for (int i = 0; i < nb_wl; ++i) {
    double t = (double)i / nb_wl;
    c += s[i] * toXYZ((1 - t) * min_wl + t * max_wl);
  }
  return (1 / sumY) * c;
}

mat RGBtoXYZ = {
  0.49000, 0.31000, 0.20000,
  0.17697, 0.81240, 0.01063,
  0.00000, 0.01000, 0.99000
};

mat XYZtoRGB = {
  2.36461385,  -0.89654057, -0.46807328,
  -0.51516621, 1.4264081,   0.0887581,
  0.0052037,   -0.01440816, 1.00920446
};

}

int clamp(double v) {
  if (v <= 0) return 0;
  else if (v >= 1) return 255;
  else return v * 255;
}

int main() {
  int const w = 400, s = 15, g = 20;
  Spectrum::full_spectrum white = Spectrum::blackbody(5800);
  image img(w, w);
  for (int y = 0; y < w; ++y) {
    for (int x = 0; x <= y; ++x) {
      double yy = w - y;
      double xx = x / yy;
      double zz = (y - x) / yy;
      vec c { xx, 1., zz };
      c = Spectrum::XYZtoRGB * c;
      c = 1 / (std::max(std::max(c[0], c[1]), c[2])) * c;
      if (c[0] < 0 || c[1] < 0 || c[2] < 0) continue;
      img.write(x, y, clamp(c[0]), clamp(c[1]), clamp(c[2]));
    }
  }
  std::array<vec, s> xyz;
  for (int i = 0; i < s; ++i) {
    Spectrum::full_spectrum sp;
    for (int j = 0; j < Spectrum::nb_wl; ++j) {
      int k = j * s / Spectrum::nb_wl;
      sp[j] = (k == i) ? 1. : 0.;
    }
    xyz[i] = Spectrum::toXYZ(sp);
  }
  std::array<int, g * g> grid;
  grid.fill(0);
  for (int i = 1; i < (1 << (2 * s)); ++i) {
    vec c { 0., 0., 0. };
    for (int j = 0; j < s; ++j) {
      c += (3 & (i >> (2 * j))) * xyz[j];
    }
    double cc = (g - 1) / (c[0] + c[1] + c[2]);
    int x = c[0] * cc, y = g - c[1] * cc;
    assert(0 <= x && x < g && 0 <= y && y < g);
    int &v = grid[y * g + x];
    int ni = 0, nv = 0;
    for (int j = 0; j < s; ++j) {
      ni += (3 & (i >> (2 * j)));
      nv += (3 & (v >> (2 * j)));
    }
    if (ni > nv) v = i;
  }
  std::cout << "int palette[" << g << "][" << g << "] = {\n" << std::showbase << std::hex;
  for (int y = 0; y < g; ++y) {
    std::cout << "  { ";
    for (int x = 0; x < g; ++x) {
      int i = grid[y * g + x];
      std::cout << i << ", ";
      if (i == 0) continue;
      Spectrum::full_spectrum sp;
      for (int j = 0; j < Spectrum::nb_wl; ++j) {
        int k = j * s / Spectrum::nb_wl;
        sp[j] = white[j] / 3 * (3 & (i >> (2 * k)));
      }
      vec c = Spectrum::XYZtoRGB * Spectrum::toXYZ(sp);
      for (int yy = 0; yy < 10; ++yy) {
        for (int xx = 0; xx < 10; ++xx) {
          img.write(w/2 + x * 10 + xx, y * 10 + yy, clamp(c[0]), clamp(c[1]), clamp(c[2]));
        }
      }
    }
    std::cout << "},\n";
  }
  std::cout << "};\n";
  for (int i = 0; i < g * 10; ++i) {
    img.write(w/2, i, 255, 255, 255);
    img.write(w/2 + g * 5, i, 255, 255, 255);
    img.write(w/2 + g * 10 - 1, i, 255, 255, 255);
    img.write(w/2 + i, 0, 255, 255, 255);
    img.write(w/2 + i, g * 5, 255, 255, 255);
    img.write(w/2 + i, g * 10, 255, 255, 255);
  }
  img.save("palette.ppm");
  return 0;
}
