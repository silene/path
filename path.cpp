#include <algorithm>
#include <atomic>
#include <cassert>
#include <fstream>
#include <functional>
#include <iostream>
#include <random>
#include <sstream>
#include <string>
#include <syncstream>
#include <thread>
#include <vector>

namespace Settings {

enum skind { None, Plain, Weighted };

#define SETTINGS
#include USERDATA
#undef SETTINGS

}

namespace Debug {

std::atomic<int64_t> samples = 0;
std::atomic<int64_t> rays = 0;
std::atomic<int64_t> irays = 0;
std::atomic<int64_t> boxes = 0;
std::atomic<int64_t> solids = 0;

}

struct debug {
  int samples;
  int rays;
  int irays;
  int boxes;
  int solids;
};

struct image {
  int width, height;
  std::string data;
  image(int w, int h): width(w), height(h), data(w * h * 3, '\0') {}
  void write(int x, int y, int r, int g, int b) {
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

thread_local std::mt19937_64 *rng;
thread_local debug *dbg;

void spawn(std::function<void(int)> const &f, int n) {
  std::atomic<int> idx = n - 1;
  auto worker = [&]() {
    std::mt19937_64 gen;
    debug d;
    rng = &gen;
    dbg = &d;
    for (;;) {
      int i = idx.fetch_sub(1);
      if (i < 0) break;
      std::osyncstream(std::cout) << ((n - 1 - i) * 100 / n) << "%\r" << std::flush;
      d = { 0, 0, 0, 0, 0 };
      f(i);
      Debug::samples += d.samples;
      Debug::rays += d.rays;
      Debug::irays += d.irays;
      Debug::boxes += d.boxes;
      Debug::solids += d.solids;
    }
  };
  int c = std::thread::hardware_concurrency();
  if (c == 0) worker();
  else {
    std::vector<std::thread> threads;
    for (int i = 0; i < c; ++i) { threads.emplace_back(worker); }
    for (auto &t: threads) { t.join(); };
  }
}

template<class T>
using biased = std::pair<T, double>;

// Computations should be invariant wrt this arbitrarily large value.
// Taken as 1 to avoid numerical issues.
double const Dirac = 1.;

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

vec &operator*=(vec &u, double a) {
  for (int i = 0; i < 3; ++i) { u[i] *= a; }
  return u;
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

struct sphere_uniform_sampler {
  biased<vec> sample() const;
  double pdf(vec const &) const;
};

biased<vec> sphere_uniform_sampler::sample() const {
  std::uniform_real_distribution dis(-1., 1.);
  for (;;) {
    vec w;
    for (int i = 0; i < 3; ++i) { w[i] = dis(*rng); }
    double l = norm(w);
    if (l > 1 || l < 1e-6) continue;
    return { (1 / l) * w, 0.25 * M_1_PI };
  }
}

double sphere_uniform_sampler::pdf(vec const &) const {
  return 0.25 * M_1_PI;
}

struct hemisphere_uniform_sampler {
  vec u;
  hemisphere_uniform_sampler(vec const &u_): u(u_) {}
  biased<vec> sample() const;
  double pdf(vec const &) const;
};

biased<vec> hemisphere_uniform_sampler::sample() const {
  std::uniform_real_distribution dis(-1., 1.);
  for (;;) {
    vec w;
    for (int i = 0; i < 3; ++i) { w[i] = dis(*rng); }
    double l = norm(w);
    if (l > 1 || l < 1e-6) continue;
    if ((u | w) < 0) l = -l;
    return { (1 / l) * w, 0.5 * M_1_PI };
  }
}

double hemisphere_uniform_sampler::pdf(vec const &v) const {
  if ((u | v) <= 0) return 0.;
  return 0.5 * M_1_PI;
}

struct hemisphere_linear_sampler {
  vec u;
  hemisphere_linear_sampler(vec const &u_): u(u_) {}
  biased<vec> sample() const;
  double pdf(vec const &) const;
};

biased<vec> hemisphere_linear_sampler::sample() const {
  std::uniform_real_distribution dis(-1., 1.);
  for (;;) {
    vec w;
    for (int i = 0; i < 3; ++i) { w[i] = dis(*rng); }
    double l = norm(w);
    if (l > 1 || l < 1e-6) continue;
    w = (1 / l) * w;
    double s = (u | w);
    if (s < 0) { w = -w; s = -s; }
    // w - s * u is orthogonal to u, of norm sqrt(1 - s²)
    double ss = sqrt(s);
    double f = 1 / sqrt(1 + s);
    // return a unit vector r such that r|u = ss = sqrt(w|u)
    return { f * w + (ss - f * s) * u, ss * M_1_PI };
  }
}

double hemisphere_linear_sampler::pdf(vec const &v) const {
  if ((u | v) <= 0) return 0.;
  return (u | v) * M_1_PI;
}

struct hemisphere_power_sampler {
  vec u;
  double p;
  hemisphere_power_sampler(vec const &u_, double p_): u(u_), p(p_) {}
  biased<vec> sample() const;
  double pdf(vec const &) const;
};

biased<vec> hemisphere_power_sampler::sample() const {
  std::uniform_real_distribution dis(-1., 1.);
  for (;;) {
    vec w;
    for (int i = 0; i < 3; ++i) { w[i] = dis(*rng); }
    double l = norm(w);
    if (l > 1 || l < 1e-6) continue;
    w = (1 / l) * w;
    double s = (u | w);
    if (s < 0) { w = -w; s = -s; }
    l = log(s) / (p + 1);
    double ss = exp(l);
    double f = sqrt((1 - ss * ss) / (1 - s * s));
    // return a unit vector r such that r|u = ss = (w|u)^(1/(p+1))
    return { f * w + (ss - f * s) * u, (p + 1) * exp(l * p) * 0.5 * M_1_PI };
  }
}

double hemisphere_power_sampler::pdf(vec const &v) const {
  if ((u | v) <= 0) return 0.;
  return (p + 1) * exp(log(u | v) * p) * 0.5 * M_1_PI;
}

struct cone_uniform_sampler {
  vec u;
  double cmax;
  cone_uniform_sampler(vec const &u_, double c_): u(u_), cmax(c_) {}
  biased<vec> sample() const;
  double pdf(vec const &) const;
};

biased<vec> cone_uniform_sampler::sample() const {
  std::uniform_real_distribution dis(-1., 1.);
  std::uniform_real_distribution dis2(0., 1.);
  for (;;) {
    vec w;
    for (int i = 0; i < 3; ++i) { w[i] = dis(*rng); }
    vec v = cross(u, w);
    double d = norm(v);
    if (d < 1e-6) continue;
    double c = 1 + dis2(*rng) * (cmax - 1);
    v *= sqrt(1 - c * c) / d;
    return { v + c * u, 1 / (2 * M_PI * (1 - cmax)) };
  }
}

double cone_uniform_sampler::pdf(vec const &v) const {
  if ((u | v) <= cmax) return 0.;
  return 1 / (2 * M_PI * (1 - cmax));
}

}

using vec = Vector::vec;

#if 0
int main() {
  std::mt19937_64 gen;
  rng = &gen;
  int nb = 100000000;
  double v = 0;
  vec u { 1., 0., 0. };
  Vector::hemisphere_uniform_sampler t(u);
  Vector::hemisphere_power_sampler s(u, 2.1);
  //Vector::hemisphere_linear_sampler s(u);
  for (int i = 0; i < nb; ++i) {
    auto [vv, pp] = t.sample();
    v += s.pdf(vv) / pp;
  }
  std::cout << v / nb << '\n';
}
#endif

namespace Box {

struct box {
  Vector::small p1, p2;
};

box const empty { { INFINITY, INFINITY, INFINITY }, { -INFINITY, -INFINITY, -INFINITY } };

box const whole { { -INFINITY, -INFINITY, -INFINITY }, { INFINITY, INFINITY, INFINITY } };

bool intersect(vec const &pos, vec const &dir, box const &b, double dmax) {
  ++dbg->boxes;
  double v1 = 0, v2 = dmax;
  for (int i = 0; i < 3; ++i) {
    double d = 1 / dir[i];
    double t1 = (b.p1[i] - pos[i]) * d;
    double t2 = (b.p2[i] - pos[i]) * d;
    if (t2 < t1) std::swap(t1, t2);
    if (t1 > v1) v1 = t1;
    if (t2 < v2) v2 = t2;
    if (v2 < v1) return false;
  }
  return true;
}

void inflate(box &b, vec const &v) {
  for (int i = 0; i < 3; ++i) { b.p1[i] = std::min<float>(b.p1[i], v[i]); }
  for (int i = 0; i < 3; ++i) { b.p2[i] = std::max<float>(b.p2[i], v[i]); }
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

using box = Box::box;

namespace Ball {

struct ball {
  vec center;
  double radius;
};

ball const empty { vec(), -INFINITY };

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

using ball = Ball::ball;

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

namespace Transform {

struct iso {
  vec center;
  double scale, iscale;
  mat rotate;
  iso(vec const &c, double sc, vec const &a, double r):
    center(c), scale(sc), iscale(1 / sc), rotate(Matrix::rotation(a, r)) {}

  vec position(vec const &pos) const {
    return iscale * rotate * (pos - center);
  }
};

using ptr = iso const *;

}

namespace Spectrum {

int const min_wl = 360, max_wl = 750;
int const nb_wl = 390;
int const nb_s = 4;

struct sampled_spectrum: std::array<float, nb_s> {
  sampled_spectrum() = default;
  sampled_spectrum(double d) { fill(d); }
  bool zero() const {
    for (int i = 0; i < nb_wl; ++i) {
      if ((*this)[i]) return false;
    }
    return true;
  }
};

sampled_spectrum &operator+=(sampled_spectrum &u, sampled_spectrum const &v) {
  for (int i = 0; i < nb_s; ++i) { u[i] += v[i]; }
  return u;
}

sampled_spectrum operator+(sampled_spectrum const &u, sampled_spectrum const &v) {
  sampled_spectrum w(u);
  return w += v;
}

sampled_spectrum &operator*=(sampled_spectrum &u, double a) {
  for (int i = 0; i < nb_s; ++i) { u[i] *= a; }
  return u;
}

sampled_spectrum operator*(double a, sampled_spectrum const &u) {
  sampled_spectrum w;
  for (int i = 0; i < nb_s; ++i) { w[i] = a * u[i]; }
  return w;
}

sampled_spectrum operator*(sampled_spectrum const &u, sampled_spectrum const &v) {
  sampled_spectrum w;
  for (int i = 0; i < nb_s; ++i) { w[i] = u[i] * v[i]; }
  return w;
}

struct sampled_wl {
  std::array<double, nb_s> lambda, pdf;
  sampled_wl(double s) {
    for (int i = 0; i < nb_s; ++i) {
      double l = s + i * (double)(max_wl - min_wl) / nb_s;
      if (l >= max_wl) l = l - max_wl + min_wl;
      assert(min_wl <= l && l < max_wl);
      lambda[i] = l;
      pdf[i] = 1. / (max_wl - min_wl);
    }
  }
};

struct full_spectrum: std::array<float, nb_wl> {
  full_spectrum() { fill(0.); }
};

void add(full_spectrum &fs, sampled_wl const &wl, sampled_spectrum const &s) {
  for (int i = 0; i < nb_s; ++i) {
    int l = (wl.lambda[i] - min_wl) / (max_wl - min_wl) * nb_wl;
    assert(0 <= l && l < nb_wl);
    fs[l] += s[i] / wl.pdf[i];
  }
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

double toY(double wl) {
  double v =
    0.821  * wss(wl - 568.8, 0.0213, 0.0247) +
    0.286  * wss(wl - 530.9, 0.0613, 0.0322);
  return v / sumY;
}

vec toXYZ(full_spectrum const &s) {
  vec c { 0., 0., 0. };
  for (int i = 0; i < nb_wl; ++i) {
    double t = (double)i / nb_wl;
    c += s[i] * toXYZ((1 - t) * min_wl + t * max_wl);
  }
  return (1 / sumY) * c;
}

struct spectrum_base {
  virtual sampled_spectrum sample(sampled_wl const &) const = 0;
  virtual ~spectrum_base() = default;
};

struct blackbody: spectrum_base {
  double strength;
  int temp;
  //full_spectrum fsp;

  blackbody(double f, int t): temp(t) {
    int ml = 2.8977721e6 / temp;
    if (ml > max_wl) ml = max_wl;
    strength = f / aux(ml);
  }

  double aux(double wl) const {
    double c = 299792458.;
    double h = 6.62606957e-34;
    double kb = 1.3806488e-23;
    double l = wl * 1e-9;
    double l2 = l * l, l5 = l2 * l2 * l;
    return 2 * h * c * c / (l5 * (exp((h * c) / (l * kb * temp)) - 1));
  }

  sampled_spectrum sample(sampled_wl const &l) const {
    sampled_spectrum s;
    for (int i = 0; i < nb_s; ++i) {
      s[i] = aux(l.lambda[i]) * strength;
    }
    return s;
  }
};

struct uniform: spectrum_base {
  double strength;
  uniform(double s): strength(s) {}
  sampled_spectrum sample(sampled_wl const &) const {
    return sampled_spectrum(strength);
  }
};

struct from_wl: spectrum_base {
  double strength, fact;
  double lambda, width;
  from_wl(double s, double l, double w)
    : strength(s), lambda(l), width(w) {
    fact = 1.; // / (width * sqrt(2 * M_PI));
  }
  sampled_spectrum sample(sampled_wl const &wl) const {
    sampled_spectrum s;
    double w = max_wl - min_wl;
    for (int i = 0; i < nb_s; ++i) {
      double l = wl.lambda[i];
      if (l - lambda > 0.5 * w) { l -= w; }
      else if (lambda - l > 0.5 * w) { l += w; }
      double f = (l - lambda) / width;
      s[i] = fact * exp(-0.5 * f * f);
    }
    return s;
  }
};

// 15*2-bit palette generated by the palette.cpp helper.
int palette[20][20] = {
  { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, },
  { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, },
  { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, },
  { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, },
  { 0, 0, 0x3c003003, 0x31002000, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, },
  { 0, 0, 0x3f003403, 0x3e00300b, 0x3f007003, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, },
  { 0, 0x3c003c0b, 0x3e00380f, 0x3f00740f, 0x3f00b00f, 0x3e00f00f, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, },
  { 0, 0x3e002c0f, 0x3f403c0f, 0x3f807c0f, 0x3fc0fc0f, 0x3fc0f80f, 0x3f41f00f, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, },
  { 0x3c000c03, 0x3f001c0f, 0x3f007d0f, 0x3f00fd0f, 0x3fc0fd0f, 0x3ee0fc0f, 0x3fc3fc0f, 0x3fd2f40f, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, },
  { 0x3d000c0b, 0x3c003d0f, 0x3e803d0f, 0x3f40fe0f, 0x3fc1fe0f, 0x3f43fe0f, 0x3fe3fd0f, 0x3febfc0f, 0x3fdffc0b, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, },
  { 0x3c000c0f, 0x3c003e0f, 0x3f803e0f, 0x3f80ff0f, 0x3fa0ff0f, 0x3fd3ff0f, 0x3fcbff0f, 0x3ffbfe0f, 0x3ffffc1f, 0x3fffa80f, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, },
  { 0, 0x3c003f0b, 0x3f007f1f, 0x3fc0ff2f, 0x3fc2ff2f, 0x3fe3ff2f, 0x3fdfff2f, 0x3ffffd3f, 0x3fffbc2f, 0x3fff6c0f, 0x3fff480f, 0, 0, 0, 0, 0, 0, 0, 0, 0, },
  { 0, 0x3d002f0f, 0x3f007e3f, 0x3fc0fe7f, 0x3fc3ffaf, 0x3fdbffaf, 0x3fffff7f, 0x3fffff3f, 0x3fffcd3f, 0x3fff882f, 0x3fff0c0f, 0x3ffe040f, 0, 0, 0, 0, 0, 0, 0, 0, },
  { 0, 0x3f001f0f, 0x3f803f3f, 0x3f90ffef, 0x3fe2ffff, 0x3fdfffff, 0x3fffffff, 0x3fffdf7f, 0x3fff5d3f, 0x3fff083f, 0x3fff002f, 0x3ffd040f, 0x3ffd000f, 0, 0, 0, 0, 0, 0, 0, },
  { 0, 0x3c001f1f, 0x3f007fbf, 0x3f90ffff, 0x3ff1ffff, 0x3feebfff, 0x3ffebfff, 0x3fffcbbf, 0x3fff0e3f, 0x3ffe083f, 0x3ffe002f, 0x3ffd001f, 0x3ffc000f, 0x3ff00007, 0, 0, 0, 0, 0, 0, },
  { 0, 0x3e000f1f, 0x3f803fff, 0x3fa07fff, 0x3ff47fff, 0x3ff93fff, 0x3ffe2fff, 0x3fff4aff, 0x3fff01bf, 0x3ffd013f, 0x3ffc002f, 0x3ff8001f, 0x3ee0000f, 0, 0, 0, 0, 0, 0, 0, },
  { 0, 0x3c000f4f, 0x3f802fff, 0x3fe02fff, 0x3fe42fff, 0x3ffc1fff, 0x3ffc4fff, 0x3fff03ff, 0x3ffd00bf, 0x3ffc006f, 0x3ff0001f, 0, 0, 0, 0, 0, 0, 0, 0, 0, },
  { 0, 0x30000001, 0x3f801fff, 0x3fa00fff, 0x3fb40fff, 0x3ff40fff, 0x3ffc07ff, 0x3ffc01ff, 0x3ff0003f, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, },
  { 0x30000003, 0x30000002, 0x3f400bff, 0x3f900bff, 0x3ef007ff, 0x3fe403ff, 0x3ff000bf, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, },
  { 0x3, 0, 0x3f0003ff, 0x3ec003ff, 0x3f9000ff, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, },
};

struct from_palette: spectrum_base {
  double strength;
  int val;
  from_palette(double s, int x, int y)
    : strength(s), val(palette[y][x]) {}
  sampled_spectrum sample(sampled_wl const &wl) const {
    sampled_spectrum s;
    for (int i = 0; i < nb_s; ++i) {
      int l = 15 * (wl.lambda[i] - min_wl) / (max_wl - min_wl);
      s[i] = (3 & (val >> (2 * l))) * strength / 3;
    }
    return s;
  }
};

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

using ptr = Spectrum::spectrum_base const *;

}

using sampled_wl = Spectrum::sampled_wl;
using sampled_spectrum = Spectrum::sampled_spectrum;

namespace Solid {

struct base {
  virtual double distance(vec const &, vec const &, int) const = 0;
  virtual vec normal(vec const &, int) const = 0;
  virtual int subparts() const { return 1; }
  virtual box bounds(int, Transform::ptr) const = 0;
  virtual ball sbounds(int, Transform::ptr) const = 0;
  virtual bool inside(vec const &) const = 0;
  virtual bool valid(vec const &, int p) const
  { assert(p == 0); return true; }
  virtual ~base() = default;
};

using ptr = base const *;

struct sphere: base {
  vec center;
  double radius;
  sphere(vec const &c, double r):
    center(c), radius(r) {}

  double distance(vec const &pos, vec const &dir, int) const {
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

  vec normal(vec const &pos, int) const {
    return normalize(pos - center);
  }

  box bounds(int, Transform::ptr t) const {
    assert(t == NULL);
    vec r { radius, radius, radius };
    return { center - r, center + r };
  }

  ball sbounds(int, Transform::ptr t) const {
    assert(t == NULL);
    return { center, radius };
  }

  bool inside(vec const &pos) const {
    vec p = pos - center;
    return (p | p) - radius * radius <= 1e-10;
  }
};

struct plane: base {
  vec normal_;
  double dist;
  plane(vec const &n, double d): normal_(n), dist(d) {}

  double distance(vec const &pos, vec const &dir, int) const {
    double d = (pos | normal_) + dist;
    double s = (dir | normal_);
    if (s == 0) return d == 0 ? 0 : INFINITY;
    double r = - d / s;
    return r < 0 ? INFINITY : r;
  }

  vec normal(vec const &, int) const {
    return normal_;
  }

  box bounds(int, Transform::ptr t) const {
    assert(t == NULL);
    box b = Box::empty;
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
        inflate(b, v);
      }
    }
    return b;
  }

  ball sbounds(int, Transform::ptr t) const {
    assert(false);
  }

  bool inside(vec const &pos) const {
    return (pos | normal_) + dist <= 1e-10;
  }
};

std::pair<double, double> toUV(vec p, vec q, vec r) {
  // the barycentric coordinates u, v, 1-u-v are proportional
  // to the areas of the three subtriangles opr, oqr, pqr
  vec n = cross(p, q);
  double d = (n | n);
  double u = (n | cross(r, q)) / d;
  double v = (n | cross(p, r)) / d;
  return { u, v };
}

struct mesh: base {
  std::vector<vec> vertices;
  std::vector<vec> normals;
  std::vector<std::array<int, 3>> facets;
  bool inv_normal;

  mesh(char const *name, bool = false);
  int subparts() const { return facets.size(); }

  double distance(vec const &pos, vec const &dir, int data) const {
    std::array<int, 3> const &f = facets[data];
    vec const &p0 = vertices[f[0]], &p1 = vertices[f[1]], &p2 = vertices[f[2]];
    vec p10 = p1 - p0, p20 = p2 - p0, p = pos - p0;
    vec n = cross(p10, p20);
    double d = p | n;
    double s = dir | n;
    double r;
    if (s == 0) {
      if (d != 0) return INFINITY;
      r = 0;
    } else {
      r = - d / s;
      if (r < 0) return INFINITY;
    }
    p += r * dir;
    auto [u, v] = toUV(p10, p20, p);
    if (u < 0 || v < 0 || u + v > 1) return INFINITY;
    return r;
  }

  vec normal(vec const &pos, int data) const;

  box bounds(int d, Transform::ptr t) const {
    assert(t);
    std::array<int, 3> const &f = facets[d];
    vec const &p0 = vertices[f[0]], &p1 = vertices[f[1]], &p2 = vertices[f[2]];
    mat m = t->scale * transpose(t->rotate);
    vec v = m * p0 + t->center;
    box b { v, v };
    v = m * p1 + t->center;
    inflate(b, v);
    v = m * p2 + t->center;
    inflate(b, v);
    return b;
  }

  ball sbounds(int, Transform::ptr t) const {
    assert(false);
  }

  bool inside(vec const &) const {
    assert(false);
  }
};

mesh::mesh(char const *name, bool b)
  : inv_normal(b) {
  std::ifstream file(name);
  for (std::string line; std::getline(file, line); ) {
    std::string::size_type n = line.find(' ');
    if (n == std::string::npos) continue;
    if (line[0] == 'v' && n == 1) {
      std::istringstream l(line);
      l.seekg(2);
      vertices.resize(vertices.size() + 1);
      vec &v = vertices.back();
      l >> v[0] >> v[1] >> v[2];
    } else if (line[0] == 'f' && n == 1) {
      std::istringstream l(line);
      l.seekg(2);
      facets.resize(facets.size() + 1);
      auto &f = facets.back();
      l >> f[0] >> f[1] >> f[2];
      for (int i = 0; i < 3; ++i) {
        int &v = f[i];
        if (v > 0) --v;
        else v = vertices.size() + v;
        assert(0 <= v && (unsigned)v < vertices.size());
      }
      //vec const &v0 = vertices[f[0]];
      //double n = norm(cross(vertices[f[1]] - v0, vertices[f[2]] - v0));
      //if (n <= 0.99) facets.pop_back();
    } else if (line[0] == 'v' && line[1] == 'n' && n == 2) {
      std::istringstream l(line);
      l.seekg(3);
      normals.resize(normals.size() + 1);
      vec &v = normals.back();
      l >> v[0] >> v[1] >> v[2];
    } else continue;
  }
  if (normals.size() < vertices.size()) normals.clear();
}

vec mesh::normal(vec const &pos, int data) const {
  std::array<int, 3> const &f = facets[data];
  vec const &p0 = vertices[f[0]], &p1 = vertices[f[1]], &p2 = vertices[f[2]];
  vec p10 = p1 - p0, p20 = p2 - p0;
  if (normals.empty()) {
    vec n = normalize(cross(p10, p20));
    return inv_normal ? -n : n;
  }
  auto [u, v] = toUV(p10, p20, pos - p0);
  vec const &n0 = normals[f[0]], &n1 = normals[f[1]], &n2 = normals[f[2]];
  return normalize((1 - u - v) * n0 + u * n1 + v * n2);
}

struct union_: base {
  ptr obj1, obj2;
  int parts1, parts;

  union_(ptr o1, ptr o2)
    : obj1(o1), obj2(o2), parts1(obj1->subparts()),
      parts(parts1 + obj2->subparts()) {}

  double distance(vec const &pos, vec const &dir, int p) const {
    if (p < parts1) return obj1->distance(pos, dir, p);
    return obj2->distance(pos, dir, p - parts1);
  }

  vec normal(vec const &pos, int p) const {
    if (p < parts1) return obj1->normal(pos, p);
    return obj2->normal(pos, p - parts1);
  }

  int subparts() const { return parts; }

  box bounds(int p, Transform::ptr t) const {
    assert(t == NULL);
    if (p < parts1) return obj1->bounds(p, t);
    return obj2->bounds(p - parts1, t);
  }

  ball sbounds(int, Transform::ptr) const {
    assert(false);
  }

  bool inside(vec const &pos) const {
    return obj1->inside(pos) || obj2->inside(pos);
  }

  bool valid(vec const &pos, int p) const {
    if (p < parts1) return obj1->valid(pos, p) && !obj2->inside(pos);
    return obj2->valid(pos, p - parts1) && !obj1->inside(pos);
  }

};

struct intersection: base {
  ptr obj1, obj2;
  int parts1, parts;
  box bnds;

  intersection(ptr o1, ptr o2)
    : obj1(o1), obj2(o2), parts1(obj1->subparts()),
      parts(parts1 + obj2->subparts())
  {
    box b = Box::empty;
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

  double distance(vec const &pos, vec const &dir, int p) const {
    if (p < parts1) return obj1->distance(pos, dir, p);
    return obj2->distance(pos, dir, p - parts1);
  }

  vec normal(vec const &pos, int p) const {
    if (p < parts1) return obj1->normal(pos, p);
    return obj2->normal(pos, p - parts1);
  }

  int subparts() const { return parts; }

  box bounds(int p, Transform::ptr t) const {
    assert(t == NULL);
    if (p < parts1) return Box::intersect(bnds, obj1->bounds(p, t));
    return Box::intersect(bnds, obj2->bounds(p - parts1, t));
  }

  ball sbounds(int, Transform::ptr) const {
    assert(false);
  }

  bool inside(vec const &pos) const {
    return obj1->inside(pos) && obj2->inside(pos);
  }

  bool valid(vec const &pos, int p) const {
    if (p < parts1) return obj1->valid(pos, p) && obj2->inside(pos);
    return obj2->valid(pos, p - parts1) && obj1->inside(pos);
  }

};

struct exclusion: base {
  ptr obj1, obj2;
  int parts1, parts;
  box bnds1;

  exclusion(ptr o1, ptr o2)
    : obj1(o1), obj2(o2), parts1(obj1->subparts()),
      parts(parts1 + obj2->subparts())
  {
    box b = Box::empty;
    for (int n = 0; n < parts1; ++n) {
      b = Box::merge(b, obj1->bounds(n, NULL));
    }
    bnds1 = b;
  }

  double distance(vec const &pos, vec const &dir, int p) const {
    if (p < parts1) return obj1->distance(pos, dir, p);
    return obj2->distance(pos, dir, p - parts1);
  }

  vec normal(vec const &pos, int p) const {
    if (p < parts1) return obj1->normal(pos, p);
    return -obj2->normal(pos, p - parts1);
  }

  int subparts() const { return parts; }

  box bounds(int p, Transform::ptr t) const {
    assert(t == NULL);
    if (p < parts1) return obj1->bounds(p, t);
    return Box::intersect(bnds1, obj2->bounds(p - parts1, t));
  }

  ball sbounds(int, Transform::ptr) const {
    assert(false);
  }

  bool inside(vec const &pos) const {
    return obj1->inside(pos) && !obj2->inside(pos);
  }

  bool valid(vec const &pos, int p) const {
    if (p < parts1) return obj1->valid(pos, p) && !obj2->inside(pos);
    return obj2->valid(pos, p - parts1) && obj1->inside(pos);
  }

};

}

struct intersection {
  vec pos, normal, out;
};

struct path_point {
  intersection pt;
  vec inc;
  double pdf;
  bool already_illuminated;
};

namespace Light {

struct ray {
  vec dir;
  double dist;
};

struct base {
  bool has_pdf, surrounding;
  base(bool b, bool s): has_pdf(b), surrounding(s) {}
  virtual biased<ray> sample(vec const &pos, vec const &n) const = 0;
  virtual double pdf(vec const &pos, vec const &n, vec const &dir) const { assert(false); }
  virtual sampled_spectrum get_sp(vec const &pos, vec const &dir, sampled_wl const &wl) const = 0;
  virtual ~base() = default;
};

struct point: base {
  Spectrum::ptr sp;
  vec pos;

  point(Spectrum::ptr s, vec const &p)
    : base(false, false), sp(s), pos(p) {}

  biased<ray> sample(vec const &p, vec const &) const {
    vec d = pos - p;
    double dist = norm(d);
    d = (1 / dist) * d;
    return { { d, dist }, Dirac };
  }

  sampled_spectrum get_sp(vec const &p, vec const &, sampled_wl const &wl) const {
    vec d = pos - p;
    double f = 1 / (d | d);
    return Dirac * f * sp->sample(wl);
  }
};

struct spot: base {
  Spectrum::ptr sp;
  vec pos, dir;
  double angle1, angle2;

  spot(Spectrum::ptr s, vec const &p, vec const &d, double a1, double a2)
    : base(false, false), sp(s), pos(p), dir(d), angle1(a1), angle2(a2) {}

  biased<ray> sample(vec const &p, vec const &) const {
    vec d = pos - p;
    double dist = norm(d);
    d = (1 / dist) * d;
    return { { d, dist }, Dirac };
  }

  sampled_spectrum get_sp(vec const &p, vec const &, sampled_wl const &wl) const {
    vec d = pos - p;
    double dist = norm(d);
    d = (1 / dist) * d;
    double g = - (d | dir);
    if (g <= angle2) return sampled_spectrum(0.);
    double f = 1 / (dist * dist);
    if (g <= angle1) f = f * (g - angle2) / (angle1 - angle2);
    return Dirac * f * sp->sample(wl);
  }
};

struct directional: base {
  Spectrum::ptr sp;
  vec dir;

  directional(Spectrum::ptr s, vec const &d)
    : base(false, false), sp(s), dir(-d) {}

  biased<ray> sample(vec const &, vec const &) const {
    return { { dir, INFINITY }, Dirac };
  }

  sampled_spectrum get_sp(vec const &, vec const &, sampled_wl const &wl) const {
    return Dirac * sp->sample(wl);
  }
};

struct multidirectional: base {
  Spectrum::ptr sp;
  vec dir;
  double cmax, inv_area;

  multidirectional(Spectrum::ptr s, vec const &d, double a)
    : base(true, true), sp(s), dir(d), cmax(cos(a)),
      inv_area(1 / (2 * M_PI * (1 - cmax))) {}

  biased<ray> sample(vec const &, vec const &) const {
    Vector::cone_uniform_sampler s(dir, cmax);
    auto [dir, p] = s.sample();
    return { { dir, INFINITY }, p };
  }

  double pdf(vec const &, vec const &, vec const &d) const {
    Vector::cone_uniform_sampler s(dir, cmax);
    return s.pdf(d);
  }

  sampled_spectrum get_sp(vec const &, vec const &d, sampled_wl const &wl) const {
    if ((dir | d) <= cmax) return { 0. };
    return inv_area * sp->sample(wl);
  }
};

struct uniform: base {
  Spectrum::ptr sp;

  uniform(Spectrum::ptr s): base(true, true), sp(s) {}

  biased<ray> sample(vec const &, vec const &n) const {
    Vector::hemisphere_linear_sampler s(n);
    auto [d, p] = s.sample();
    return { { d, INFINITY }, p };
  }

  double pdf(vec const &, vec const &n, vec const &d) const {
    Vector::hemisphere_linear_sampler s(n);
    return s.pdf(d);
  }

  sampled_spectrum get_sp(vec const &, vec const &, sampled_wl const &wl) const {
    return sp->sample(wl);
  }
};

struct spherical: base {
  Spectrum::ptr sp;
  Solid::sphere *sph;
  double inv_area;

  spherical(Spectrum::ptr s, vec const &c, double r)
    : base(true, false), sp(s), sph(new Solid::sphere(c, r)),
      inv_area(0.25 * M_1_PI / (r * r)) {}

  biased<ray> sample(vec const &pos, vec const &) const {
    vec v = sph->center - pos;
    double cmax = sqrt(1 - sph->radius * sph->radius / (v | v));
    v = normalize(v);
    Vector::cone_uniform_sampler s(v, cmax);
    auto [dir,p] = s.sample();
    return { { dir, sph->distance(pos, dir, 0) }, p };
  }

  double pdf(vec const &pos, vec const &, vec const &d) const {
    vec v = sph->center - pos;
    double cmax = sqrt(1 - sph->radius * sph->radius / (v | v));
    v = normalize(v);
    Vector::cone_uniform_sampler s(v, cmax);
    return s.pdf(d);
  }

  sampled_spectrum get_sp(vec const &pos, vec const &dir, sampled_wl const &wl) const {
    return inv_area * sp->sample(wl);
  }
};

using ptr = base const *;

}

namespace Material {

enum skind { Specular, SemiSpecular, Diffuse };

struct ray {
  sampled_spectrum sp;
  vec dir;
  skind specular;
};

enum mkind { Solid, Emissive, Transmitive };

struct material_base {
  mkind kind;
  bool has_bxdf, has_pdf;
  material_base(mkind k, bool b, bool p)
    : kind(k), has_bxdf(b), has_pdf(p) {}
  virtual sampled_spectrum bxdf(intersection const &, vec const &inc, sampled_wl const &) const { assert(false); }
  virtual double pdf(intersection const &, vec const &inc) const { assert(false); }
  virtual biased<ray> sample(intersection const &, sampled_wl const &) const = 0;
  virtual material_base const *regularize() const { return this; }
  virtual ~material_base() = default;
};

vec reflect(vec const &normal, vec const &out) {
  return 2 * (out | normal) * normal - out;
}

std::pair<vec, double> refract(vec const &normal, vec const &out, double eta) {
  double cr = out | normal;
  double s2t = (1 - cr * cr) / (eta * eta);
  if (s2t >= 1) return { vec(), 0. };
  double ct = sqrt(1 - s2t);
  double f1 = (eta * cr - ct) / (eta * cr + ct);
  double f2 = (cr - eta * ct) / (cr + eta * ct);
  double f = 0.5 * (f1 * f1 + f2 * f2);
  if (f >= 1) return { vec(), 0. };
  return { (-1/eta) * out + (cr / eta - ct) * normal, 1 - f };
}

struct emissive: material_base {
  Light::ptr light;

  emissive(Light::ptr l)
    : material_base(Emissive, false, false), light(l) {}

  biased<ray> sample(intersection const &, sampled_wl const &wl) const {
    assert(false);
  }
};

struct lambertian: material_base {
  Spectrum::ptr sp;
  lambertian(Spectrum::ptr s)
    : material_base(Solid, true, true), sp(s) {}

  sampled_spectrum bxdf(intersection const &pt, vec const &inc, sampled_wl const &wl) const {
    double f = (inc | pt.normal) * M_1_PI;
    return f * sp->sample(wl);
  }

  double pdf(intersection const &pt, vec const &inc) const {
    Vector::hemisphere_linear_sampler s(pt.normal);
    return s.pdf(inc);
  }

  biased<ray> sample(intersection const &pt, sampled_wl const &wl) const {
    Vector::hemisphere_linear_sampler s(pt.normal);
    auto [inc, pdf] = s.sample();
    double f = (inc | pt.normal) * M_1_PI;
    return { { f * sp->sample(wl), inc, Diffuse }, pdf };
  }
};

struct rough: material_base {
  Spectrum::ptr sp;
  double roughness, alpha, scaling;
  bool specular;

  rough(Spectrum::ptr s, double r)
    : material_base(Solid, r >= 1e-3, r >= 1e-3), sp(s), roughness(r), specular(!has_pdf)
  {
    alpha = (1 - roughness) / roughness;
    alpha = alpha * alpha * 15;
    scaling = M_1_PI * (alpha + 2) / (4 * (1 - exp(-0.5 * M_LN2 * (alpha + 2))));
  }

  sampled_spectrum bxdf(intersection const &pt, vec const &inc, sampled_wl const &wl) const {
    assert(!specular);
    double ci = inc | pt.normal;
    double cm = normalize(inc + pt.out) | pt.normal;
    assert (cm > 0);
    double s = scaling * exp(log(cm) * alpha);
    double v = (1 - roughness) + roughness * ci;
    return v * s * sp->sample(wl);
  }

  double pdf(intersection const &pt, vec const &inc) const {
    assert(!specular);
    Vector::hemisphere_power_sampler samp(pt.normal, alpha);
    vec m = normalize(inc + pt.out);
    return samp.pdf(m) / 4 * (m | pt.out);
  }

  biased<ray> sample(intersection const &pt, sampled_wl const &wl) const {
    if (specular) {
      vec inc = reflect(pt.normal, pt.out);
      return { { sp->sample(wl), inc, Specular }, 1. };
    }
    Vector::hemisphere_power_sampler samp(pt.normal, alpha);
    auto [m, pdf] = samp.sample();
    double c = pt.out | m;
    pdf /= 4 * c;
    if (c < 1e-6) return { { sampled_spectrum(0.), vec() }, 0. };
    vec inc = 2 * c * m - pt.out;
    double ci = inc | pt.normal;
    if (ci < 1e-6) return { { sampled_spectrum(0.), vec() }, 0. };
    double cm = m | pt.normal;
    double s = scaling * exp(log(cm) * alpha);
    double v = (1 - roughness) + roughness * ci;
    return { { v * s * sp->sample(wl), inc, Diffuse }, pdf };
  }

  material_base const *regularize() const {
    double r = roughness;
    r = 0.4 + r * (0.2 + r * 0.4);
    return new rough(sp, r);
  }
};

struct reflective: material_base {
  Spectrum::ptr sp;
  reflective(Spectrum::ptr s): material_base(Solid, false, false), sp(s) {}

  biased<ray> sample(intersection const &pt, sampled_wl const &wl) const {
    vec inc = reflect(pt.normal, pt.out);
    return { { Dirac * sp->sample(wl), inc, Specular }, Dirac };
  }

  material_base const *regularize() const {
    return new rough(sp, 0.4);
  }
};

struct refractive: material_base {
  double eta;
  refractive(double e): material_base(Transmitive, false, false), eta(e) {}

  biased<ray> sample(intersection const &pt, sampled_wl const &wl) const {
    vec n = pt.normal;
    double eta = this->eta;
    double cr = pt.out | n;
    if (cr < 0) { cr = -cr; n = -n; eta = 1/eta; }
    auto [dir, f] = refract(n, pt.out, eta);
    ray r { sampled_spectrum(Dirac), 2 * cr * n - pt.out, Specular };
    if (f == 0.) return { r, Dirac };
    std::uniform_real_distribution dis(0., 1.);
    if (dis(*rng) <= f) { r.dir = dir; }
    else { f = 1 - f; }
    r.sp *= f;
    return { r, Dirac * f };
  }
};

struct thin_refractive: material_base {
  double eta;
  material_base const *mat;
  thin_refractive(double e, material_base const *m)
    : material_base(Solid, m->has_bxdf, m->has_pdf), eta(e), mat(m) {
    assert(m->kind == Solid);
  }

  sampled_spectrum bxdf(intersection const &pt, vec const &inc, sampled_wl const &wl) const {
    auto [_, f] = refract(pt.normal, pt.out, eta);
    assert(f);
    f = f / (2 - f);
    return f * mat->bxdf(pt, inc, wl);
  }

  double pdf(intersection const &pt, vec const &inc) const {
    auto [_, f] = refract(pt.normal, pt.out, eta);
    assert(f);
    f = f / (2 - f);
    return f * mat->pdf(pt, inc);
  }

  biased<ray> sample(intersection const &pt, sampled_wl const &wl) const {
    auto [_, f] = refract(pt.normal, pt.out, eta);
    ray r { sampled_spectrum(Dirac), reflect(pt.normal, pt.out), Specular };
    if (f == 0.) return { r, Dirac };
    f = f / (2 - f);
    std::uniform_real_distribution dis(0., 1.);
    double pdf;
    if (dis(*rng) <= f) {
      std::tie(r, pdf) = mat->sample(pt, wl);
      r.sp *= f;
      pdf *= f;
    } else {
      pdf = Dirac * (1 - f);
      r.sp = sampled_spectrum(pdf);
      r.specular = SemiSpecular;
    }
    return { r, pdf };
  }
};

using ptr = Material::material_base const *;

}

struct object {
  Solid::ptr solid;
  Material::ptr material;
  Transform::ptr transf;
};

#define SCENE
#include USERDATA
#undef SCENE

namespace Light {

struct light_sampler {
  biased<std::pair<Light::ptr, Light::ray>> sample(intersection const &from) {
    int nbl = Scene::lights.size();
    std::uniform_int_distribution dis(0, nbl - 1);
    Light::ptr l = Scene::lights[dis(*rng)];
    auto [r, pdf] = l->sample(from.pos, from.normal);
    pdf /= nbl;
    return { { l, r }, pdf };
  }

  double pdf(Light::ptr l, intersection const &from, vec const &dir) {
    double p = l->pdf(from.pos, from.normal, dir);
    return p / Scene::lights.size();
  }
};

light_sampler lights;

}

namespace Solver {

struct subobject {
  int obj, data;
};

std::vector<subobject> subobjects;

struct boxed_subobject {
  int obj, data;
  box bo;
};

struct split {
  int axis;
  int before, after;
  struct child {
    box bo;
    int sp;
  };
  child left, center, right;
};

std::vector<split> splits;

struct bs_cmp {
  int axis;
  bool operator()(boxed_subobject const &b1, boxed_subobject const &b2) const {
    return b1.bo.p2[axis] < b2.bo.p2[axis];
  }
};

struct bs_part {
  double value;
  int axis;
  bool operator()(boxed_subobject const &b) const {
    return b.bo.p1[axis] < value;
  }
};

split::child split_objs(std::vector<boxed_subobject> &bs, int ib, int ie, int ax, int axm) {
  if (ie - ib <= 1) {
    box b = Box::empty;
    for (int i = ib; i < ie; ++i) {
      b = merge(b, bs[i].bo);
    }
    return { b, -1 };
  }
  std::sort(&bs[ib], &bs[ie], bs_cmp { ax });
  int im = (ib + ie) / 2;
  double v = bs[im - 1].bo.p2[ax];
  int in = std::partition(&bs[im], &bs[ie], bs_part { v, ax }) - &bs[0];
  if (ax != axm && 4 * in >= ib + 3 * ie) {
    if (axm < 0) axm = ax;
    if (++ax == 3) ax = 0;
    return split_objs(bs, ib, ie, ax, axm);
  }
  splits.push_back(split { ax, im, in });
  int is = splits.size() - 1;
  if (++ax == 3) ax = 0;
  splits[is].left = split_objs(bs, ib, im, ax, -1);
  splits[is].center = split_objs(bs, im, in, ax, -1);
  splits[is].right = split_objs(bs, in, ie, ax, -1);
  return { merge(merge(splits[is].left.bo, splits[is].center.bo), splits[is].right.bo), is };
}

std::vector<Light::ptr> distant_lights;

void prepare_lights() {
  for (Light::ptr p: Scene::lights) {
    if (p->surrounding) distant_lights.push_back(p);
    Light::spherical const *l = dynamic_cast<Light::spherical const *>(p);
    if (!l) continue;
    Scene::objects.push_back({ l->sph, new Material::emissive(l), NULL });
  }
}

void prepare_bounds() {
  std::vector<boxed_subobject> bs;
  for (auto const &o: Scene::objects) {
    int on = &o - &Scene::objects[0];
    int n = o.solid->subparts();
    for (int i = 0; i < n; ++i) {
      box b = o.solid->bounds(i, o.transf);
      bs.push_back(boxed_subobject { on, i, b });
    }
  }
  split_objs(bs, 0, bs.size(), 0, -1);
  for (auto const &b: bs) {
    subobjects.push_back(subobject { b.obj, b.data });
  }
}

struct contact {
  double dist;
  object const *obj;
  int data;
};

contact find_range(vec const &pos, vec const &dir, int ib, int ie) {
  object const *bo = NULL;
  double bd = INFINITY;
  int bi = 0;
  for (int i = ib; i < ie; ++i) {
    subobject const &o = subobjects[i];
    object const &obj = Scene::objects[o.obj];
    ++dbg->solids;
    vec pos2 = pos, dir2 = dir;
    if (obj.transf) {
      pos2 = obj.transf->position(pos);
      dir2 = obj.transf->rotate * dir;
    }
    double d, d2 = 0.;
    for (;;) {
      d2 += 1e-6;
      double dd = obj.solid->distance(pos2 + d2 * dir2, dir2, o.data);
      if (dd <= 1e-10) continue;
      if (dd == INFINITY) { d = INFINITY; break; }
      d2 += dd;
      d = d2;
      if (obj.transf) d *= obj.transf->scale;
      if (d >= bd) break;
      if (obj.solid->valid(pos2 + d2 * dir2, o.data)) break;
    }
    if (d >= bd) continue;
    bo = &obj;
    bd = d;
    bi = o.data;
  }
  return { bd, bo, bi };
}

struct boxed_indices {
  int ib, ie;
  split::child const *sc;
};

contact find_split(vec const &pos, vec const &dir, boxed_indices const &b0, double dmax) {
  if (b0.ib == b0.ie || !Box::intersect(pos, dir, b0.sc->bo, dmax))
    return { INFINITY, NULL, 0 };
  if (b0.sc->sp < 0) return find_range(pos, dir, b0.ib, b0.ie);
  split const &s = splits[b0.sc->sp];
  boxed_indices
    b1 { b0.ib, s.before, &s.left },
    b2 { s.before, s.after, &s.center },
    b3 { s.after, b0.ie, &s.right } ;
  if (dir[s.axis] < 0) { std::swap(b1, b3); }
  contact r1 = find_split(pos, dir, b1, dmax);
  if (r1.obj) { dmax = std::min(dmax, r1.dist); }
  contact r2 = find_split(pos, dir, b2, dmax);
  if (r2.obj) { dmax = std::min(dmax, r2.dist); }
  if (!r1.obj) { r1 = find_split(pos, dir, b3, dmax); }
  if (!r1.obj) return r2;
  return r2.obj && r2.dist < r1.dist ? r2 : r1;
}

contact find_contact(vec const &pos, vec const &dir, double dmax) {
  int nb = subobjects.size();
  if (splits.empty())
    return find_range(pos, dir, 0, nb);
  split::child init { Box::whole, 0 };
  boxed_indices b { 0, nb, &init };
  return find_split(pos, dir, b, dmax);
}

double mis_weight(double x, double y) {
  return x / (x + y);
}

sampled_spectrum handle_light(sampled_wl const &wl, path_point const &pt, Light::ptr l, sampled_spectrum const &fact) {
  if (!pt.already_illuminated)
    return fact * l->get_sp(pt.pt.pos, pt.inc, wl);
  assert(Settings::shadows == Settings::Weighted);
  // Lights without a pdf have already been fully processed.
  if (!l->has_pdf) return { 0. };
  sampled_spectrum sp = l->get_sp(pt.pt.pos, pt.inc, wl);
  if (sp.zero()) return sp;
  double l_pdf = Light::lights.pdf(l, pt.pt, pt.inc);
  double w = mis_weight(pt.pdf, l_pdf);
  return w * fact * sp;
}

sampled_spectrum path(vec const &pos, vec const &dir, sampled_wl const &wl) {
  ++dbg->samples;
  sampled_spectrum color(0.), fact(1.);
  bool all_specular = true;
  std::uniform_real_distribution dis(0., 1.);
  path_point prev { { pos, vec(), vec() }, dir, 0., false };
  for (int step = 0; step < Settings::max_steps; ++step) {
    ++dbg->rays;
    contact co = find_contact(prev.pt.pos, prev.inc, INFINITY);
    if (!co.obj) {
      if (Settings::shadows != Settings::Weighted && prev.already_illuminated) break;
      for (Light::ptr l: distant_lights) {
        color += handle_light(wl, prev, l, fact);
      }
      break;
    };
    object const &obj = *co.obj;
    //if (step == 0) { color += sampled_spectrum(1.); break; }
    /*
    double attn = log(dis(*rng)) / -0.01;
    if (attn < co.dist) {
      pos += attn * dir;
      dir = half_random(dir);
      continue;
    }
    */
    path_point curr { { prev.pt.pos + co.dist * prev.inc, vec(), -prev.inc }, vec(), 0., false };
    if (obj.transf) {
      vec p = obj.transf->position(curr.pt.pos);
      vec n = obj.solid->normal(p, co.data);
      curr.pt.normal = transpose(obj.transf->rotate) * n;
    } else {
      curr.pt.normal = obj.solid->normal(curr.pt.pos, co.data);
    }
    Material::ptr mat = obj.material;
    if (mat->kind != Material::Transmitive && (curr.pt.out | curr.pt.normal) < 1e-10) break;
    if (mat->kind == Material::Emissive) {
      if (Settings::shadows != Settings::Weighted && prev.already_illuminated) break;
      Material::emissive const *me = dynamic_cast<Material::emissive const *>(mat);
      Light::ptr l = me->light;
      color += handle_light(wl, prev, l, fact);
      break;
    }
    if (!all_specular && Settings::regularize)
      mat = mat->regularize();
    auto [r, pdf] = mat->sample(curr.pt, wl);
    if (!pdf) break;
    curr.inc = r.dir;
    curr.pdf = pdf;
    if (Settings::shadows == Settings::None || r.specular == Material::Specular || !mat->has_bxdf) {
      no_illumination:
      (void)0;
    } else {
      // If the ray is semi-specular, the diffuse component of the
      // material is unrelated to the current ray, so perform a
      // partial illumination now (wrt the bxdf) and a full lighting
      // next (wrt the ray).
      if (r.specular == Material::SemiSpecular)
        curr.already_illuminated = false;
      else
        curr.already_illuminated = true;
      auto [lr, pdf] = Light::lights.sample(curr.pt);
      auto &[l, r] = lr;
      if ((r.dir | curr.pt.normal) <= 1e-10 || !pdf)
        goto no_illumination;
      sampled_spectrum sp = l->get_sp(curr.pt.pos, r.dir, wl);
      if (sp.zero()) goto no_illumination;
      ++dbg->irays;
      contact co = find_contact(curr.pt.pos, r.dir, r.dist);
      if (co.obj) {
        Material::ptr mo = co.obj->material;
        if (mo->kind != Material::Emissive) goto no_illumination;
        Material::emissive const *me = dynamic_cast<Material::emissive const *>(mo);
        if (me->light != l) goto no_illumination;
      }
      double w = 1.;
      // In case of "plain" shadowing, any partial illumination is
      // fully processed now and ignored during next iteration.
      // Same thing for lights without a pdf.
      if (Settings::shadows == Settings::Weighted && l->has_pdf) {
        assert(mat->has_pdf);
        w = mis_weight(pdf, mat->pdf(curr.pt, r.dir));
      }
      color += (w / pdf) * fact * sp * mat->bxdf(curr.pt, r.dir, wl);
    }
    fact = (1 / pdf) * fact * r.sp;
    prev = curr;
    prev.inc = r.dir;
    all_specular &= r.specular != Material::Diffuse;
    if (mat != obj.material) delete mat;
    if (step < Settings::min_steps) continue;
    double mag = *std::max_element(fact.begin(), fact.end());
    assert(mag >= 0);
    if (mag >= 1) continue;
    if (mag <= dis(*rng)) break;
    fact = (1 / mag) * fact;
  }
  return color;
}

sampled_spectrum ray(double fx, double fy, sampled_wl const &wl) {
  auto const &m = Camera::dir;
  vec dir = {
    m[0] * fx + m[1] * fy + m[2],
    m[3] * fx + m[4] * fy + m[5],
    m[6] * fx + m[7] * fy + m[8],
  };
  dir = normalize(dir);
  return path(Camera::pos, dir, wl);
}

}

struct sampler {
  std::vector<int> values;
  int index, nb;
  sampler(int n): index(0), nb(n) {
    values.reserve(nb);
    for (int i = 0; i < nb; ++i) values.push_back(i);
    std::shuffle(values.begin(), values.end(), *rng);
  }
  int next() {
    int v = values[index];
    if (++index == nb) {
      index = 0;
      std::shuffle(values.begin(), values.end(), *rng);
    }
    return v;
  }
};

struct stats {
  double s1, s2;
  int sn;
  stats(): s1(0), s2(0), sn(0) {}
  stats &operator+=(double v)
  { s1 += v; s2 += v * v; ++sn; return *this; }
  bool accurate(double tha, double thr) const {
    double m = s1 / sn;
    double sv = s2 / sn - m * m;
    return sqrt(sv / sn) <= 0.4 * std::max(tha, thr * m);
  }
};

void pixel(image &img, int x, int y) {
  int cs = sqrt(Settings::min_samples);
  if (cs * cs < Settings::min_samples) ++cs;
  sampler cells(cs * cs);
  stats st;
  Spectrum::full_spectrum sp;
  std::uniform_real_distribution dis(0., 1.);
  auto sample = [&]() {
    int i = cells.next();
    double dx = dis(*rng) + i / cs;
    double dy = dis(*rng) + i % cs;
    dx /= cs;
    dy /= cs;
    double t = dis(*rng);
    sampled_wl wl { (1 - t) * Spectrum::min_wl + t * Spectrum::max_wl };
    int wh = std::max(img.width, img.height);
    sampled_spectrum s = Solver::ray((x + dx - img.width / 2) / wh, - (y + dy - img.height / 2) / wh, wl);
    Spectrum::add(sp, wl, s);
    for (int i = 0; i < Spectrum::nb_s; ++i)
      st += Spectrum::toY(wl.lambda[i]) * s[i] / wl.pdf[i];
  };
  for (int i = 0; i < Settings::min_samples; ++i) sample();
  while (!st.accurate(0.25 * Settings::variance, Settings::variance) && st.sn < Settings::max_samples) {
    sample();
  }
  vec color = (1. / st.sn) * Spectrum::XYZtoRGB * toXYZ(sp);
  auto clamp = [](double v) -> int {
    if (v <= 0) return 0;
    else if (v >= 1) return 255;
    else return v * 255;
  };
  img.write(x, y, clamp(color[0]), clamp(color[1]), clamp(color[2]));
  //int n = st.sn * (765. / Settings::max_samples);
  //img.write(x, y, n >= 255 ? 255 : n, n >= 510 ? 255 : (n >= 256 ? n - 255 : 0), n >= 765 ? 255 : (n >= 511 ? n - 510 : 0));
}

int main() {
  Solver::prepare_lights();
  Solver::prepare_bounds();
  int w = Settings::width, h = Settings::height;
  int b = 8;
  image img(w, h);
  int bw = (w + b - 1) / b, bh = (h + b - 1) / b;
  auto block = [=, &img](int n) {
    int u = (n % bw) * b, v = (n / bw) * b;
    for (int y = v; y < v + b && y < h; ++y) {
      for (int x = u; x < u + b && x < w; ++x) {
        pixel(img, x, y);
      }
    }
  };
  spawn(block, bw * bh);
  img.save("foo.ppm");
  int64_t samples = Debug::samples.load(), rays = Debug::rays.load(), irays = Debug::irays.load(), boxes = Debug::boxes.load(), solids = Debug::solids.load();
  std::cout << "Samples: " << samples << " (" << (double)samples / (img.width * img.height) << " per pixel)\n"
    "Bounces: " << rays - samples << " (" << (double)(rays - samples) / samples << " per samples)\n"
    "Shadows: " << irays << " (" << (double)irays / samples << " per samples)\n"
    "Bound tests: " << boxes << " (" << (double)boxes / (rays + irays) << " per ray)\n"
    "Geometry tests: " << solids << " (" << (double)solids / (rays + irays) << " per ray)\n";
  return 0;
}
