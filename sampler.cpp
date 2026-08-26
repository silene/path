#include <algorithm>
#include <cassert>

#include "sampler.hpp"

thread_local std::mt19937_64 *rng;

namespace Sampler {

using vec = Vector::vec;

discrete::discrete(int nb, float const *d)
  : probas(d, d + nb) {
  sum = 0.;
  for (float f: probas) sum += f;
  if (!sum) return;
  float is = 1. / sum;
  for (float &f: probas) f *= is;
  std::vector<std::pair<int, float>> under, over;
  aliases.resize(nb);
  for (int i = 0; i != nb; ++i) {
    float f = probas[i] * nb;
    if (f < 1.f) under.emplace_back(i, f);
    else over.emplace_back(i, f);
  }
  while (!under.empty()) {
    std::pair<int, float> uv = under.back();
    under.pop_back();
    if (over.empty()) {
      // should not happen, except for rounding errors
      aliases[uv.first] = std::make_pair(1.f, -1);
      continue;
    }
    std::pair<int, float> ov = over.back();
    over.pop_back();
    aliases[uv.first] = std::make_pair(uv.second, ov.first);
    ov.second -= 1.f - uv.second;
    if (ov.second < 1.f) under.push_back(ov);
    else over.push_back(ov);
  }
  for (auto ov: over) {
    aliases[ov.first] = std::make_pair(1.f, -1);
  }
}

biased<int> discrete::sample() const {
  assert(sum);
  std::uniform_real_distribution dis(0., 1.);
  double ir = dis(*rng) * probas.size();
  int i = std::min<int>(ir, probas.size() - 1);
  if (ir - i >= aliases[i].first) i = aliases[i].second;
  return { i, probas[i] };
}

discrete2D::discrete2D(int w, int h, float const *d) {
  std::vector<float> r;
  rows.reserve(h);
  r.reserve(h);
  for (int i = 0; i < h; ++i) {
    rows.emplace_back(w, d + i * w);
    r.push_back(rows.back().sum);
  }
  data = discrete(h, &r[0]);
}

double discrete2D::pdf(int x, int y) const {
    return data.pdf(y) * rows[y].pdf(x);
}

biased<std::pair<int,int>> discrete2D::sample() const {
  auto [y, py] = data.sample();
  auto [x, px] = rows[y].sample();
  return { { x, y }, px * py };
}

discrete_uniform::discrete_uniform(int n)
  : index(0), nb(n) {
  values.reserve(nb);
  for (int i = 0; i < nb; ++i) values.push_back(i);
  std::shuffle(values.begin(), values.end(), *rng);
}

int discrete_uniform::sample() {
  int v = values[index];
  if (++index == nb) {
    index = 0;
    std::shuffle(values.begin(), values.end(), *rng);
  }
  return v;
}

/*
triangle O - (1,Y) - (1,Y+dY) with Y = y/x
gets mapped to
triangle O - (cos t, sin t) - (cos (t + dt), sin (t + dt)) for some t

equality of relative areas -> dY/8 = dt/2pi
t = pi/4 Y = pi/4 y/x
*/

point2 disk_uniform() {
  std::uniform_real_distribution dis(-1., 1.);
  double x = dis(*rng), y = dis(*rng);
  if (x == 0 && y == 0) return { 0., 0. };
  double r, t;
  if (std::abs(y) > std::abs(x)) {
    r = y;
    t = 2. - x / y;
  } else {
    r = x;
    t = y / x;
  }
  t *= M_PI * 0.25;
  return { r * cos(t), r * sin(t) };
}

biased<vec> sphere_uniform::sample() const {
  std::uniform_real_distribution dis(0., 1.);
  double z = 2. * dis(*rng) - 1.;
  double r = sqrt(1. - z * z);
  double t = 2. * M_PI * dis(*rng);
  vec v { r * cos(t), r * sin(t), z };
  return { v, 0.25 * M_1_PI };
}

double sphere_uniform::pdf(vec const &) const {
  return 0.25 * M_1_PI;
}

biased<vec> hemisphere_uniform::sample() const {
  vec uu = std::abs(u[0]) < 0.5 ? vec { 1., 0., 0. } : vec { 0., 1., 0. };
  Matrix::mat rot = Matrix::rotation_zy(u, uu);
  std::uniform_real_distribution dis(0., 1.);
  double z = dis(*rng);
  double r = sqrt(1. - z * z);
  double t = 2. * M_PI * dis(*rng);
  vec v { r * cos(t), r * sin(t), z };
  return { rot * v, 0.5 * M_1_PI };
}

double hemisphere_uniform::pdf(vec const &v) const {
  if ((u | v) <= 0) return 0.;
  return 0.5 * M_1_PI;
}

/*
z = sqrt(1 - (x² + y²)) with x,y ~ disk

P(z >= Z) = P_disk(r <= sqrt(1 - Z²)) = 1 - Z²
A(z >= Z) = int_(cos t >= Z) 2pi sin t dt = 2pi (1 - Z)
p = dP/dA = dP/dZ / dA/dZ = Z/pi
*/

biased<vec> hemisphere_linear::sample() const {
  vec uu = std::abs(u[0]) < 0.5 ? vec { 1., 0., 0. } : vec { 0., 1., 0. };
  Matrix::mat rot = Matrix::rotation_zy(u, uu);
  point2 xy = disk_uniform();
  double z = sqrt(1. - xy[0]*xy[0] - xy[1]*xy[1]);
  return { rot * vec { xy[0], xy[1], z }, z * M_1_PI };
}

double hemisphere_linear::pdf(vec const &v) const {
  if ((u | v) <= 0) return 0.;
  return (u | v) * M_1_PI;
}

#if 0
biased<vec> hemisphere_power::sample() const {
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

double hemisphere_power::pdf(vec const &v) const {
  if ((u | v) <= 0) return 0.;
  return (p + 1) * exp(log(u | v) * p) * 0.5 * M_1_PI;
}
#endif

/*
z ~ U(cmax, 1)
r = sqrt(1 - z²)
x, y = r * (cos t, sin t) with t ~ U(0, 2pi)

P(z >= Z) = (1 - Z) / (1 - cmax)
A(z >= Z) = int_(cos t >= Z) 2pi sin t dt = 2pi (1 - Z)
p = dP/dA = dP/dZ / dA/dZ = 1 / (2pi (1 - cmax))
*/
biased<vec> cone_uniform::sample() const {
  vec uu = std::abs(u[0]) < 0.5 ? vec { 1., 0., 0. } : vec { 0., 1., 0. };
  Matrix::mat rot = Matrix::rotation_zy(u, uu);
  std::uniform_real_distribution dis(0., 1.);
  double z = 1. + dis(*rng) * (cmax - 1.);
  double r = sqrt(1. - z * z);
  double t = 2. * M_PI * dis(*rng);
  vec v { r * cos(t), r * sin(t), z };
  return { rot * v, 1. / (2. * M_PI * (1. - cmax)) };
}

double cone_uniform::pdf(vec const &v) const {
  if ((u | v) <= cmax) return 0.;
  return 1 / (2 * M_PI * (1 - cmax));
}

}

#if 0
int main() {
  std::mt19937_64 gen;
  rng = &gen;
  int nb = 100000000;
  double v = 0;
  vec u { 1., 0., 0. };
  Sampler::hemisphere_uniform t(u);
  Sampler::hemisphere_power s(u, 2.1);
  //Sampler::hemisphere_linear s(u);
  for (int i = 0; i < nb; ++i) {
    auto [vv, pp] = t.sample();
    v += s.pdf(vv) / pp;
  }
  std::cout << v / nb << '\n';
}
#endif
