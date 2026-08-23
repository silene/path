#include "color.hpp"
#include "light.hpp"
#include "sampler.hpp"
#include "solid.hpp"

namespace Light {

biased<ray> point::sample(vec const &p, vec const &) const {
  vec d = pos - p;
  double dist = norm(d);
  d = (1 / dist) * d;
  return { { d, dist }, Dirac };
}

sampled_spectrum point::get_sp(vec const &p, vec const &, sampled_wl const &wl) const {
  vec d = pos - p;
  double f = 1 / (d | d);
  return Dirac * f * sp->sample({}, wl);
}

biased<ray> spot::sample(vec const &p, vec const &) const {
  vec d = pos - p;
  double dist = norm(d);
  d = (1 / dist) * d;
  return { { d, dist }, Dirac };
}

sampled_spectrum spot::get_sp(vec const &p, vec const &, sampled_wl const &wl) const {
  vec d = pos - p;
  double dist = norm(d);
  d = (1 / dist) * d;
  double g = - (d | dir);
  if (g <= angle2) return sampled_spectrum(0.);
  double f = 1 / (dist * dist);
  if (g <= angle1) f = f * (g - angle2) / (angle1 - angle2);
  return Dirac * f * sp->sample({}, wl);
}

biased<ray> directional::sample(vec const &, vec const &) const {
  return { { dir, INFINITY }, Dirac };
}

sampled_spectrum directional::get_sp(vec const &, vec const &, sampled_wl const &wl) const {
  return Dirac * sp->sample({}, wl);
}

multidirectional::multidirectional(Spectrum::ptr s, vec const &d, double a)
  : base(true, true), sp(s), dir(d), cmax(cos(a)),
    inv_area(1 / (2 * M_PI * (1 - cmax)))
{}

biased<ray> multidirectional::sample(vec const &, vec const &) const {
  Sampler::cone_uniform s(dir, cmax);
  auto [dir, p] = s.sample();
  return { { dir, INFINITY }, p };
}

double multidirectional::pdf(vec const &, vec const &, vec const &d) const {
  Sampler::cone_uniform s(dir, cmax);
  return s.pdf(d);
}

sampled_spectrum multidirectional::get_sp(vec const &, vec const &d, sampled_wl const &wl) const {
  if ((dir | d) <= cmax) return { 0. };
  return inv_area * sp->sample({}, wl);
}

biased<ray> uniform::sample(vec const &, vec const &n) const {
  Sampler::hemisphere_linear s(n);
  auto [d, p] = s.sample();
  return { { d, INFINITY }, p };
}

double uniform::pdf(vec const &, vec const &n, vec const &d) const {
  Sampler::hemisphere_linear s(n);
  return s.pdf(d);
}

sampled_spectrum uniform::get_sp(vec const &, vec const &, sampled_wl const &wl) const {
  return sp->sample({}, wl);
}

environment::environment(Image::base const *i, double s, double a)
  : base(true, true), img(i), strength(s)
  , area(img->width * img->height * M_1_PI * 0.25)
  , rot(Matrix::rotation({ 0., 1., 0. }, a)) {
  if (Settings::shadows == Settings::None) return;
  std::vector<float> lum;
  lum.reserve(img->width * img->height);
  double sum = 0.;
  for (int y = 0; y < img->height; ++y) {
    for (int x = 0; x < img->width; ++x) {
      double v = Color::RGBtoY | img->read(x, y);
      sum += v;
      lum.push_back(v);
    }
  }
  if (Settings::shadows == Settings::Weighted) {
    // When the environment can be reached in two different ways,
    // the density function can ignore the darker parts.
    sum /= img->width * img->height;
    int nb = 0;
    for (float &l: lum) {
      if (l < sum) ++nb;
      l = std::max(0.f, l - (float)sum);
    }
  }
  samp = Sampler::discrete2D(img->width, img->height, &lum[0]);
}

biased<ray> environment::sample(vec const &, vec const &n) const {
  auto [xy, p] = samp.sample();
  double u = (xy.first + 0.5) / img->width, v = (xy.second + 0.5) / img->height;
  vec d = transpose(rot) * Vector::to_sphere(point2 { u, v });
  return { { d, INFINITY }, p * area };
}

double environment::pdf(vec const &, vec const &n, vec const &d) const {
  auto [u, v] = Vector::from_sphere(rot * d);
  int x = std::min<int>(u * img->width, img->width - 1);
  int y = std::min<int>(v * img->height, img->height - 1);
  return samp.pdf(x, y) * area;
}

sampled_spectrum environment::get_sp(vec const &, vec const &dir, sampled_wl const &wl) const {
  auto [u, v] = Vector::from_sphere(rot * dir);
  int x = std::min<int>(u * img->width, img->width - 1);
  int y = std::min<int>(v * img->height, img->height - 1);
  vec c = img->read(x, y);
  sampled_spectrum s = strength * fromXYZ(Color::RGBtoXYZ * c, wl);
  return s;
}

spherical::spherical(Spectrum::ptr s, vec const &c, double r)
  : base(true, false), sp(s), sph(new Solid::sphere(c, r)),
    strength(M_1_PI / (r * r))
{}

biased<ray> spherical::sample(vec const &pos, vec const &) const {
  vec v = sph->center - pos;
  double cmax = sqrt(1 - sph->radius * sph->radius / (v | v));
  v = normalize(v);
  Sampler::cone_uniform s(v, cmax);
  auto [dir,p] = s.sample();
  return { { dir, sph->distance(pos, dir) }, p };
}

double spherical::pdf(vec const &pos, vec const &, vec const &d) const {
  vec v = sph->center - pos;
  double cmax = sqrt(1 - sph->radius * sph->radius / (v | v));
  v = normalize(v);
  Sampler::cone_uniform s(v, cmax);
  return s.pdf(d);
}

sampled_spectrum spherical::get_sp(vec const &pos, vec const &dir, sampled_wl const &wl) const {
  return strength * sp->sample({}, wl);
}

}
