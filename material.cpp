#include <cassert>
#include <complex>

#include "material.hpp"
#include "sampler.hpp"

namespace Material {

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

double metal(vec const &normal, vec const &out, std::complex<double> eta) {
  using complex = std::complex<double>;
  double cr = out | normal;
  complex s2t = (1 - cr * cr) / (eta * eta);
  complex ct = sqrt(1. - s2t);
  complex f1 = (eta * cr - ct) / (eta * cr + ct);
  complex f2 = (cr - eta * ct) / (cr + eta * ct);
  return 0.5 * (std::norm(f1) + std::norm(f2));
}

biased<ray> emissive::sample(interaction const &, sampled_wl const &wl) const {
  assert(false);
}

sampled_spectrum lambertian::bxdf(interaction const &pt, vec const &inc, sampled_wl const &wl) const {
  return M_1_PI * sp->sample(pt.uv, wl);
}

double lambertian::pdf(interaction const &pt, vec const &inc) const {
  Sampler::hemisphere_linear s(pt.normal);
  return s.pdf(inc);
}

biased<ray> lambertian::sample(interaction const &pt, sampled_wl const &wl) const {
  Sampler::hemisphere_linear s(pt.normal);
  auto [inc, pdf] = s.sample();
  return { { M_1_PI * sp->sample(pt.uv, wl), inc, Diffuse }, pdf };
}

#if 0
struct rough: base {
  Spectrum::ptr sp;
  double roughness, alpha, scaling;
  bool specular;

  rough(Spectrum::ptr s, double r)
    : base(Solid, r >= 1e-3, r >= 1e-3), sp(s), roughness(r), specular(!has_pdf)
  {
    alpha = (1 - roughness) / roughness;
    alpha = alpha * alpha * 15;
    scaling = M_1_PI * (alpha + 2) / (4 * (1 - exp(-0.5 * M_LN2 * (alpha + 2))));
  }

  sampled_spectrum bxdf(interaction const &pt, vec const &inc, sampled_wl const &wl) const {
    assert(!specular);
    double ci = inc | pt.normal;
    double cm = normalize(inc + pt.out) | pt.normal;
    assert(cm > 0);
    double s = scaling * exp(log(cm) * alpha);
    double v = s * (1 - roughness) + roughness * ci;
    return v * sp->sample(pt.uv, wl);
  }

  double pdf(interaction const &pt, vec const &inc) const {
    assert(!specular);
    Sampler::hemisphere_power samp(pt.normal, alpha);
    vec m = normalize(inc + pt.out);
    return samp.pdf(m) / (4 * (m | pt.out));
  }

  biased<ray> sample(interaction const &pt, sampled_wl const &wl) const {
    if (specular) {
      vec inc = reflect(pt.normal, pt.out);
      return { { sp->sample(pt.uv, wl), inc, Specular }, 1. };
    }
    Sampler::hemisphere_power samp(pt.normal, alpha);
    auto [m, pdf] = samp.sample();
    double c = pt.out | m;
    pdf /= 4 * c;
    if (c < 1e-6) return { { sampled_spectrum(0.), vec() }, 0. };
    vec inc = 2 * c * m - pt.out;
    double ci = inc | pt.normal;
    if (ci < 1e-6) return { { sampled_spectrum(0.), vec() }, 0. };
    double cm = m | pt.normal;
    double s = scaling * exp(log(cm) * alpha);
    double v = s * (1 - roughness) + roughness * ci;
    return { { v * sp->sample(pt.uv, wl), inc, Diffuse }, pdf };
  }

  base const *regularize() const {
    double r = roughness;
    r = 0.4 + r * (0.2 + r * 0.4);
    return new rough(sp, r);
  }
};
#endif

biased<ray> reflective::sample(interaction const &pt, sampled_wl const &wl) const {
  vec inc = reflect(pt.normal, pt.out);
  sampled_spectrum eta_ = eta->sample({}, wl), ext = extinct->sample({}, wl);
  sampled_spectrum sp;
  for (int i = 0; i < Spectrum::nb_s; ++i) {
    sp[i] = metal(pt.normal, pt.out, std::complex(eta_[i], ext[i]));
  }
  return { { Dirac * sp, inc, Specular }, Dirac };
}

biased<ray> refractive::sample(interaction const &pt, sampled_wl const &wl) const {
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

thin_refractive::thin_refractive(double e, base const *m)
  : base(Solid, m->has_bxdf, m->has_pdf), eta(e), mat(m) {
  assert(m->kind == Solid);
}

sampled_spectrum thin_refractive::bxdf(interaction const &pt, vec const &inc, sampled_wl const &wl) const {
  auto [_, f] = refract(pt.normal, pt.out, eta);
  assert(f);
  f = f / (2 - f);
  return f * mat->bxdf(pt, inc, wl);
}

double thin_refractive::pdf(interaction const &pt, vec const &inc) const {
  auto [_, f] = refract(pt.normal, pt.out, eta);
  assert(f);
  f = f / (2 - f);
  return f * mat->pdf(pt, inc);
}

biased<ray> thin_refractive::sample(interaction const &pt, sampled_wl const &wl) const {
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

}
