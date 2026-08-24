#include <algorithm>
#include <atomic>
#include <cassert>
#include <functional>
#include <iostream>
#include <random>
#include <string>
#include <syncstream>
#include <thread>
#include <vector>

#include "base.hpp"
#include "color.hpp"
#include "debug.hpp"
#include "geometry.hpp"
#include "image.hpp"
#include "light.hpp"
#include "linalg.hpp"
#include "material.hpp"
#include "path.hpp"
#include "sampler.hpp"

void spawn(std::function<void(int)> const &f, int n) {
  std::atomic<int> idx = n - 1;
  auto worker = [&]() {
    std::mt19937_64 gen;
    Debug::debug d;
    rng = &gen;
    dbg = &d;
    for (;;) {
      int i = idx.fetch_sub(1);
      if (i < 0) break;
      std::osyncstream(std::cout) << ((n - 1 - i) * 100 / n) << "%\r" << std::flush;
      d = { 0, 0, 0, 0, 0 };
      f(i);
      Debug::merge(d);
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

using Vector::vec;

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

using Box::box;
using Ball::ball;
using Matrix::mat;
using Spectrum::sampled_wl;
using Spectrum::sampled_spectrum;

struct path_point {
  Material::interaction pt;
  vec inc;
  double pdf;
  bool already_illuminated;
};

namespace Light {

struct light_sampler {
  biased<std::pair<Light::ptr, Light::ray>> sample(Material::interaction const &from) {
    int nbl = Scene::lights.size();
    std::uniform_int_distribution dis(0, nbl - 1);
    Light::ptr l = Scene::lights[dis(*rng)];
    auto [r, pdf] = l->sample(from.pos, from.normal);
    pdf /= nbl;
    return { { l, r }, pdf };
  }

  double pdf(Light::ptr l, Material::interaction const &from, vec const &dir) {
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
  // If axis is negative, [left,right) is the range of subobjects.
  // Otherwise, left, center, and right are the indices of the children.
  int axis;
  int left, center, right;
  box bo;
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

int split_objs(std::vector<boxed_subobject> &bs, int ib, int ie, int ax, int axm) {
  if (ib == ie) return -1;
  if (ie - ib <= 1) {
    box b = ie > ib ? bs[ib].bo : Box::empty;
    splits.push_back(split { -1, ib, -1, ie, b } );
    return splits.size() - 1;
  }
  std::sort(&bs[ib], &bs[ie], bs_cmp { ax });
  int im = (ib + ie) / 2;
  double v = bs[im - 1].bo.p2[ax];
  int in = std::partition(&bs[im], &bs[ie], bs_part { v, ax }) - &bs[0];
  if (im > ib + 1) {
    im = std::max(ib + 1, (3 * im - in) / 2);
    v = bs[im - 1].bo.p2[ax];
    in = std::partition(&bs[im], &bs[in], bs_part { v, ax }) - &bs[0];
  }
  if (ax != axm && 4 * in >= ib + 3 * ie) {
    if (axm < 0) axm = ax;
    if (++ax == 3) ax = 0;
    return split_objs(bs, ib, ie, ax, axm);
  }
  //if (ax == axm && in == ie) goto no_split;
  int axn = ax;
  if (++axn == 3) axn = 0;
  int il = split_objs(bs, ib, im, axn, -1);
  int ic = split_objs(bs, im, in, axn, -1);
  int ir = split_objs(bs, in, ie, axn, -1);
  box b = Box::empty;
  if (il >= 0) b = merge(b, splits[il].bo);
  if (ic >= 0) b = merge(b, splits[ic].bo);
  if (ir >= 0) b = merge(b, splits[ir].bo);
  splits.push_back(split { ax, il, ic, ir, b });
  return splits.size() - 1;
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

#if 0
void print_splits(int si, int d) {
  std::string indent(d, ' ');
  if (si < 0) {
    std::cout << indent << "{},\n";
    return;
  }
  split const &s = splits[si];
  if (s.axis < 0) {
    std::cout << indent << "{ ";
    for (int i = s.left; i < s.right; ++i) std::cout << i << ", ";
    std::cout << "},\n";
    return;
  }
  std::cout << indent << "{ axis = " << s.axis << ",\n";
  print_splits(s.left, d + 2);
  print_splits(s.center, d + 2);
  print_splits(s.right, d + 2);
  std::cout << indent << "},\n";
}
#endif

void prepare() {
  prepare_lights();
  prepare_bounds();
  //print_splits(Solver::splits.size() - 1, 0);
}

struct contact {
  Solid::contact co;
  double dist = INFINITY;
  object const *obj = NULL;
  int data;
};

struct contact_finder {
  vec pos, dir;
  Box::checker checker;
  contact_finder(vec const &p, vec const &d)
    : pos(p), dir(d), checker(p, d) {}
  bool check_range(int ib, int ie, contact &) const;
  bool check_split(int sp, contact &) const;
  contact operator()(double dmax) const;
};

bool contact_finder::check_range(int ib, int ie, contact &bco) const {
  bool res = false;
  for (int i = ib; i < ie; ++i) {
    subobject const &o = subobjects[i];
    object const &obj = Scene::objects[o.obj];
    ++dbg->solids;
    vec pos2 = pos, dir2 = dir;
    if (obj.transf) {
      pos2 = obj.transf->to_relative(pos);
      dir2 = obj.transf->to_relative_d(dir);
    }
    Solid::contact co;
    double d, d2 = 0.;
    for (;;) {
      d2 += 1e-6;
      double dd = obj.solid->distance(pos2 + d2 * dir2, dir2, co, o.data);
      if (dd <= 1e-10) continue;
      if (dd == INFINITY) { d = INFINITY; break; }
      d2 += dd;
      d = d2;
      if (obj.transf) d *= obj.transf->scale;
      if (d >= bco.dist) break;
      co.pos = pos2 + d2 * dir2;
      if (obj.solid->complete(co, o.data)) break;
    }
    if (d >= bco.dist) continue;
    bco.co = co;
    bco.dist = d;
    bco.obj = &obj;
    bco.data = o.data;
    res = true;
  }
  return res;
}

bool contact_finder::check_split(int sp, contact &co) const {
  if (sp < 0) return false;
  split const &s = splits[sp];
  if (!checker(s.bo, co.dist)) return false;
  if (s.axis < 0)
    return check_range(s.left, s.right, co);
  int sl = s.left, sr = s.right;
  if (dir[s.axis] < 0) { std::swap(sl, sr); }
  bool bl = check_split(sl, co);
  bool br = check_split(s.center, co);
  if (bl) return true;
  return check_split(sr, co) | br;
}

contact contact_finder::operator()(double dmax) const {
  contact co;
  co.dist = dmax;
  check_split(splits.size() - 1, co);
  return co;
}

double mis_weight(double x, double y) {
  return x / (x + y);
}

sampled_spectrum handle_light(sampled_wl const &wl, path_point const &pt, Light::ptr l) {
  if (!pt.already_illuminated)
    return l->get_sp(pt.pt.pos, pt.inc, wl);
  assert(Settings::shadows == Settings::Weighted);
  // Lights without a pdf have already been fully processed.
  if (!l->has_pdf) return { 0. };
  sampled_spectrum sp = l->get_sp(pt.pt.pos, pt.inc, wl);
  if (sp.zero()) return sp;
  double l_pdf = Light::lights.pdf(l, pt.pt, pt.inc);
  double w = mis_weight(pt.pdf, l_pdf);
  return w * sp;
}

sampled_spectrum path(vec const &pos, vec const &dir, sampled_wl const &wl) {
  ++dbg->samples;
  sampled_spectrum color(0.), fact(1.);
  bool all_specular = true;
  std::uniform_real_distribution dis(0., 1.);
  path_point prev { { pos }, dir, 0., false };
  for (int step = 0; step < Settings::max_steps; ++step) {
    ++dbg->rays;
    contact co = contact_finder(prev.pt.pos, prev.inc)(INFINITY);
    if (!co.obj) {
      if (Settings::shadows != Settings::Weighted && prev.already_illuminated) break;
      for (Light::ptr l: distant_lights) {
        color += fact * handle_light(wl, prev, l);
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
    path_point curr {
      { prev.pt.pos + co.dist * prev.inc,
        obj.solid->snormal(co.co),
        -prev.inc },
      vec(), 0., false };
    if (obj.transf)
      curr.pt.normal = obj.transf->of_relative_d(curr.pt.normal);
    curr.pt.uv = obj.solid->uv(co.co);
    Material::ptr mat = obj.material;
    if (mat->kind != Material::Transmitive && (curr.pt.out | curr.pt.normal) < 1e-10) break;
    if (mat->kind == Material::Emissive) {
      if (Settings::shadows != Settings::Weighted && prev.already_illuminated) break;
      Material::emissive const *me = dynamic_cast<Material::emissive const *>(mat);
      Light::ptr l = me->light;
      color += fact * handle_light(wl, prev, l);
      break;
    }
    if (!all_specular && Settings::regularize)
      mat = mat->regularize();
    auto [r, pdf] = mat->sample(curr.pt, wl);
    if (!pdf) break;
    curr.inc = r.dir;
    curr.pdf = pdf;
    if ((curr.inc | curr.pt.normal) < 0)
      curr.pt.pos += -1e-6 * curr.pt.normal;
    else
      curr.pt.pos += 1e-6 * curr.pt.normal;
    if (r.specular == Material::Diffuse)
      r.sp *= (curr.inc | curr.pt.normal);
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
      contact co = contact_finder(curr.pt.pos, r.dir)(r.dist);
      if (co.obj && co.dist < r.dist) {
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
      sampled_spectrum isp = (r.dir | curr.pt.normal) * mat->bxdf(curr.pt, r.dir, wl);
      color += (w / pdf) * fact * sp * isp;
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

}

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

integrator::integrator(int w, int h)
  : img(w, h) {}

void integrator::pixel(int x, int y) {
  int cs = sqrt(Settings::min_samples);
  if (cs * cs < Settings::min_samples) ++cs;
  Sampler::discrete_uniform cells(cs * cs);
  stats st;
  Spectrum::full_spectrum sp;
  std::uniform_real_distribution dis(0., 1.);
  auto sample = [&]() {
    int i = cells.sample();
    double dx = dis(*rng) + i / cs;
    double dy = dis(*rng) + i % cs;
    dx /= cs;
    dy /= cs;
    double t = dis(*rng);
    sampled_wl wl { (1 - t) * Spectrum::min_wl + t * Spectrum::max_wl };
    int wh = std::max(img.width, img.height);
    auto [camp, camd] = Camera::camera->get((x + dx - img.width / 2) / wh, - (y + dy - img.height / 2) / wh);
    sampled_spectrum s = Solver::path(camp, camd, wl);
    Spectrum::add(sp, wl, s);
    for (int i = 0; i < Spectrum::nb_s; ++i)
      st += Color::toY(wl.lambda[i]) * s[i] / wl.pdf[i];
  };
  for (int i = 0; i < Settings::min_samples; ++i) sample();
  while (!st.accurate(0.25 * Settings::variance, Settings::variance) && st.sn < Settings::max_samples) {
    sample();
  }
  vec color = (1. / st.sn) * (Color::XYZtoRGB * toXYZ(sp));
  img.write(x, y, color);
  //int n = st.sn * (765. / Settings::max_samples);
  //img.write(x, y, n >= 255 ? 255 : n, n >= 510 ? 255 : (n >= 256 ? n - 255 : 0), n >= 765 ? 255 : (n >= 511 ? n - 510 : 0));
}

void integrator::tiled() {
  int w = img.width, h = img.height;
  Debug::pixels += w * h;
  int b = 8;
  int bw = (w + b - 1) / b, bh = (h + b - 1) / b;
  auto block = [=, this](int n) {
    int u = (n % bw) * b, v = (n / bw) * b;
    for (int y = v; y < v + b && y < h; ++y) {
      for (int x = u; x < u + b && x < w; ++x) {
        pixel(x, y);
      }
    }
  };
  spawn(block, bw * bh);
}
