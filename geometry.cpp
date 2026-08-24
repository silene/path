#include <algorithm>
#include <vector>

#include "base.hpp"
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

namespace Geometry {

using Box::box;

struct subobject {
  int obj, data;
};

static std::vector<subobject> subobjects;

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

static std::vector<split> splits;

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

contact_finder::contact_finder(vec const &p, vec const &d)
  : pos(p), dir(d), checker(p, d) {}

bool contact_finder::check_range(int ib, int ie, full_contact &bco) const {
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
    contact co;
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

bool contact_finder::check_split(int sp, full_contact &co) const {
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

full_contact contact_finder::operator()(double dmax) const {
  full_contact co;
  co.dist = dmax;
  check_split(splits.size() - 1, co);
  return co;
}

}
