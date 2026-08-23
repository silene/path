#include <algorithm>
#include <cassert>
#include <fstream>
#include <iostream>
#include <sstream>

#include "mesh.hpp"

namespace Solid {

Box::box mesh::bounds(int d, Transform::ptr t) const {
  assert(t);
  std::array<int, 3> const &f = facets[d];
  vec const &p0 = vertices[f[0]], &p1 = vertices[f[1]], &p2 = vertices[f[2]];
  Box::box b = Box::merge(t->of_relative(p0), t->of_relative(p1));
  return Box::merge(b, t->of_relative(p2));
}

Ball::ball mesh::sbounds(int, Transform::ptr t) const {
  assert(false);
}

bool mesh::inside(vec const &) const {
  assert(false);
}

mesh::mesh(char const *name, bool b1, bool b2)
  : inv_normal(b1), auto_normal(b2) {
  std::ifstream file(name);
  //double ymin = INFINITY;
  bool has_n = false, has_uv = false;
  for (std::string line; std::getline(file, line); ) {
    std::string::size_type n = line.find(' ');
    if (n == std::string::npos) continue;
    if (line[0] == 'v' && n == 1) {
      std::istringstream l(line);
      l.ignore(2);
      vertices.resize(vertices.size() + 1);
      vec &v = vertices.back();
      l >> v[0] >> v[1] >> v[2];
      //if (v[1] < ymin) ymin = v[1];
    } else if (line[0] == 'f' && n == 1) {
      std::istringstream l(line);
      l.ignore(2);
      facets.resize(facets.size() + 1);
      facets_n.resize(facets_n.size() + 1);
      facets_uv.resize(facets_uv.size() + 1);
      for (int i = 0; i < 3; ++i) {
        int &v = facets.back()[i], &n = facets_n.back()[i], &u = facets_uv.back()[i];
        l >> v;
        if (v > 0) --v; else v += vertices.size();
        assert(0 <= v && (unsigned)v < vertices.size());
        n = -1;
        u = -1;
        if (l.peek() != '/') continue;
        l.ignore(1);
        if (l.peek() == '/') goto read_n;
        has_uv = true;
        l >> u;
        if (u > 0) --u; else u += textures.size();
        assert(0 <= u && (unsigned)u < textures.size());
        if (l.peek() != '/') continue;
      read_n:
        has_n = true;
        l.ignore(1);
        l >> n;
        if (n > 0) --n; else n += normals.size();
        if (!(0 <= n && (unsigned)n < normals.size())) {
          std::cerr << '"' << line << "\" " << u << ' ' << n << '\n';
        }
        assert(0 <= n && (unsigned)n < normals.size());
      }
    } else if (line[0] == 'v' && line[1] == 'n' && n == 2) {
      std::istringstream l(line);
      l.ignore(3);
      normals.resize(normals.size() + 1);
      vec &v = normals.back();
      l >> v[0] >> v[1] >> v[2];
    } else if (line[0] == 'v' && line[1] == 't' && n == 2) {
      std::istringstream l(line);
      l.ignore(3);
      textures.resize(textures.size() + 1);
      point2 &v = textures.back();
      l >> v[0] >> v[1];
    } else continue;
  }
  //std::cout << ymin << '\n';
  if (!has_uv) facets_uv.clear();
  if (!has_n) {
    facets_n.clear();
    if (auto_normal) generate_normal();
  } else auto_normal = false;
}

void mesh::generate_normal() {
  assert(normals.empty());
  normals.resize(vertices.size(), vec { 0., 0., 0. });
  for (auto const &f: facets) {
    vec const &p0 = vertices[f[0]], &p1 = vertices[f[1]], &p2 = vertices[f[2]];
    vec n = cross(p1 - p0, p2 - p0);
    if (inv_normal) n = -n;
    for (int i = 0; i < 3; ++i) {
      normals[f[i]] += n;
    }
  }
  for (vec &n: normals) {
    n = normalize(n);
  }
}

double mesh::distance(vec const &pos, vec const &dir, contact &co, int data) const {
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
  co.pos = p;
  if (inv_normal) n = -n;
  co.normal = normalize(n);
  co.uv = { u, v };
  return r;
}

vec mesh::snormal(contact const &co) const {
  std::array<int, 3> idx;
  if (auto_normal) {
    idx = facets[co.data];
  } else if (!facets_n.empty()) {
    idx = facets_n[co.data];
    if (idx[0] < 0 || idx[1] < 0 || idx[2] < 0) return co.normal;
  } else
    return co.normal;
  vec const &n0 = normals[idx[0]], &n1 = normals[idx[1]],
    &n2 = normals[idx[2]];
  auto [u, v] = co.uv;
  return normalize((1 - u - v) * n0 + u * n1 + v * n2);
}

point2 mesh::uv(contact const &co) const {
  if (facets_uv.empty()) return { 0., 0. };
  std::array<int, 3> const &ft = facets_uv[co.data];
  int i0 = ft[0], i1 = ft[1], i2 = ft[2];
  if (i0 < 0 || i1 < 0 || i2 < 0) return { 0., 0. };
  point2 const &t0 = textures[i0], &t1 = textures[i1], &t2 = textures[i2];
  auto [u, v] = co.uv;
  return {
    std::clamp((1 - u - v) * t0[0] + u * t1[0] + v * t2[0], 0., 1.),
    std::clamp((1 - u - v) * t0[1] + u * t1[1] + v * t2[1], 0., 1.) };
}

}
