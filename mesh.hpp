#ifndef MESH_HPP
#define MESH_HPP

#include <array>
#include <vector>

#include "solid.hpp"

namespace Solid {

struct mesh: base {
  std::vector<vec> vertices;
  std::vector<vec> normals;
  std::vector<point2> textures;
  std::vector<std::array<int, 3>> facets;
  std::vector<std::array<int, 3>> facets_n;
  std::vector<std::array<int, 3>> facets_uv;
  bool inv_normal, auto_normal;

  mesh(char const *name, bool = false, bool = true);
  int subparts() const { return facets.size(); }
  double distance(vec const &pos, vec const &dir, contact &, int data) const;
  vec snormal(contact const &) const;
  point2 uv(contact const &) const;
  void generate_normal();

  bool complete(contact &co, int d) const {
    co.data = d;
    return true;
  }

  Box::box bounds(int d, Transform::ptr t) const;
  Ball::ball sbounds(int, Transform::ptr t) const;
  bool inside(vec const &) const;
};

}

#endif
