#include "camera.hpp"

namespace Camera {

simple::simple(vec const &p, mat const &r, double z)
  : pos(p), rot(r), zoom(1. / (3. * z))
{}

simple::simple(vec const &p, vec const &t, double z, vec const &y)
  : simple(p, Matrix::rotation_zy(t - p, y), z)
{}

simple::simple(vec const &p, vec const &t, double z)
  : simple(p, t, z, vec { 0., 1., 0. })
{}

std::pair<vec, vec> simple::get(double fx, double fy) const {
  vec d { fx * zoom, fy * zoom, 1. };
  return { pos, normalize(rot * d) };
}

simple_lens::simple_lens(vec const &p, mat const &r, double z, double l, double f)
  : simple(p, r, z), lens(l), focal(f)
{}

simple_lens::simple_lens(vec const &p, vec const &t, double z, double l, vec const &y)
  : simple(p, t, z, y), lens(l), focal(norm(t - p))
{}

simple_lens::simple_lens(vec const &p, vec const &t, double z, double l)
  : simple_lens(p, t, z, l, vec { 0., 1., 0. })
{}

std::pair<vec, vec> simple_lens::get(double fx, double fy) const {
  // Take a random point on the lens as the start of the ray,
  // and direct the ray toward the point on the focal plane
  // that would have been targeted if there was no lens.
  point2 l = Sampler::disk_uniform();
  vec s = lens * vec { l[0], l[1], 0. };
  vec d { fx * zoom, fy * zoom, 1. };
  d = focal * d - s;
  return { pos + rot * s, normalize(rot * d) };
}

}
