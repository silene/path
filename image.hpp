#ifndef IMAGE_HPP
#define IMAGE_HPP

#include <string>
#include <vector>

#include "linalg.hpp"

namespace Image {

struct base {
  int width, height;
  virtual Vector::vec read(int x, int y) const = 0;
  ~base() = default;
};

struct ppm: base {
  std::string data;
  ppm(char const *name);
  Vector::vec read(int x, int y) const;
};

struct pfm: base {
  std::vector<float> data;
  pfm(char const *name);
  Vector::vec read(int x, int y) const;
};

struct buffer {
  int width, height;
  std::vector<float> data;
  buffer(int w, int h)
    : width(w), height(h), data(w * h * 3, 0.f) {}
  void save(char const *) const;
  void write(int x, int y, Vector::vec const &c);
};

using ptr = Image::base const *;

}

#endif
