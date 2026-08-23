#include <cassert>
#include <fstream>
#include <string>

#include "image.hpp"

namespace Image {

using Vector::vec;

ppm::ppm(char const *name) {
  std::ifstream file(name);
  char header[2];
  file.read(header, 2);
  assert(header[0] == 'P' && header[1] == '6');
  int m;
  file >> width >> height >> m;
  file.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
  int s = width * height * 3;
  data.resize(s);
  if (m == 255) {
    file.read(&data[0], s);
  } else {
    for (int n = 0; n < s; ++n) {
      unsigned char u = file.get();
      data[n] = (int)(255. * u / m);
    }
  }
}

vec ppm::read(int x, int y) const {
  assert(0 <= x && x < width && 0 <= y && y < height);
  int i = (y * width + x) * 3;
  vec v;
  for (int j = 0; j < 3; ++j) {
    unsigned char c = data[i + j];
    v[j] = c * (1. / 255.);
  }
  return v;
}

pfm::pfm(char const *name) {
  std::ifstream file(name);
  char header[2];
  file.read(header, 2);
  assert(header[0] == 'P' && header[1] == 'F');
  double m;
  file >> width >> height >> m;
  assert(m < 0.);
  file.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
  int s = width * height * 3;
  data.resize(s);
  file.read((char *)&data[0], s * 12);
}

vec pfm::read(int x, int y) const {
  assert(0 <= x && x < width && 0 <= y && y < height);
  int i = ((height - 1 - y) * width + x) * 3;
  vec v;
  for (int j = 0; j < 3; ++j) v[j] = data[i + j];
  return v;
}

void buffer::write(int x, int y, vec const &c) {
  assert(0 <= x && x < width && 0 <= y && y < height);
  int i = (y * width + x) * 3;
  data[i + 0] = c[0];
  data[i + 1] = c[1];
  data[i + 2] = c[2];
}

void buffer::save(char const *name) const {
  std::string name_(name);
  assert(name_.size() >= 5);
  name_ = name_.substr(name_.size() - 4);
  std::ofstream out(name, std::ios::binary);
  if (name_ == ".ppm") {
    out << "P6 " << width << ' ' << height << " 255\n";
    int s = width * height * 3;
    std::string d(s, '\0');
    for (int i = 0; i < s; ++i) {
      float v = data[i];
      unsigned char c;
      if (v <= 0.) c = 0;
      else if (v >= 1.) c = 255;
      else c = data[i] * 255.;
      d[i] = c;
    }
    out << d << '\n';
  } else if (name_ == ".pfm") {
    out << "PF\n" << width << ' ' << height << "\n-1.0\n";
    for (int y = height - 1; y >= 0; --y) {
      out.write((char const *)&data[y * width * 3], width * 12);
    }
  }
}

}
