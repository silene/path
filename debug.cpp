#include <iostream>

#include "debug.hpp"

namespace Debug {

std::atomic<int64_t> samples = 0;
std::atomic<int64_t> rays = 0;
std::atomic<int64_t> irays = 0;
std::atomic<int64_t> boxes = 0;
std::atomic<int64_t> solids = 0;
int pixels = 0;

void merge(debug const &d) {
  samples += d.samples;
  rays += d.rays;
  irays += d.irays;
  boxes += d.boxes;
  solids += d.solids;
}

void dump() {
  int64_t samples = Debug::samples.load(), rays = Debug::rays.load(), irays = Debug::irays.load(), boxes = Debug::boxes.load(), solids = Debug::solids.load();
  std::cout << "Samples: " << samples << " (" << (double)samples / pixels << " per pixel)\n"
    "Bounces: " << rays - samples << " (" << (double)(rays - samples) / samples << " per samples)\n"
    "Shadows: " << irays << " (" << (double)irays / samples << " per samples)\n"
    "Bound tests: " << boxes << " (" << (double)boxes / (rays + irays) << " per ray)\n"
    "Geometry tests: " << solids << " (" << (double)solids / (rays + irays) << " per ray)\n";
}

}

thread_local Debug::debug *dbg;
