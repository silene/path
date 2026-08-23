#ifndef DEBUG_HPP
#define DEBUG_HPP

#include <atomic>

namespace Debug {

extern std::atomic<int64_t> samples;
extern std::atomic<int64_t> rays;
extern std::atomic<int64_t> irays;
extern std::atomic<int64_t> boxes;
extern std::atomic<int64_t> solids;
extern int pixels;

struct debug {
  int samples;
  int rays;
  int irays;
  int boxes;
  int solids;
};

void merge(debug const &);
void dump();

}

extern thread_local Debug::debug *dbg;

#endif
