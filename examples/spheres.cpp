#include "camera.hpp"
#include "light.hpp"
#include "material.hpp"
#include "path.hpp"
#include "solid.hpp"

namespace Settings {

double max_depth = 100;
skind shadows = Weighted;
int min_steps = 6, max_steps = 20;
int min_samples = 900, max_samples = 10000;
double variance = 0.08;
bool regularize = false;

}

// Werner, Glantschnig, Ambrosch-Draxl, JPCRD 2009
Spectrum::sampled gold_n { {
  { 354.241, 1.8444 },
  { 381.490, 1.8920 },
  { 413.281, 1.8570 },
  { 450.852, 1.7933 },
  { 495.937, 1.6936 },
  { 551.041, 1.4173 },
  { 619.921, 0.8199 },
  { 708.481, 0.4391 },
  { 826.561, 0.3347 } } };

Spectrum::sampled gold_k { {
  { 354.241, 2.0676 },
  { 381.490, 2.0938 },
  { 413.281, 2.1072 },
  { 450.852, 2.1932 },
  { 495.937, 2.2562 },
  { 551.041, 2.3358 },
  { 619.921, 2.7124 },
  { 708.481, 3.6965 },
  { 826.561, 4.8406 } } };

Spectrum::sampled silver_n { {
  { 354.241, 1.4685 },
  { 381.490, 1.2298 },
  { 413.281, 0.6230 },
  { 450.852, 0.1930 },
  { 495.937, 0.1394 },
  { 551.041, 0.1338 },
  { 619.921, 0.1466 },
  { 708.481, 0.1747 },
  { 826.561, 0.2242 } } };

Spectrum::sampled silver_k { {
  { 354.241, 1.4048 },
  { 381.490, 1.3234 },
  { 413.281, 1.3282 },
  { 450.852, 2.0570 },
  { 495.937, 2.6785 },
  { 551.041, 3.2689 },
  { 619.921, 3.9130 },
  { 708.481, 4.6754 },
  { 826.561, 5.6373 } } };

namespace Scene {

std::vector<Light::ptr> lights {
  new Light::spherical {
    new Spectrum::blackbody { 12.5, 5800 },
    { 2.2, 2., -2.2 }, 0.4 },
  new Light::uniform {
    new Spectrum::blackbody { 0.6, 4300 } },
};

std::vector<object> objects {
  { new Solid::plane { { 0., 1., 0. }, 1. },
    new Material::lambertian { new Spectrum::xyY { 0.32, 0.53, 0.6 } },
    NULL },
  { new Solid::sphere { { 0.2, 0., -1.2 }, 1. },
    new Material::reflective { &gold_n, &gold_k },
    NULL },
  { new Solid::sphere { { 2.3, -0.1, -1. }, 0.9 },
    new Material::thin_refractive { 1.3,
      new Material::lambertian { new Spectrum::xyY { 0.15, 0.25, 0.9 } }, },
    NULL },
  { new Solid::sphere { { 1.1, 0., 0.7 }, 1. },
    new Material::reflective { &silver_n, &silver_k },
    NULL },
  { new Solid::sphere { { 1.3, -0.5, -2. }, 0.5 },
    new Material::refractive { 1.2 },
    NULL },
  { new Solid::sphere { { 0.1, -0.5, -2.5 }, 0.5 },
    new Material::lambertian { new Spectrum::xyY { 0.4, 0.5, 0.9 } },
    NULL },
  { new Solid::sphere { { 2.4, -0.5, -2.5 }, 0.5 },
    new Material::lambertian { new Spectrum::xyY { 0.64, 0.33, 0.9 } },
    NULL },
};

}

namespace Camera {

Camera::ptr camera = new simple {
  vec { 0.8, 5.8, -13. },
  vec { 1.2, 0., -1. },
  1.
};

}

int main() {
  Solver::prepare();
  integrator worker(400, 400);
  worker.tiled();
  worker.img.save("foo.ppm");
  Debug::dump();
  return 0;
}
