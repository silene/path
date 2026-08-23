#include "camera.hpp"
#include "light.hpp"
#include "material.hpp"
#include "mesh.hpp"
#include "path.hpp"
#include "solid.hpp"

namespace Settings {

double max_depth = 100;
skind shadows = Weighted;
int min_steps = 4, max_steps = 20;
int min_samples = 100, max_samples = 10000;
double variance = 0.08;
bool regularize = false;

}

namespace Scene {

std::vector<Light::ptr> lights {
  new Light::directional {
    new Spectrum::blackbody { 2.5, 5500 },
    Vector::normalize({ -4., 4., -1. }) },
  new Light::uniform {
    new Spectrum::blackbody { 0.6, 6500 } },
};

/***
 * Both the model and the texture for Spot comes from Keenan Crane's
 * repository:
 * https://www.cs.cmu.edu/~kmcrane/Projects/ModelRepository/
 */

std::vector<object> objects {
  { new Solid::plane { { 0., 1., 0. }, 1. },
    new Material::lambertian { new Spectrum::xyY { 0.25, 0.55, 0.4 } },
    NULL },
  { new Solid::mesh { "objs/spot.obj" },
      new Material::lambertian {
        new Spectrum::from_texture(new Image::ppm("objs/spot.ppm")) },
    new Transform::iso { { 0., -1. - 1.8 * -0.736784, -1.2 }, 1.8, { 0., 1., 0. }, -0.7 } },
};

}

namespace Camera {

Camera::ptr camera = new simple {
  vec { 0., 1.8, -11. },
  vec { 0., 0.4, 0. },
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
