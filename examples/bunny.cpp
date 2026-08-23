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

/***
 * The model is from the Stanford repository of scanned models
 * and the .obj version was recovered from Alex Jacobson:
 * https://graphics.stanford.edu/data/3Dscanrep/
 * https://github.com/alecjacobson/common-3d-test-models
 *
 * The environment map is from Bernhard Vogel and was converted
 * by Benedikt Bitterli:
 * https://dativ.at/lightprobes/
 * https://benedikt-bitterli.me/resources/
 */

namespace Scene {

std::vector<Light::ptr> lights {
  /*
  new Light::multidirectional {
    new Spectrum::blackbody { 2.5, 5800 },
    Vector::normalize({ 1., 2., -1. }), 4.65e-3 },
  */
  new Light::environment { new Image::pfm("objs/envmap.pfm"), 0.8, 0.3 }
};

std::vector<object> objects {
  { new Solid::plane { { 0., 1., 0. }, 1. },
    new Material::lambertian { new Spectrum::xyY { 0.38, 0.42, 0.8 } },
    NULL },
  { new Solid::mesh { "objs/stanford-bunny.obj" },
    new Material::thin_refractive { 1.3,
      new Material::lambertian { new Spectrum::uniform { 0.9 } }, },
    new Transform::iso { { 1.5, -1. - 19. * 0.032987, -1.2 }, 19., { 0., 1., 0. }, 1. } },
};

}

namespace Camera {

Camera::ptr camera = new simple {
  vec { 0.9, 5.3, -10. },
  vec { 0.9, 0.8, -2.5 },
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
