#ifdef SETTINGS
int const width = 400, height = 400;
double const max_depth = 100;
skind const shadows = Weighted;
int const min_steps = 4;
int const max_steps = 20;
int const min_samples = 100;
int const max_samples = 10000;
double const variance = 0.08;
bool const regularize = false;
#endif

#ifdef SCENE
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

camera const *camera = new simple {
  vec { 0., 1.8, -11. },
  vec { 0., 0.4, 0. },
  1.
};

}
#endif
