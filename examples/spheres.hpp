#ifdef SETTINGS
int const width = 400, height = 400;
double const max_depth = 100;
skind const shadows = Weighted;
int const min_steps = 4;
int const max_steps = 20;
int const min_samples = 400;
int const max_samples = 10000;
double const variance = 0.08;
bool const regularize = false;
#endif

#ifdef SCENE
namespace Scene {

std::vector<Light::ptr> lights {
  new Light::spherical {
    new Spectrum::blackbody { 500., 5800 },
    { 5., 5., -5. }, 1. },
  new Light::uniform {
    new Spectrum::blackbody { 0.5, 4300 } },
};

std::vector<object> objects {
  { new Solid::plane { { 0., 1., 0. }, 1. },
    new Material::lambertian { new Spectrum::from_palette { 0.3, 5, 9 } },
    NULL },
  { new Solid::sphere { { 0.2, 0., -1.2 }, 1. },
    new Material::refractive { 1.2 },
    NULL },
  { new Solid::sphere { { 2.3, 0., -1. }, 1. },
    new Material::thin_refractive { 1.3,
      new Material::lambertian { new Spectrum::from_palette { 0.9, 2, 15 } }, },
    NULL },
  { new Solid::sphere { { 1.1, 0., 0.7 }, 1. },
    new Material::lambertian { new Spectrum::from_palette { 0.9, 8, 10 } },
    NULL },
};

}

namespace Camera {

vec pos { 0.8, 2., -5. };

std::array<double, 9> dir = {
  1., 0., 0.,
  0., cos(0.4), -sin(0.4),
  0., sin(0.4), cos(0.4)
};

}
#endif
