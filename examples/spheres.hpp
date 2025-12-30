#ifdef SETTINGS
int const width = 400, height = 400;
double const max_depth = 100;
skind const shadows = Weighted;
int const min_steps = 6;
int const max_steps = 20;
int const min_samples = 900;
int const max_samples = 10000;
double const variance = 0.08;
bool const regularize = false;
#endif

#ifdef SCENE
namespace Scene {

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

std::vector<Light::ptr> lights {
  new Light::spherical {
    new Spectrum::blackbody { 50., 5800 },
    { 2.2, 2., -2.2 }, 0.4 },
  new Light::uniform {
    new Spectrum::blackbody { 0.6, 4300 } },
};

std::vector<object> objects {
  { new Solid::plane { { 0., 1., 0. }, 1. },
    new Material::lambertian { new Spectrum::from_palette { 0.4, 5, 9 } },
    NULL },
  { new Solid::sphere { { 0.2, 0., -1.2 }, 1. },
    new Material::reflective { &gold_n, &gold_k },
    NULL },
  { new Solid::sphere { { 2.3, -0.1, -1. }, 0.9 },
    new Material::thin_refractive { 1.3,
      new Material::lambertian { new Spectrum::from_palette { 0.9, 2, 15 } }, },
    NULL },
  { new Solid::difference {
      new Solid::sphere { { 1.1, 0., 0.7 }, 1. },
      new Solid::union_ {
        new Solid::sphere { { 0.9, 0.5, -0.2 }, 0.2 },
        new Solid::sphere { { 1.3, 0.5, -0.2 }, 0.2 } } },
    new Material::lambertian { new Spectrum::from_palette { 0.9, 8, 10 } },
    NULL },
  { new Solid::sphere { { 1.3, -0.5, -2. }, 0.5 },
    new Material::refractive { 1.2 },
    NULL },
};

}

namespace Camera {

vec pos { 1., 2., -5. };

std::array<double, 9> dir = {
  1., 0., 0.,
  0., cos(0.45), -sin(0.45),
  0., sin(0.45), cos(0.45)
};

}
#endif
