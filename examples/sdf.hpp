#ifdef SETTINGS
int const width = 400, height = 400;
double const max_depth = 100;
skind const shadows = Weighted;
int const min_steps = 3;
int const max_steps = 20;
int const min_samples = 100;
int const max_samples = 10000;
double const variance = 0.08;
bool const regularize = false;
#endif

#ifdef SCENE
namespace Scene {

// Werner, Glantschnig, Ambrosch-Draxl, JPCRD 2009
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

std::vector<Light::ptr> lights {
  { new Light::spot {
      new Spectrum::blackbody { 500., 5800 },
      { 10., 10., -10. },
      Vector::normalize({ -1., -1., 1. }),
      0.997, 0.99 } },
  { new Light::directional {
      new Spectrum::blackbody { 1.5, 4300 },
      Vector::normalize({ 0.2, -0.8, 0.3 }) } },
};

std::vector<object> objects {
  { new Solid::plane { { 0., 1., 0. }, 1. },
    new Material::lambertian { new Spectrum::from_palette { 0.4, 5, 9 } },
    NULL },
  { new Solid::plane { Vector::normalize ({ 1., 0., -1. }), 2.5 },
    new Material::lambertian { new Spectrum::uniform { 0.9 } },
    NULL },
  { new Solid::plane { Vector::normalize ({ -1., 0., -1. }), 2.5 },
    new Material::reflective { &silver_n, &silver_k },
    NULL },
  { new Solid::sdf {
      new SDF::smooth_union { 0.25,
        new SDF::sphere { { 2., -0.4, -1.2 }, 0.6 },
        new SDF::sphere { { 2., 0.4, -1.2 }, 0.3 } } },
    new Material::lambertian { new Spectrum::from_palette { 0.9, 7, 11 } },
    NULL },
  { new Solid::sdf {
      new SDF::difference {
        new SDF::inflated { 0.3,
          new SDF::box { { 0., 0., 0. }, { 0.7, 0.7, 0.7 } } },
        new SDF::union_ {
          new SDF::sphere { { 0., 0., -1.1 }, 0.15 },
          new SDF::sphere { { 0.3, 1.1, 0.3 }, 0.15 },
          new SDF::sphere { { -0.3, 1.1, -0.3 }, 0.15 },
          new SDF::sphere { { 1.1, -0.4, 0.4 }, 0.15 },
          new SDF::sphere { { 1.1, 0., 0. }, 0.15 },
          new SDF::sphere { { 1.1, 0.4, -0.4 }, 0.15 },
          new SDF::sphere { { -0.4, -0.2, 1.1 }, 0.15 },
          new SDF::sphere { { 0., -0.2, 1.1 }, 0.15 },
          new SDF::sphere { { 0.4, -0.2, 1.1 }, 0.15 },
          new SDF::sphere { { -0.4, 0.2, 1.1 }, 0.15 },
          new SDF::sphere { { 0., 0.2, 1.1 }, 0.15 },
          new SDF::sphere { { 0.4, 0.2, 1.1 }, 0.15 }, } }, },
    new Material::lambertian { new Spectrum::from_palette { 0.9, 10, 11 } },
    new Transform::iso { { 0., 0., 0. }, 1., { 0., 1., 0. }, -1.05 } },
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

