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
  new Light::multidirectional {
    new Spectrum::blackbody { 2.5, 5800 },
    Vector::normalize({ 1., 2., -1. }), 4.65e-3 },
};

std::vector<object> objects {
  { new Solid::plane { { 0., 1., 0. }, 1. },
    new Material::lambertian { new Spectrum::from_palette { 0.9, 7, 13 } },
    NULL },
  { new Solid::mesh { "objs/stanford-bunny.obj" },
      new Material::lambertian { new Spectrum::uniform { 0.9 } },
    new Transform::iso { { 1.5, -1. - 19. * 0.032987, -1.2 }, 19., { 0., 1., 0. }, 1. } },
};

}

namespace Camera {

vec pos { 0.9, 2., -5. };

std::array<double, 9> dir = {
  1., 0., 0.,
  0., cos(0.45), -sin(0.45),
  0., sin(0.45), cos(0.45)
};

}
#endif
