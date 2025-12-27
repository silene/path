# A Toy Pathtracer

This project is born from my curiosity about what pathtracing is and
how it differs from raytracing (or rather what I knew of raytracing 30
years ago) but also from my intent to experiment with Monte-Carlo
methods and to hone my rusted skills when it comes to calculus (damn
density functions).

It has been heavily influenced by the book [Physically Based
Rendering: From Theory To Implementation](https://pbr-book.org/) by
Matt Pharr, Wenzel Jakob, and Greg Humphreys, which I heartily
recommend. Despite its seemingly narrow subject, this book covers a
wide variety of topics, from physic models to algorithmic concerns,
from mathematics to implementation details.

In particular, without this book, I would never have fully understood
the following points:

- Handling specular rays does not have to be a hack. They neatly fit
  into a Monte-Carlo framework, despite their Dirac distributions, as
  long as importance sampling is used.

- Sampling both the material bsdf and the light sources along a path
  is a mathematically-grounded way to perform illumination, thanks to
  multiple importance sampling.

- Sampling the whole visible spectrum rather than just three color
  channels does not make things much harder or slower. Moreover, it
  makes sense to process multiple, equally spaced wavelengths along a
  single path.

Note: A slightly modern C++ compiler is needed, as the code relies on
several C++20 features.

Disclaimer: Most samplers are naive, rejection-based ones, and should
not be taken as example. Also, the monolithic structure of the code is
not a good practice either.
