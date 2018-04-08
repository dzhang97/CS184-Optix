#include "sphere.h"

#include <cmath>

#include  "../bsdf.h"
#include "../misc/sphere_drawing.h"

namespace CGL { namespace StaticScene {

bool Sphere::test(const Ray& r, double& t1, double& t2) const {

  // TODO (Part 1.4):
  // Implement ray - sphere intersection test.
  // Return true if there are intersections and writing the
  // smaller of the two intersection times in t1 and the larger in t2.
  double a = dot(r.d, r.d);
  double b = dot(2*(r.o-o), r.d);
  double c = dot((r.o-o), (r.o-o)) - r2;
  if ((b*b - 4*a*c) >= 0)
  {
    t1 = min((-b+sqrt(b*b-4*a*c)) / 2*a, (-b-sqrt(b*b-4*a*c)) / 2*a);
    t2 = max((-b+sqrt(b*b-4*a*c)) / 2*a, (-b-sqrt(b*b-4*a*c)) / 2*a);
    if (t1 > r.max_t || t2 < r.min_t)
      return false;
    return true;
  }
  return false;
}

bool Sphere::intersect(const Ray& r) const {

  // TODO (Part 1.4):
  // Implement ray - sphere intersection.
  // Note that you might want to use the the Sphere::test helper here.
  double t1;
  double t2;
  if (test(r, t1, t2))
  {
    if (t1 < 0)
      r.max_t = t2;
    else
      r.max_t = t1;
    return true;
  }
  return false;
}

bool Sphere::intersect(const Ray& r, Intersection *i) const {

  // TODO (Part 1.4):
  // Implement ray - sphere intersection.
  // Note again that you might want to use the the Sphere::test helper here.
  // When an intersection takes place, the Intersection data should be updated
  // correspondingly.
  double t1;
  double t2;
  if (test(r, t1, t2))
  {
    if (t1 < 0)
      t1 = t2;
    i->t = t1;
    i->n = normal(t1*r.d+r.o);
    i->primitive = this;
    i->bsdf = get_bsdf();
    r.max_t = t1;
    return true;
  }
  return false;
}

void Sphere::draw(const Color& c) const {
  Misc::draw_sphere_opengl(o, r, c);
}

void Sphere::drawOutline(const Color& c) const {
    //Misc::draw_sphere_opengl(o, r, c);
}


} // namespace StaticScene
} // namespace CGL
