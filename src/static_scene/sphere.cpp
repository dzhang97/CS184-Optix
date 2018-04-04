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
  double alpha = dot(r.d, r.d);
  double beta = dot(2 * (r.o - o), r.d);
  double gamma = dot(r.o - o, r.o - o) - r2;
  double dual = sqrt(beta*beta - 4*alpha*gamma);
  t1 = (-beta - dual) / (2*alpha);
  t2 = (-beta + dual) / (2*alpha);
  if (t2 < t1) {
    double temp = t1;
    t1 = t2;
    t2 = temp;
  }
  return true;
}

bool Sphere::intersect(const Ray& r) const {

  // TODO (Part 1.4):
  // Implement ray - sphere intersection.
  // Note that you might want to use the the Sphere::test helper here.
  double t1, t2;
  test(r, t1, t2);
  if (t1 < r.max_t && t1 > r.min_t) {
    r.max_t = t1;
    return true;
  } else if (t2 < r.max_t && t2 > r.min_t) {
    r.max_t = t2;
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
  if (intersect(r)) {
    i->t = r.max_t;
    i->primitive = this;
    i->bsdf = get_bsdf();
    i->n = normal(r.o + r.d * r.max_t);
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
