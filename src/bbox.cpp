#include "bbox.h"

#include "GL/glew.h"

#include <algorithm>
#include <iostream>

namespace CGL {

bool BBox::intersect(const Ray& r, double& t0, double& t1) const {

  // TODO (Part 2.2):
  // Implement ray - bounding box intersection test
  // If the ray intersected the bouding box within the range given by
  // t0, t1, update t0 and t1 with the new intersection times.
  double t0x, t1x, t0y, t1y, t0z, t1z, tmin, tmax;
  t0x = (min.x - r.o.x) / r.d.x;
  t1x = (max.x - r.o.x) / r.d.x;
  t0y = (min.y - r.o.y) / r.d.y;
  t1y = (max.y - r.o.y) / r.d.y;
  t0z = (min.z - r.o.z) / r.d.z;
  t1z = (max.z - r.o.z) / r.d.z;
  if (t0x > t1x) {
    double temp = t0x;
    t0x = t1x;
    t1x = temp;
  }
  if (t0y > t1y) {
    double temp = t0y;
    t0y = t1y;
    t1y = temp;
  }
  if (t0z > t1z) {
    double temp = t0z;
    t0z = t1z;
    t1z = temp;
  }
  tmin = std::max(t0x, std::max(t0y, t0z));
  tmax = std::min(t1x, std::min(t1y, t1z));
  if (tmin > tmax) {
    return false;
  }
  t0 = tmin;
  t1 = tmax;
  return true;
}

void BBox::draw(Color c) const {

  glColor4f(c.r, c.g, c.b, c.a);

	// top
	glBegin(GL_LINE_STRIP);
	glVertex3d(max.x, max.y, max.z);
  glVertex3d(max.x, max.y, min.z);
  glVertex3d(min.x, max.y, min.z);
  glVertex3d(min.x, max.y, max.z);
  glVertex3d(max.x, max.y, max.z);
	glEnd();

	// bottom
	glBegin(GL_LINE_STRIP);
  glVertex3d(min.x, min.y, min.z);
  glVertex3d(min.x, min.y, max.z);
  glVertex3d(max.x, min.y, max.z);
  glVertex3d(max.x, min.y, min.z);
  glVertex3d(min.x, min.y, min.z);
	glEnd();

	// side
	glBegin(GL_LINES);
	glVertex3d(max.x, max.y, max.z);
  glVertex3d(max.x, min.y, max.z);
	glVertex3d(max.x, max.y, min.z);
  glVertex3d(max.x, min.y, min.z);
	glVertex3d(min.x, max.y, min.z);
  glVertex3d(min.x, min.y, min.z);
	glVertex3d(min.x, max.y, max.z);
  glVertex3d(min.x, min.y, max.z);
	glEnd();

}

std::ostream& operator<<(std::ostream& os, const BBox& b) {
  return os << "BBOX(" << b.min << ", " << b.max << ")";
}

} // namespace CGL
