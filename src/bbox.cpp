#include "bbox.h"

#include "GL/glew.h"

#include <algorithm>
#include <iostream>

namespace CGL {

bool BBox::intersect(const Ray& r, double& t0, double& t1) const {

  // TODO (Part 2.2):
  // Implement ray - bounding box intersection test
  // If the ray intersected the bounding box within the range given by
  // t0, t1, update t0 and t1 with the new intersection times.

  double min_x = std::min((min[0]-r.o[0])/r.d[0], (max[0]-r.o[0])/r.d[0]);
  double max_x = std::max((min[0]-r.o[0])/r.d[0], (max[0]-r.o[0])/r.d[0]);

  double min_y = std::min((min[1]-r.o[1])/r.d[1], (max[1]-r.o[1])/r.d[1]);
  double max_y = std::max((min[1]-r.o[1])/r.d[1], (max[1]-r.o[1])/r.d[1]);

  double min_z = std::min((min[2]-r.o[2])/r.d[2], (max[2]-r.o[2])/r.d[2]);
  double max_z = std::max((min[2]-r.o[2])/r.d[2], (max[2]-r.o[2])/r.d[2]);


  double start_t = std::max(std::max(min_x, min_y), min_z);
  double end_t = std::min(std::min(max_x, max_y), max_z);

  if (start_t <= end_t)
  {
    t0 = start_t;
    t1 = end_t;
    return true;
  }
  return false;
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
