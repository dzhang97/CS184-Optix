#include "triangle.h"

#include "CGL/CGL.h"
#include "GL/glew.h"

namespace CGL { namespace StaticScene {

Triangle::Triangle(const Mesh* mesh, size_t v1, size_t v2, size_t v3) :
    mesh(mesh), v1(v1), v2(v2), v3(v3) { }

BBox Triangle::get_bbox() const {

  Vector3D p1(mesh->positions[v1]), p2(mesh->positions[v2]), p3(mesh->positions[v3]);
  BBox bb(p1);
  bb.expand(p2); 
  bb.expand(p3);
  return bb;

}

bool Triangle::intersect(const Ray& r) const {
  
  // Part 1, Task 3: implement ray-triangle intersection
  Vector3D p0(mesh->positions[v1]), p1(mesh->positions[v2]), p2(mesh->positions[v3]);
  Vector3D e1 = p1 - p0;
  Vector3D e2 = p2 - p0;
  Vector3D s = r.o - p0;
  Vector3D s1 = cross(r.d, e2);
  Vector3D s2 = cross(s, e1);
  Vector3D out = (Vector3D(dot(s2, e2), dot(s1, s), dot(s2, r.d)) / (dot(s1, e1)));
  if (out.x > r.max_t || out.x < r.min_t || out.y < 0 || out.z < 0 || out.y + out.z > 1) return false;
  r.max_t = out.x;
  return true;
}

bool Triangle::intersect(const Ray& r, Intersection *isect) const {
  
  // Part 1, Task 3: 
  // implement ray-triangle intersection. When an intersection takes
  // place, the Intersection data should be updated accordingly
  Vector3D p0(mesh->positions[v1]), p1(mesh->positions[v2]), p2(mesh->positions[v3]);
  Vector3D n1(mesh->normals[v1]), n2(mesh->normals[v2]), n3(mesh->normals[v3]);
  Vector3D e1 = p1 - p0;
  Vector3D e2 = p2 - p0;
  Vector3D s = r.o - p0;
  Vector3D s1 = cross(r.d, e2);
  Vector3D s2 = cross(s, e1);
  Vector3D out = (Vector3D(dot(s2, e2), dot(s1, s), dot(s2, r.d)) / (dot(s1, e1)));
  if (out.x > r.max_t || out.x < r.min_t || out.y < 0 || out.z < 0 || out.y + out.z > 1) return false;
  r.max_t = out.x;
  if (intersect(r)) {
    isect->t = r.max_t;
    isect->primitive = this;
    isect->bsdf = get_bsdf();
    isect->n = (1-out.y-out.z) * n1 + out.y * n2 + out.z * n3;
    return true;
  }  
  return false;
}

void Triangle::draw(const Color& c) const {
  glColor4f(c.r, c.g, c.b, c.a);
  glBegin(GL_TRIANGLES);
  glVertex3d(mesh->positions[v1].x,
             mesh->positions[v1].y,
             mesh->positions[v1].z);
  glVertex3d(mesh->positions[v2].x,
             mesh->positions[v2].y,
             mesh->positions[v2].z);
  glVertex3d(mesh->positions[v3].x,
             mesh->positions[v3].y,
             mesh->positions[v3].z);
  glEnd();
}

void Triangle::drawOutline(const Color& c) const {
  glColor4f(c.r, c.g, c.b, c.a);
  glBegin(GL_LINE_LOOP);
  glVertex3d(mesh->positions[v1].x,
             mesh->positions[v1].y,
             mesh->positions[v1].z);
  glVertex3d(mesh->positions[v2].x,
             mesh->positions[v2].y,
             mesh->positions[v2].z);
  glVertex3d(mesh->positions[v3].x,
             mesh->positions[v3].y,
             mesh->positions[v3].z);
  glEnd();
}



} // namespace StaticScene
} // namespace CGL
