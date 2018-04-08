#include "bvh.h"

#include "CGL/CGL.h"
#include "static_scene/triangle.h"

#include <iostream>
#include <stack>

using namespace std;

namespace CGL { namespace StaticScene {

BVHAccel::BVHAccel(const std::vector<Primitive *> &_primitives,
                   size_t max_leaf_size) {

  root = construct_bvh(_primitives, max_leaf_size);

}

BVHAccel::~BVHAccel() {
  if (root) delete root;
}

BBox BVHAccel::get_bbox() const {
  return root->bb;
}

void BVHAccel::draw(BVHNode *node, const Color& c) const {
  if (node->isLeaf()) {
    for (Primitive *p : *(node->prims))
      p->draw(c);
  } else {
    draw(node->l, c);
    draw(node->r, c);
  }
}

void BVHAccel::drawOutline(BVHNode *node, const Color& c) const {
  if (node->isLeaf()) {
    for (Primitive *p : *(node->prims))
      p->drawOutline(c);
  } else {
    drawOutline(node->l, c);
    drawOutline(node->r, c);
  }
}

void split(BVHNode *parent, BBox centroid_bbox, size_t max_leaf_size) {
  if (parent->prims->size() <= max_leaf_size)
    return;

  vector<Primitive *> *prims_left = new vector<Primitive *>();
  vector<Primitive *> *prims_right = new vector<Primitive *>();
  BBox centroid_box_left, bbox_left, centroid_box_right, bbox_right;

  int axis = 0;
  if (centroid_bbox.extent[1] > centroid_bbox.extent[0] && centroid_bbox.extent[1] > centroid_bbox.extent[2])
    axis = 1;
  else if (centroid_bbox.extent[2] > centroid_bbox.extent[0] && centroid_bbox.extent[2] > centroid_bbox.extent[1])
    axis = 2;
  double s = centroid_bbox.min[axis] + (centroid_bbox.extent[axis] / 2);

  for (Primitive *p : *(parent->prims)) {
    BBox bb = p->get_bbox();
    Vector3D c = bb.centroid();
    if (c[axis] <= s) {
      bbox_left.expand(bb);
      centroid_box_left.expand(c);
      prims_left->push_back(p);
    } else {
      bbox_right.expand(bb);
      centroid_box_right.expand(c);
      prims_right->push_back(p);
    }
  }

  BVHNode *l = new BVHNode(bbox_left);
  l->prims = prims_left;
  BVHNode *r = new BVHNode(bbox_right);
  r->prims = prims_right;

  parent->prims = NULL;
  split(l, centroid_box_left, max_leaf_size);
  split(r, centroid_box_right, max_leaf_size);

  parent->l = l;
  parent->r = r;
  return;
}

BVHNode *BVHAccel::construct_bvh(const std::vector<Primitive*>& prims, size_t max_leaf_size) {
  
  // TODO (Part 2.1):
  // Construct a BVH from the given vector of primitives and maximum leaf
  // size configuration. The starter code build a BVH aggregate with a
  // single leaf node (which is also the root) that encloses all the
  // primitives.

  BBox centroid_box, bbox;

  for (Primitive *p : prims) {
    BBox bb = p->get_bbox();
    bbox.expand(bb);
    Vector3D c = bb.centroid();
    centroid_box.expand(c);
  }

  BVHNode *node = new BVHNode(bbox);
  node->prims = new vector<Primitive *>(prims);
  split(node, centroid_box, max_leaf_size);
  return node;
}

bool BVHAccel::intersect(const Ray& ray, BVHNode *node) const {
  // TODO (Part 2.3):
  // Fill in the intersect function.
  // Take note that this function has a short-circuit that the
  // Intersection version cannot, since it returns as soon as it finds
  // a hit, it doesn't actually have to find the closest hit.
  if (node == NULL)
    node = root;
  double t0, t1;
  if (!(node->bb.intersect(ray, t0, t1)))
    return false;
  if (node->isLeaf()) {
    for (Primitive *p : *(node->prims)) {
      total_isects++;
      if (p->intersect(ray))
        return true;
    }
    return false;
  }
  return intersect(ray, node->l) || intersect(ray, node->r);
}

bool BVHAccel::intersect(const Ray& ray, Intersection* i, BVHNode *node) const {
  // TODO (Part 2.3):
  // Fill in the intersect function.
  if (node == NULL)
    node = root;
  double t0, t1;
  if (!(node->bb.intersect(ray, t0, t1)))
    return false;
  if (node->isLeaf()) {
    bool hit = false;
    for (Primitive *p : *(node->prims)) {
      total_isects++;
      if (p->intersect(ray, i))
        hit = true;
    }
    return hit;
  }

  bool temp = intersect(ray, i, node->l);
  return intersect(ray, i, node->r) || temp;
}

}  // namespace StaticScene
}  // namespace CGL
