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

BVHNode *BVHAccel::construct_bvh(const std::vector<Primitive*>& prims, size_t max_leaf_size) {
  
  // TODO (Part 2.1):
  // Construct a BVH from the given vector of primitives and maximum leaf
  // size configuration. The starter code build a BVH aggregate with a
  // single leaf node (which is also the root) that encloses all the
  // primitives.
  //cout << prims.size() << endl;
  BBox bbox;

  std::vector<double> xs, ys, zs;

  for (Primitive *p : prims) {
    BBox bb = p->get_bbox();
    bbox.expand(bb);
    xs.push_back(bb.centroid().x);
    ys.push_back(bb.centroid().y);
    zs.push_back(bb.centroid().z);
  }

  if (prims.size() <= max_leaf_size) {
    BVHNode *node = new BVHNode(bbox);
    node->prims = new vector<Primitive *>(prims);
    return node;
  }

  double dx = bbox.max.x - bbox.min.x;
  double dy = bbox.max.y - bbox.min.y;
  double dz = bbox.max.z - bbox.min.z;
  double mid;
  int ind;

  if (dx >= dy && dx >= dz) {
    std::sort(xs.begin(), xs.end());
    mid = xs[xs.size() / 2];
    ind = 0;
  } else if (dy >= dx && dy >= dz) {
    std::sort(ys.begin(), ys.end());
    mid = ys[ys.size() / 2];
    ind = 1;
  } else {
    std::sort(zs.begin(), zs.end());
    mid = zs[zs.size() / 2];
    ind = 2;
  }

  vector<Primitive*> *left = new vector<Primitive *>();
  vector<Primitive*> *right = new vector<Primitive *>();

  for (Primitive *p : prims) {
    if (p->get_bbox().centroid()[ind] <= mid) {
      left->push_back(p);
    } else {
      right->push_back(p);
    }
  }

  if (left->size() == 0 || right->size() == 0) {
    vector<Primitive*> *big, *small;
    if (left->size() == 0) {
      small = left;
      big = right;
    } else {
      small = right;
      big = left;
    }
    for (int i = 0; i < big->size() / 2; i++) {
      small->push_back(big->back());
      big->pop_back();
    }
  }

  //cout << left->size() << " " << right->size() << endl;

  BVHNode *node = new BVHNode(bbox);
  node->prims = new vector<Primitive *>(prims);
  node->l = construct_bvh(*left, max_leaf_size);
  node->r = construct_bvh(*right, max_leaf_size);
  return node;
}


bool BVHAccel::intersect(const Ray& ray, BVHNode *node) const {
  // TODO (Part 2.3):
  // Fill in the intersect function.
  // Take note that this function has a short-circuit that the
  // Intersection version cannot, since it returns as soon as it finds
  // a hit, it doesn't actually have to find the closest hit.
  double t0, t1;
  if (!node->bb.intersect(ray, t0, t1)) {
    return false;
  }
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
  double t0, t1;
  if (!node->bb.intersect(ray, t0, t1)) {
    return false;
  }
  if (node->isLeaf()) {
    //cout << (*(node->prims)).size() << endl;
    bool hit = false;
    for (Primitive *p : *(node->prims)) {
      total_isects++;
      if (p->intersect(ray, i)) 
        hit = true;
    }
    return hit;
  }
  bool b0 = intersect(ray, i, node->l);
  bool b1 = intersect(ray, i, node->r);
  return b0 || b1;
}

}  // namespace StaticScene
}  // namespace CGL
