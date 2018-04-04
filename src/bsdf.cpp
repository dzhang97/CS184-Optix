#include "bsdf.h"

#include <iostream>
#include <algorithm>
#include <utility>

using std::min;
using std::max;
using std::swap;

namespace CGL {

void make_coord_space(Matrix3x3& o2w, const Vector3D& n) {

    Vector3D z = Vector3D(n.x, n.y, n.z);
    Vector3D h = z;
    if (fabs(h.x) <= fabs(h.y) && fabs(h.x) <= fabs(h.z)) h.x = 1.0;
    else if (fabs(h.y) <= fabs(h.x) && fabs(h.y) <= fabs(h.z)) h.y = 1.0;
    else h.z = 1.0;

    z.normalize();
    Vector3D y = cross(h, z);
    y.normalize();
    Vector3D x = cross(z, y);
    x.normalize();

    o2w[0] = x;
    o2w[1] = y;
    o2w[2] = z;
}

// Diffuse BSDF //

Spectrum DiffuseBSDF::f(const Vector3D& wo, const Vector3D& wi) {
  // TODO (Part 3.1): 
  // This function takes in both wo and wi and returns the evaluation of
  // the BSDF for those two directions.
  return reflectance / PI;
}

Spectrum DiffuseBSDF::sample_f(const Vector3D& wo, Vector3D* wi, float* pdf) {
  // TODO (Part 3.1): 
  // This function takes in only wo and provides pointers for wi and pdf,
  // which should be assigned by this function.
  // After sampling a value for wi, it returns the evaluation of the BSDF
  // at (wo, *wi).
  *wi = sampler.get_sample(pdf);
  return reflectance / PI;
}

// Mirror BSDF //

Spectrum MirrorBSDF::f(const Vector3D& wo, const Vector3D& wi) {
  return Spectrum();
}

Spectrum MirrorBSDF::sample_f(const Vector3D& wo, Vector3D* wi, float* pdf) {
  // TODO:
  // Implement MirrorBSDF
  reflect(wo, wi);
  *pdf = 1.0;
  return reflectance / abs_cos_theta(*wi);
}

// Microfacet BSDF //

double MicrofacetBSDF::G(const Vector3D& wo, const Vector3D& wi) {
    return 1.0 / (1.0 + Lambda(wi) + Lambda(wo));
}

double MicrofacetBSDF::D(const Vector3D& h) {
    // TODO: proj3-2, part 2
    // Compute Beckmann normal distribution function (NDF) here.
    // You will need the roughness alpha.
    return exp(-pow(sin_theta(h)/abs_cos_theta(h), 2.0) / pow(alpha, 2.0)) / (PI * pow(alpha,2.0) * pow(abs_cos_theta(h), 4.0));
}

Spectrum MicrofacetBSDF::F(const Vector3D& wi) {
    // TODO: proj3-2, part 2
    // Compute Fresnel term for reflection on dielectric-conductor interface.
    // You will need both eta and etaK, both of which are Spectrum.
    Spectrum Rs, Rp;
    Rs = ((eta*eta + k*k) - 2.0*eta*abs_cos_theta(wi) + pow(abs_cos_theta(wi), 2.0)) / ((eta*eta + k*k) + 2.0*eta*abs_cos_theta(wi) + pow(abs_cos_theta(wi), 2.0));
    Rp = ((eta*eta + k*k)*pow(abs_cos_theta(wi), 2.0) - 2.0*eta*abs_cos_theta(wi) + 1.0) / ((eta*eta + k*k)*pow(abs_cos_theta(wi), 2.0) + 2.0*eta*abs_cos_theta(wi) + 1.0);
    return (Rs + Rp) / 2.0;
}

Spectrum MicrofacetBSDF::f(const Vector3D& wo, const Vector3D& wi) {
    // TODO: proj3-2, part 2
    // Implement microfacet model here.
    if (wo.z < 0 || wi.z < 0) {
      return Spectrum(0,0,0);
    }
    return (F(wi) * G(wo, wi) * D((wo+wi).unit())) / (4.0 * wo.z * wi.z);
}

Spectrum MicrofacetBSDF::sample_f(const Vector3D& wo, Vector3D* wi, float* pdf) {
    // TODO: proj3-2, part 2
    // *Importance* sample Beckmann normal distribution function (NDF) here.
    // Note: You should fill in the sampled direction *wi and the corresponding *pdf,
    //       and return the sampled BRDF value.
    Vector2D uv = sampler.get_sample();
    double theta = atan(sqrt(-pow(alpha,2.0)*log(1.0-uv.x)));
    double phi = 2.0 * PI * uv.y;
    Vector3D h = Vector3D(sin(theta)*cos(phi), sin(theta)*sin(phi), cos(theta)).unit();
    *wi = (-wo + 2.0 * dot(wo, h) * h).unit();
    *pdf = exp(-pow(tan(theta),2.0)/pow(alpha,2.0)) / (4.0 * dot(*wi, h) * PI * pow(alpha,2.0) * pow(cos(theta),3.0));
    return MicrofacetBSDF::f(wo, *wi);
}

// Refraction BSDF //

Spectrum RefractionBSDF::f(const Vector3D& wo, const Vector3D& wi) {
  return Spectrum();
}

Spectrum RefractionBSDF::sample_f(const Vector3D& wo, Vector3D* wi, float* pdf) {

  // TODO:
  // Implement RefractionBSDF
  return Spectrum();
}

// Glass BSDF //

Spectrum GlassBSDF::f(const Vector3D& wo, const Vector3D& wi) {
  return Spectrum();
}

Spectrum GlassBSDF::sample_f(const Vector3D& wo, Vector3D* wi, float* pdf) {
  // TODO: 3-2 Part 1 Task 4
  // Compute Fresnel coefficient and either reflect or refract based on it.
  double n, n1, n2;
  if (cos_theta(wo) < 0) {
    n = ior;
    n1 = ior;
    n2 = 1.0;
  }
  else {
    n = 1.0 / ior;
    n1 = 1.0;
    n2 = ior;
  }

  if (!refract(wo, wi, ior)) {
    reflect(wo, wi);
    *pdf = 1.0;
    return reflectance / abs_cos_theta(*wi);
  }

  double R0 = pow((n1 - n2) / (n1 + n2), 2);
  double R = R0 + (1 - R0) * pow((1 - abs_cos_theta(wo)), 5);

  if (coin_flip(R)) {
    reflect(wo, wi);
    *pdf = R;
    return R * reflectance / abs_cos_theta(*wi);
  } else {
    *pdf = 1 - R;
    return (1 - R) * transmittance / abs_cos_theta(*wi) / (n * n);
  }
}


void BSDF::reflect(const Vector3D& wo, Vector3D* wi) {
  // TODO: 3-2 Part 1 Task 1
  // Implement reflection of wo about normal (0,0,1) and store result in wi.
  wi->x = -wo.x;
  wi->y = -wo.y;
  wi->z = wo.z;
}

bool BSDF::refract(const Vector3D& wo, Vector3D* wi, float ior) {
  // TODO: 3-2 Part 1 Task 3
  // Use Snell's Law to refract wo surface and store result ray in wi.
  // Return false if refraction does not occur due to total internal reflection
  // and true otherwise. When wo.z is positive, then wo corresponds to a
  // ray entering the surface through vacuum.
  double n;
  if (cos_theta(wo) < 0) n = ior;
  else n = 1.0 / ior;

  if (1 - n*n * (1 - cos_theta(wo)*cos_theta(wo)) < 0) return false;

  wi->x = -n * wo.x;
  wi->y = -n * wo.y;
  wi->z = sqrt(1 - n*n * (1 - wo.z*wo.z));
  if (wo.z > 0) wi->z = -wi->z;
  return true;
}

// Emission BSDF //

Spectrum EmissionBSDF::f(const Vector3D& wo, const Vector3D& wi) {
  return Spectrum();
}

Spectrum EmissionBSDF::sample_f(const Vector3D& wo, Vector3D* wi, float* pdf) {
  *pdf = 1.0 / PI;
  *wi  = sampler.get_sample(pdf);
  return Spectrum();
}

} // namespace CGL
