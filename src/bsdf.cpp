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
  *pdf = 1;
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
    double top = exp(-pow(sin_theta(h) / abs_cos_theta(h),2) / pow(alpha, 2));
    double bottom = PI * pow(alpha,2) * pow(abs_cos_theta(h), 4);
    return top / bottom;
}

Spectrum MicrofacetBSDF::F(const Vector3D& wi) {
    // TODO: proj3-2, part 2
    // Compute Fresnel term for reflection on dielectric-conductor interface.
    // You will need both eta and etaK, both of which are Spectrum.
    double costheta = abs_cos_theta(wi);
    Spectrum Rs_top = (eta*eta + k*k) - 2*eta*costheta + costheta*costheta;
    Spectrum Rs_bot = (eta*eta + k*k) + 2*eta*costheta + costheta*costheta;
    Spectrum Rp_top = (eta*eta + k*k)*costheta*costheta - 2*eta*costheta + 1;
    Spectrum Rp_bot = (eta*eta + k*k)*costheta*costheta + 2*eta*costheta + 1;
    Spectrum Rs = Rs_top / Rs_bot;
    Spectrum Rp = Rp_top / Rp_bot;
    return (Rs + Rp) / 2;
}

Spectrum MicrofacetBSDF::f(const Vector3D& wo, const Vector3D& wi) {
    Vector3D n(0,0,1);
    double dni = dot(wi, n);
    double dno = dot(wo,n);
    if (dni <= 0 || dno <= 0)
      return Spectrum(0,0,0);

    Vector3D h = (wi + wo).unit();
    return F(wi) * G(wo, wi) * D(h) / (4*dni*dno);
}

Spectrum MicrofacetBSDF::sample_f(const Vector3D& wo, Vector3D* wi, float* pdf) {
    // TODO: proj3-2, part 2
    // *Importance* sample Beckmann normal distribution function (NDF) here.
    // Note: You should fill in the sampled direction *wi and the corresponding *pdf,
    //       and return the sampled BRDF value.
    Vector2D rand2 = sampler.get_sample();
    double r1 = rand2.x;
    double r2 = rand2.y;
    double theta = atan(sqrt(-alpha*alpha*log(1-r1)));
    double phi = 2*PI*r2;

    Vector3D h(sin(theta)*cos(phi), sin(theta)*sin(phi), cos(theta));
    *wi = 2*dot(wo, h) * h - wo;

    double p_theta_top = 2 * sin(theta) * exp(-pow(tan(theta),2) / (alpha*alpha));
    double p_theta_bottom = alpha * alpha * pow(cos(theta), 3);
    double p_theta = p_theta_top / p_theta_bottom;

    double p_phi = 1/(2*PI);

    *pdf = p_theta * p_phi / (sin(theta) * 4 * dot(*wi, h));

//    *wi = cosineHemisphereSampler.get_sample(pdf); //placeholder
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
  if (!refract(wo,wi, ior)) {
    reflect(wo, wi);
    *pdf = 1;
    return reflectance / abs_cos_theta(*wi);
  }
  float R0 = pow((1.-ior)/(1.+ior), 2);
  double R = R0 + (1.-R0)*pow(1.-abs_cos_theta(wo), 5);
  if (coin_flip(R)) {
    reflect(wo, wi);
    *pdf = R;
    return R * reflectance / abs_cos_theta(*wi);
  } else {
    *pdf = 1-R;
    refract(wo,wi, ior);
    float eta;
    (wo.z > 0) ? eta = 1. / ior : eta = ior;
    return (1-R) * transmittance / abs_cos_theta(*wi) / (eta*eta);
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

  float eta;

  (wo.z > 0) ? eta = 1.0 / ior : eta = ior;

  if (1 - pow(eta,2)*(1 - pow(wo.z,2)) < 0)
    return false;
  wi->x = -eta*wo.x;
  wi->y = -eta*wo.y;
  (wo.z > 0) ? wi->z = -sqrt(1 - pow(eta,2)*(1 - pow(wo.z,2))) : wi->z = sqrt(1 - pow(eta,2)*(1 - pow(wo.z,2)));

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
