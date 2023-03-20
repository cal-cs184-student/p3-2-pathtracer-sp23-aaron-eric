#include "bsdf.h"

#include <algorithm>
#include <iostream>
#include <utility>

#include "application/visual_debugger.h"

using std::max;
using std::min;
using std::swap;

namespace CGL {

// Mirror BSDF //

    Vector3D MirrorBSDF::f(const Vector3D wo, const Vector3D wi) {
        return Vector3D();
    }

    Vector3D MirrorBSDF::sample_f(const Vector3D wo, Vector3D* wi, double* pdf) {

        // TODO Project 3-2: Part 1
        // Implement MirrorBSDF
        *pdf = 1;
        reflect(wo, wi);
        return reflectance / abs_cos_theta(*wi);
    }

    void MirrorBSDF::render_debugger_node()
    {
        if (ImGui::TreeNode(this, "Mirror BSDF"))
        {
            DragDouble3("Reflectance", &reflectance[0], 0.005);
            ImGui::TreePop();
        }
    }

// Microfacet BSDF //

    double MicrofacetBSDF::G(const Vector3D wo, const Vector3D wi) {
        return 1.0 / (1.0 + Lambda(wi) + Lambda(wo));
    }

    double MicrofacetBSDF::D(const Vector3D h) {
        // TODO Project 3-2: Part 2
        // Compute Beckmann normal distribution function (NDF) here.
        // You will need the roughness alpha.
        return 1.0;
    }

    Vector3D MicrofacetBSDF::F(const Vector3D wi) {
        // TODO Project 3-2: Part 2
        // Compute Fresnel term for reflection on dielectric-conductor interface.
        // You will need both eta and etaK, both of which are Vector3D.

        return Vector3D();
    }

    Vector3D MicrofacetBSDF::f(const Vector3D wo, const Vector3D wi) {
        // TODO Project 3-2: Part 2
        // Implement microfacet model here.

        return Vector3D();
    }

    Vector3D MicrofacetBSDF::sample_f(const Vector3D wo, Vector3D* wi, double* pdf) {
        // TODO Project 3-2: Part 2
        // *Importance* sample Beckmann normal distribution function (NDF) here.
        // Note: You should fill in the sampled direction *wi and the corresponding *pdf,
        //       and return the sampled BRDF value.

        *wi = cosineHemisphereSampler.get_sample(pdf);
        return MicrofacetBSDF::f(wo, *wi);
    }

    void MicrofacetBSDF::render_debugger_node()
    {
        if (ImGui::TreeNode(this, "Micofacet BSDF"))
        {
            DragDouble3("eta", &eta[0], 0.005);
            DragDouble3("K", &k[0], 0.005);
            DragDouble("alpha", &alpha, 0.005);
            ImGui::TreePop();
        }
    }

// Refraction BSDF //

    Vector3D RefractionBSDF::f(const Vector3D wo, const Vector3D wi) {
        return Vector3D();
    }

    Vector3D RefractionBSDF::sample_f(const Vector3D wo, Vector3D* wi, double* pdf) {
        // TODO Project 3-2: Part 1
        // Implement RefractionBSDF
        bool hit = refract(wo, wi, ior);
        double eta = wo.z > 0 ? 1.0 / ior : ior;
        if (!hit) {
            return Vector3D();
        }
        *pdf = 1;
        return (transmittance / abs_cos_theta(*wi)) / (eta*eta);
    }

    void RefractionBSDF::render_debugger_node()
    {
        if (ImGui::TreeNode(this, "Refraction BSDF"))
        {
            DragDouble3("Transmittance", &transmittance[0], 0.005);
            DragDouble("ior", &ior, 0.005);
            ImGui::TreePop();
        }
    }

// Glass BSDF //

    Vector3D GlassBSDF::f(const Vector3D wo, const Vector3D wi) {
        return Vector3D();
    }

    Vector3D GlassBSDF::sample_f(const Vector3D wo, Vector3D* wi, double* pdf) {

        // TODO Project 3-2: Part 1
        // Compute Fresnel coefficient and either reflect or refract based on it.

        // compute Fresnel coefficient and use it as the probability of reflection
        // - Fundamentals of Computer Graphics page 305

        bool hit = refract(wo, wi, ior);
        if (!hit) {
            reflect(wo, wi);
            *pdf = 1;
            return reflectance / abs_cos_theta(*wi);
        } else {
            double R_0 = pow((1.0 - ior) / (1.0 + ior), 2);
            double R = R_0 + (1.0-R_0) * pow(1.0-abs_cos_theta(*wi), 5);
            //std::cout << 1.0-cos_theta(*wi) << std::endl;
            if (coin_flip(R)) {
                reflect(wo, wi);
                *pdf = R;
                return R * reflectance / abs_cos_theta(*wi);
            } else {
                refract(wo, wi, ior);
                *pdf = 1.0 - R;
                double eta = wo.z > 0 ? 1.0 / ior : ior;
                return (1.0-R) * (transmittance / abs_cos_theta(*wi)) / (eta * eta);
            }
        }
        return Vector3D();
    }

    void GlassBSDF::render_debugger_node()
    {
        if (ImGui::TreeNode(this, "Refraction BSDF"))
        {
            DragDouble3("Reflectance", &reflectance[0], 0.005);
            DragDouble3("Transmittance", &transmittance[0], 0.005);
            DragDouble("ior", &ior, 0.005);
            ImGui::TreePop();
        }
    }

    void BSDF::reflect(const Vector3D wo, Vector3D* wi) {

        // TODO Project 3-2: Part 1
        // Implement reflection of wo about normal (0,0,1) and store result in wi.
        *wi = {-wo.x, -wo.y, wo.z};
        return;

    }

    bool BSDF::refract(const Vector3D wo, Vector3D* wi, double ior) {

        // TODO Project 3-2: Part 1
        // Use Snell's Law to refract wo surface and store result ray in wi.
        // Return false if refraction does not occur due to total internal reflection
        // and true otherwise. When dot(wo,n) is positive, then wo corresponds to a
        // ray entering the surface through vacuum.
        double n;
        double z;
        if (wo.z > 0) {
            //entering the material
            n = 1.0/ior;
            if (1.0 - n*n * (1.0 - wo.z*wo.z) < 0) {
                return false;
            }
            wi->z = -std::sqrt(1.0 - n*n * (1.0 - wo.z*wo.z));
            z = -std::sqrt(1.0 - n*n * (1.0 - wo.z*wo.z));
        } else {
            //exiting the material
            n = ior;
            if (1.0 - n*n * (1.0 - wo.z*wo.z) < 0) {
                return false;
            }
            wi->z = std::sqrt(1.0 - n*n*(1.0 - wo.z*wo.z));
            z = std::sqrt(1.0 - n*n * (1.0 - wo.z*wo.z));
        }
        wi->x = -n * wo.x;
        wi->y = -n * wo.y;
        *wi = Vector3D(-n * wo.x, -n * wo.y, z);
        return true;

    }

} // namespace CGL
