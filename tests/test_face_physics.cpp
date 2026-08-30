#include "mfem.hpp"

#include "../src/face_physics.hpp"

#include <algorithm>
#include <cmath>
#include <iostream>
#include <stdexcept>
#include <string>

using namespace mfem;

namespace
{

void Require(bool condition, const std::string &message)
{
   if (!condition) { throw std::runtime_error(message); }
}

bool Near(real_t actual, real_t expected, real_t tolerance = 1e-13)
{
   return std::abs(actual - expected) <= tolerance *
          std::max(real_t(1.0), std::max(std::abs(actual), std::abs(expected)));
}

FaceElementTransformations &InteriorFace(Mesh &mesh)
{
   for (int f = 0; f < mesh.GetNumFaces(); ++f)
   {
      FaceElementTransformations *transformations =
         mesh.GetFaceElementTransformations(f);
      if (transformations && transformations->Elem2No >= 0)
      {
         return *transformations;
      }
   }
   throw std::runtime_error("test mesh has no interior face");
}

Vector EulerState(real_t density, real_t u, real_t v, real_t pressure,
                  real_t gamma)
{
   Vector state(4);
   state(0) = density;
   state(1) = density * u;
   state(2) = density * v;
   state(3) = pressure / (gamma - 1.0) +
              0.5 * density * (u * u + v * v);
   return state;
}

void TestBurgersSpeeds(FaceElementTransformations &transformations)
{
   BurgersFacePhysics physics;
   Vector left(1), right(1), normal(2);
   left(0) = -2.0;
   right(0) = 3.0;
   normal(0) = 0.6;
   normal(1) = 0.8;

   const FacePhysicsSample sample =
      physics.Evaluate(left, right, normal, transformations);
   Require(Near(sample.beta1, 2.8), "wrong left Burgers spectral radius");
   Require(Near(sample.beta2, 4.2), "wrong right Burgers spectral radius");
   Require(Near(sample.normal_transport, 0.7),
           "wrong Burgers inflow transport speed");
}

void TestEulerSpeeds(FaceElementTransformations &transformations)
{
   constexpr real_t gamma = 1.4;
   EulerFacePhysics physics(2, gamma);
   Vector left = EulerState(1.0, 2.0, -0.5, 1.0, gamma);
   Vector right = EulerState(0.5, -1.0, 0.75, 0.4, gamma);
   Vector normal(2);
   normal(0) = 1.0;
   normal(1) = 0.0;

   const FacePhysicsSample sample =
      physics.Evaluate(left, right, normal, transformations);
   Require(Near(sample.beta1, 2.0 + std::sqrt(1.4)),
           "wrong left Euler characteristic speed");
   Require(Near(sample.beta2, 1.0 + std::sqrt(1.12)),
           "wrong right Euler characteristic speed");
   Require(Near(sample.normal_transport, 0.5),
           "Euler KXRCF transport must use the averaged fluid velocity");
}

void TestEulerStateValidation()
{
   EulerPrimitiveState primitive;
   std::string reason;

   Vector negative_density(4);
   negative_density = 0.0;
   negative_density(0) = -1.0;
   negative_density(3) = 2.5;
   Require(!DecodeEulerState(negative_density, 2, 1.4, primitive, &reason),
           "negative density was accepted");
   Require(reason.find("density") != std::string::npos,
           "density failure did not provide a useful diagnostic");

   Vector negative_pressure = EulerState(1.0, 4.0, 0.0, 1.0, 1.4);
   negative_pressure(3) = 1.0;
   Require(!DecodeEulerState(negative_pressure, 2, 1.4, primitive, &reason),
           "negative pressure was accepted");
   Require(reason.find("pressure") != std::string::npos,
           "pressure failure did not provide a useful diagnostic");
}

void TestEulerReflection()
{
   const Vector interior = EulerState(2.0, 3.0, -1.0, 1.5, 1.4);
   Vector normal(2);
   normal(0) = 0.6;
   normal(1) = 0.8;
   Vector exterior;
   ReflectEulerState(interior, normal, 2, exterior);
   Require(Near(exterior(0), interior(0)) &&
           Near(exterior(3), interior(3)),
           "reflecting wall changed density or energy");
   const real_t interior_normal =
      interior(1) * normal(0) + interior(2) * normal(1);
   const real_t exterior_normal =
      exterior(1) * normal(0) + exterior(2) * normal(1);
   const real_t interior_tangent =
      -interior(1) * normal(1) + interior(2) * normal(0);
   const real_t exterior_tangent =
      -exterior(1) * normal(1) + exterior(2) * normal(0);
   Require(Near(exterior_normal, -interior_normal),
           "reflecting wall did not reverse normal momentum");
   Require(Near(exterior_tangent, interior_tangent),
           "reflecting wall changed tangential momentum");
}

} // namespace

int main()
{
   try
   {
      Mesh mesh = Mesh::MakeCartesian2D(2, 1, Element::QUADRILATERAL,
                                        true, 1.0, 1.0);
      FaceElementTransformations &transformations = InteriorFace(mesh);
      TestBurgersSpeeds(transformations);
      TestEulerSpeeds(transformations);
      TestEulerStateValidation();
      TestEulerReflection();
   }
   catch (const std::exception &error)
   {
      std::cerr << "Face-physics tests failed: " << error.what() << '\n';
      return 1;
   }

   std::cout << "Face-physics analytical tests passed.\n";
   return 0;
}
