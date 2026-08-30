#include "mfem.hpp"

#include "../src/rk4.hpp"

#include <cmath>
#include <iostream>
#include <stdexcept>

using namespace mfem;
using namespace ofdg;

namespace
{

class LinearEvolution final : public TimeDependentOperator
{
public:
   LinearEvolution() : TimeDependentOperator(1) { }
   void Mult(const Vector &x, Vector &y) const override
   {
      y.SetSize(1);
      y(0) = x(0);
   }
};

class CountingIdentityStabilizer
{
public:
   mutable int calls = 0;
   void CompDecay(const Vector &input, Vector &output, real_t) const
   {
      ++calls;
      output = input;
   }
};

void Require(bool condition, const char *message)
{
   if (!condition) { throw std::runtime_error(message); }
}

} // namespace

int main()
{
   try
   {
      LinearEvolution evolution;
      CountingIdentityStabilizer stage_stabilizer;
      OEDG_RK4Solver<CountingIdentityStabilizer> stage_solver(stage_stabilizer);
      stage_solver.Init(evolution);

      Vector stage_solution(1);
      stage_solution(0) = 1.0;
      real_t stage_time = 0.0;
      real_t stage_dt = 0.1;
      stage_solver.Step(stage_solution, stage_time, stage_dt);

      CountingIdentityStabilizer post_step_stabilizer;
      RK4Solver standard_solver;
      standard_solver.Init(evolution);
      Vector post_step_solution(1);
      post_step_solution(0) = 1.0;
      real_t post_step_time = 0.0;
      real_t post_step_dt = 0.1;
      standard_solver.Step(post_step_solution, post_step_time, post_step_dt);
      Vector stabilized;
      post_step_stabilizer.CompDecay(post_step_solution, stabilized,
                                     post_step_dt);
      post_step_solution = stabilized;

      Require(stage_stabilizer.calls == 4,
              "stage-wise RK4 must apply decay after all four completed stages");
      Require(post_step_stabilizer.calls == 1,
              "standard RK4 path must apply decay once per full step");
      Require(std::abs(stage_solution(0) - post_step_solution(0)) < 1e-14,
              "identity stabilization changed the underlying RK4 result");
      Require(std::abs(post_step_solution(0) - 1.1051708333333332) < 1e-14,
              "RK4 reference fingerprint changed");
   }
   catch (const std::exception &error)
   {
      std::cerr << "RK4 cadence test failed: " << error.what() << '\n';
      return 1;
   }

   std::cout << "RK4 stage-versus-full-step cadence test passed.\n";
   return 0;
}
