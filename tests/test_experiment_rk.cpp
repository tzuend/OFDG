#include "mfem.hpp"

#include "../src/experiment_rk.hpp"

#include <cmath>
#include <iostream>
#include <stdexcept>

using namespace mfem;

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

void Require(bool condition, const char *message)
{
   if (!condition) { throw std::runtime_error(message); }
}
}

int main()
{
   try
   {
      LinearEvolution evolution;
      int calls = 0;
      int final_calls = 0;
      ExperimentRungeKutta rk(
         ExperimentRungeKutta::Scheme::SSPRK3,
         [&](Vector &, real_t dt, bool final_stage)
         {
            Require(std::abs(dt - 0.1) < 1e-15,
                    "stage callback received a partial step size");
            ++calls;
            final_calls += final_stage ? 1 : 0;
            return true;
         });
      rk.Init(evolution);
      Vector state(1);
      state(0) = 1.0;
      real_t time = 0.0;
      Require(rk.TryStep(state, time, 0.1), "valid SSPRK3 step was rejected");
      Require(calls == 3 && final_calls == 1,
              "SSPRK3 did not call the postprocessor once per stage");
      Require(std::abs(time - 0.1) < 1e-15,
              "accepted step did not advance time");
      Require(std::abs(state(0) - 1.1051666666666666) < 1e-14,
              "SSPRK3 reference fingerprint changed");

      ExperimentRungeKutta rejecting(
         ExperimentRungeKutta::Scheme::ClassicalRK4,
         [](Vector &, real_t, bool) { return false; });
      rejecting.Init(evolution);
      Vector rejected_state(1);
      rejected_state(0) = 2.0;
      real_t rejected_time = 0.3;
      Require(!rejecting.TryStep(rejected_state, rejected_time, 0.1),
              "rejecting postprocessor accepted a step");
      Require(rejected_state(0) == 2.0 && rejected_time == 0.3,
              "rejected RK step changed state or time");
   }
   catch (const std::exception &error)
   {
      std::cerr << "Experiment RK tests failed: " << error.what() << '\n';
      return 1;
   }
   std::cout << "Experiment RK callback and rejection tests passed.\n";
   return 0;
}
