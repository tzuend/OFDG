#ifndef OEDG_RK4_SOLVER_HPP
#define OEDG_RK4_SOLVER_HPP

#include "mfem.hpp"

// #include "util.hpp"

using namespace mfem;

/**
 * Classical fourth-order explicit Runge-Kutta solver with an OE decay
 * applied after every completed RK stage.
 *
 * StabilizerType must provide:
 *
 *    void CompDecay(const Vector &input, Vector &output) const;
 *
 * If CompDecay depends on dt, call stabilizer.SetTimeStep(dt) at the
 * beginning of Step, or change CompDecay to accept dt explicitly.
 */
template <typename StabilizerType>
class OEDG_RK4Solver : public ODESolver
{
private:
   StabilizerType &stabilizer;

   Vector x0;

   Vector k1;
   Vector k2;
   Vector k3;
   Vector k4;

   Vector stage_raw;
   Vector stage_decayed;

public:
   explicit OEDG_RK4Solver(StabilizerType &stabilizer_)
      : stabilizer(stabilizer_)
   {
   }

   void Init(TimeDependentOperator &f_) override
   {
      ODESolver::Init(f_);

      const int size = f_.Height();

      x0.SetSize(size);

      k1.SetSize(size);
      k2.SetSize(size);
      k3.SetSize(size);
      k4.SetSize(size);

      stage_raw.SetSize(size);
      stage_decayed.SetSize(size);
   }

   void Step(Vector &x, real_t &t, real_t &dt) override
   {
      MFEM_VERIFY(f != nullptr,
                  "OEDG_RK4Solver must be initialized before Step().");

      // Preserve u^n throughout the entire RK step.
      x0 = x;

      // ---------------------------------------------------------------
      // k1 = L(t_n, u^n)
      // ---------------------------------------------------------------
      f->SetTime(t);
      f->Mult(x0, k1);

      // Raw first stage:
      //
      // u^(1)_h = u^n + dt/2 * k1
      add(x0, 0.5 * dt, k1, stage_raw);

      // Stabilized first stage:
      //
      // u^(1)_sigma = F_dt(u^(1)_h)
      stabilizer.CompDecay(stage_raw, stage_decayed, 0.5 * dt);

      // ---------------------------------------------------------------
      // k2 = L(t_n + dt/2, u^(1)_sigma)
      // ---------------------------------------------------------------
      f->SetTime(t + 0.5 * dt);
      f->Mult(stage_decayed, k2);

      // Raw second stage:
      //
      // u^(2)_h = u^n + dt/2 * k2
      add(x0, 0.5 * dt, k2, stage_raw);

      // Stabilized second stage.
      stabilizer.CompDecay(stage_raw, stage_decayed, 0.5 * dt);

      // ---------------------------------------------------------------
      // k3 = L(t_n + dt/2, u^(2)_sigma)
      // ---------------------------------------------------------------
      f->SetTime(t + 0.5 * dt);
      f->Mult(stage_decayed, k3);

      // Raw third stage:
      //
      // u^(3)_h = u^n + dt * k3
      add(x0, dt, k3, stage_raw);

      // Stabilized third stage.
      stabilizer.CompDecay(stage_raw, stage_decayed, dt);

      // ---------------------------------------------------------------
      // k4 = L(t_n + dt, u^(3)_sigma)
      // ---------------------------------------------------------------
      f->SetTime(t + dt);
      f->Mult(stage_decayed, k4);

      // ---------------------------------------------------------------
      // Classical RK4 output before OE decay:
      //
      // u^(n+1)_h =
      //     u^n + dt/6 * (k1 + 2 k2 + 2 k3 + k4)
      // ---------------------------------------------------------------
      stage_raw = x0;


      stage_raw.Add(dt / 6.0, k1);
      stage_raw.Add(dt / 3.0, k2);
      stage_raw.Add(dt / 3.0, k3);
      stage_raw.Add(dt / 6.0, k4);

      // Final RK stage / new timestep solution:
      //
      // u^(n+1)_sigma = F_dt(u^(n+1)_h)
      stabilizer.CompDecay(stage_raw, x, dt);

      t += dt;
      f->SetTime(t);
   }
};

#endif
