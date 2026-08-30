#pragma once

#include "mfem.hpp"

#include <functional>

/**
 * Small explicit RK driver used only by controlled numerical experiments.
 *
 * The callback is invoked after every completed RK stage.  It receives the
 * candidate state, the full step size used by the OE pseudo-time problem, and
 * whether the candidate is the final stage.  Returning false rejects the step
 * without changing the input state or time, allowing the caller to reduce dt.
 */
class ExperimentRungeKutta {
public:
    using StageCallback =
        std::function<bool(mfem::Vector &, mfem::real_t, bool)>;

    enum class Scheme {
        SSPRK3,
        ClassicalRK4
    };

private:
    Scheme scheme;
    StageCallback callback;
    mfem::TimeDependentOperator *operator_ = nullptr;
    mfem::Vector initial;
    mfem::Vector stage;
    mfem::Vector derivative;
    mfem::Vector k1, k2, k3, k4;

    bool Postprocess(mfem::Vector &candidate, mfem::real_t dt,
                     bool final_stage) const
    {
        return !callback || callback(candidate, dt, final_stage);
    }

    bool StepSSPRK3(mfem::Vector &state, mfem::real_t time,
                    mfem::real_t dt)
    {
        initial = state;

        operator_->SetTime(time);
        operator_->Mult(initial, derivative);
        add(initial, dt, derivative, stage);
        if (!Postprocess(stage, dt, false)) { return false; }

        operator_->SetTime(time + dt);
        operator_->Mult(stage, derivative);
        stage.Add(dt, derivative);
        stage *= 0.25;
        stage.Add(0.75, initial);
        if (!Postprocess(stage, dt, false)) { return false; }

        operator_->SetTime(time + 0.5 * dt);
        operator_->Mult(stage, derivative);
        stage.Add(dt, derivative);
        stage *= 2.0 / 3.0;
        stage.Add(1.0 / 3.0, initial);
        if (!Postprocess(stage, dt, true)) { return false; }

        state = stage;
        return true;
    }

    bool StepRK4(mfem::Vector &state, mfem::real_t time,
                 mfem::real_t dt)
    {
        initial = state;
        operator_->SetTime(time);
        operator_->Mult(initial, k1);

        add(initial, 0.5 * dt, k1, stage);
        if (!Postprocess(stage, dt, false)) { return false; }
        operator_->SetTime(time + 0.5 * dt);
        operator_->Mult(stage, k2);

        add(initial, 0.5 * dt, k2, stage);
        if (!Postprocess(stage, dt, false)) { return false; }
        operator_->SetTime(time + 0.5 * dt);
        operator_->Mult(stage, k3);

        add(initial, dt, k3, stage);
        if (!Postprocess(stage, dt, false)) { return false; }
        operator_->SetTime(time + dt);
        operator_->Mult(stage, k4);

        stage = initial;
        stage.Add(dt / 6.0, k1);
        stage.Add(dt / 3.0, k2);
        stage.Add(dt / 3.0, k3);
        stage.Add(dt / 6.0, k4);
        if (!Postprocess(stage, dt, true)) { return false; }

        state = stage;
        return true;
    }

public:
    ExperimentRungeKutta(Scheme scheme_, StageCallback callback_ = {})
        : scheme(scheme_), callback(std::move(callback_))
    {
    }

    void Init(mfem::TimeDependentOperator &operator_value)
    {
        operator_ = &operator_value;
        const int size = operator_value.Height();
        initial.SetSize(size);
        stage.SetSize(size);
        derivative.SetSize(size);
        k1.SetSize(size);
        k2.SetSize(size);
        k3.SetSize(size);
        k4.SetSize(size);
    }

    bool TryStep(mfem::Vector &state, mfem::real_t &time,
                 mfem::real_t dt)
    {
        MFEM_VERIFY(operator_ != nullptr,
                    "ExperimentRungeKutta must be initialized before use.");
        const bool accepted = scheme == Scheme::SSPRK3
                                  ? StepSSPRK3(state, time, dt)
                                  : StepRK4(state, time, dt);
        if (accepted) {
            time += dt;
            operator_->SetTime(time);
        }
        return accepted;
    }
};

