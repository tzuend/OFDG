#pragma once

#include <stdexcept>
#include <string>

namespace ofdg
{

enum class StudyMethod {
    DG,
    OFDG,
    OFDGKXRCF,
    OEDG2024
};

enum class FilterCadence {
    Auto,
    Step,
    Stage
};

inline StudyMethod ParseStudyMethod(const std::string &name)
{
    if (name == "dg") { return StudyMethod::DG; }
    if (name == "ofdg") { return StudyMethod::OFDG; }
    if (name == "ofdg-kxrcf") { return StudyMethod::OFDGKXRCF; }
    if (name == "oedg") { return StudyMethod::OEDG2024; }
    throw std::invalid_argument(
        "method must be dg, ofdg, ofdg-kxrcf, or oedg");
}

inline FilterCadence ParseFilterCadence(const std::string &name)
{
    if (name == "auto") { return FilterCadence::Auto; }
    if (name == "step") { return FilterCadence::Step; }
    if (name == "stage") { return FilterCadence::Stage; }
    throw std::invalid_argument("filter cadence must be auto, step, or stage");
}

inline FilterCadence ResolveFilterCadence(StudyMethod method,
                                          FilterCadence cadence)
{
    if (cadence != FilterCadence::Auto) { return cadence; }
    return method == StudyMethod::OEDG2024 ? FilterCadence::Stage
                                           : FilterCadence::Step;
}

inline const char *StudyMethodName(StudyMethod method)
{
    switch (method) {
        case StudyMethod::DG: return "dg";
        case StudyMethod::OFDG: return "ofdg";
        case StudyMethod::OFDGKXRCF: return "ofdg-kxrcf";
        case StudyMethod::OEDG2024: return "oedg";
    }
    return "unknown";
}

inline const char *FilterCadenceName(FilterCadence cadence)
{
    switch (cadence) {
        case FilterCadence::Auto: return "auto";
        case FilterCadence::Step: return "step";
        case FilterCadence::Stage: return "stage";
    }
    return "unknown";
}

} // namespace ofdg
