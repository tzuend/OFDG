// Independent one-dimensional finite-volume reference solver.
// Component-wise WENO5 reconstruction, Lax--Friedrichs flux splitting, and
// SSPRK3 time integration intentionally share no DG/OFDG/OEDG code.

#include <algorithm>
#include <array>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <stdexcept>
#include <string>
#include <vector>

namespace
{

using State = std::array<double, 3>;
constexpr double gamma_value = 1.4;

double Pressure(const State &u)
{
   return (gamma_value - 1.0) *
          (u[2] - 0.5 * u[1] * u[1] / u[0]);
}

State Flux(const State &u)
{
   const double velocity = u[1] / u[0];
   const double pressure = Pressure(u);
   return {u[1], u[1] * velocity + pressure,
           velocity * (u[2] + pressure)};
}

double Speed(const State &u)
{
   return std::abs(u[1] / u[0]) +
          std::sqrt(gamma_value * Pressure(u) / u[0]);
}

State Primitive(double density, double velocity, double pressure)
{
   return {density, density * velocity,
           pressure / (gamma_value - 1.0) +
              0.5 * density * velocity * velocity};
}

double WenoLeft(double a, double b, double c, double d, double e)
{
   const double q0 = (2.0 * a - 7.0 * b + 11.0 * c) / 6.0;
   const double q1 = (-b + 5.0 * c + 2.0 * d) / 6.0;
   const double q2 = (2.0 * c + 5.0 * d - e) / 6.0;
   const double beta0 = 13.0 / 12.0 * std::pow(a - 2.0 * b + c, 2) +
                        0.25 * std::pow(a - 4.0 * b + 3.0 * c, 2);
   const double beta1 = 13.0 / 12.0 * std::pow(b - 2.0 * c + d, 2) +
                        0.25 * std::pow(b - d, 2);
   const double beta2 = 13.0 / 12.0 * std::pow(c - 2.0 * d + e, 2) +
                        0.25 * std::pow(3.0 * c - 4.0 * d + e, 2);
   const double epsilon = 1e-12;
   const double alpha0 = 0.1 / std::pow(epsilon + beta0, 2);
   const double alpha1 = 0.6 / std::pow(epsilon + beta1, 2);
   const double alpha2 = 0.3 / std::pow(epsilon + beta2, 2);
   const double sum = alpha0 + alpha1 + alpha2;
   return (alpha0 * q0 + alpha1 * q1 + alpha2 * q2) / sum;
}

class ReferenceSolver
{
private:
   int problem;
   int cells;
   int ghost = 3;
   double scale;
   double left;
   double right;
   double dx;
   std::vector<State> state;
   std::vector<State> work;
   std::vector<State> residual;

   State Initial(double x) const
   {
      if (problem == 4)
      {
         State value = Primitive(2.0 + 2.0 * std::pow(std::sin(x), 2),
                                 1.0, 2.0);
         for (double &component : value) { component *= scale; }
         return value;
      }
      if (problem == 5)
      {
         const double pressure = x < 0.1 ? 1000.0 : (x < 0.9 ? 0.01 : 100.0);
         return Primitive(1.0, 0.0, pressure);
      }
      if (problem == 6)
      {
         return x < 0.0 ? Primitive(1.0, 0.0, 1.0)
                        : Primitive(0.125, 0.0, 0.1);
      }
      if (problem == 7)
      {
         const State value = x < 0.0 ? Primitive(0.445, 0.698, 3.528)
                                     : Primitive(0.5, 0.0, 0.571);
         State scaled = value;
         for (double &component : scaled) { component *= scale; }
         return scaled;
      }
      if (problem == 8)
      {
         return x < -4.0
                   ? Primitive(3.857143, 2.629369, 10.33333)
                   : Primitive(1.0 + 0.2 * std::sin(5.0 * x), 0.0, 1.0);
      }
      throw std::invalid_argument(
         "reference problem must be 4, 5, 6, 7, or 8");
   }

   void FillGhosts(std::vector<State> &values) const
   {
      for (int layer = 0; layer < ghost; ++layer)
      {
         if (problem == 4)
         {
            values[ghost - 1 - layer] =
               values[ghost + cells - 1 - layer];
            values[ghost + cells + layer] = values[ghost + layer];
         }
         else if (problem == 5)
         {
            values[ghost - 1 - layer] = values[ghost + layer];
            values[ghost - 1 - layer][1] *= -1.0;
            values[ghost + cells + layer] =
               values[ghost + cells - 1 - layer];
            values[ghost + cells + layer][1] *= -1.0;
         }
         else
         {
            values[ghost - 1 - layer] = values[ghost];
            values[ghost + cells + layer] = values[ghost + cells - 1];
         }
      }
   }

   void Residual(std::vector<State> &values, std::vector<State> &rhs)
   {
      FillGhosts(values);
      double alpha = 0.0;
      for (int i = ghost; i < ghost + cells; ++i)
      {
         const double pressure = Pressure(values[i]);
         if (!(values[i][0] > 0.0) || !(pressure > 0.0))
         {
            throw std::runtime_error("reference solver encountered a nonphysical state");
         }
         alpha = std::max(alpha, Speed(values[i]));
      }

      std::vector<State> plus(values.size());
      std::vector<State> minus(values.size());
      for (int i = 0; i < static_cast<int>(values.size()); ++i)
      {
         const State flux = Flux(values[i]);
         for (int c = 0; c < 3; ++c)
         {
            plus[i][c] = 0.5 * (flux[c] + alpha * values[i][c]);
            minus[i][c] = 0.5 * (flux[c] - alpha * values[i][c]);
         }
      }

      std::vector<State> interface_flux(cells + 1);
      for (int face = 0; face <= cells; ++face)
      {
         const int i = ghost - 1 + face;
         for (int c = 0; c < 3; ++c)
         {
            const double left_flux = WenoLeft(
               plus[i - 2][c], plus[i - 1][c], plus[i][c],
               plus[i + 1][c], plus[i + 2][c]);
            const double right_flux = WenoLeft(
               minus[i + 3][c], minus[i + 2][c], minus[i + 1][c],
               minus[i][c], minus[i - 1][c]);
            interface_flux[face][c] = left_flux + right_flux;
         }
      }

      for (int cell = 0; cell < cells; ++cell)
      {
         for (int c = 0; c < 3; ++c)
         {
            rhs[ghost + cell][c] =
               -(interface_flux[cell + 1][c] -
                 interface_flux[cell][c]) / dx;
         }
      }
   }

public:
   ReferenceSolver(int problem_, int cells_, double scale_)
      : problem(problem_), cells(cells_), scale(scale_),
        left(problem == 4 ? 0.0 :
             (problem == 5 ? 0.0 : (problem == 6 ? -0.5 : -5.0))),
        right(problem == 4 ? 2.0 * std::acos(-1.0) :
              (problem == 5 ? 1.0 : (problem == 6 ? 0.5 : 5.0))),
        dx((right - left) / cells), state(cells + 2 * ghost),
        work(state.size()), residual(state.size())
   {
      for (int i = 0; i < cells; ++i)
      {
         state[ghost + i] = Initial(left + (i + 0.5) * dx);
      }
      FillGhosts(state);
   }

   int Run(double final_time, double cfl)
   {
      double time = 0.0;
      int steps = 0;
      std::vector<State> first(state.size()), second(state.size());
      while (time < final_time)
      {
         double maximum_speed = 0.0;
         for (int i = ghost; i < ghost + cells; ++i)
         {
            maximum_speed = std::max(maximum_speed, Speed(state[i]));
         }
         const double dt = std::min(cfl * dx / maximum_speed,
                                    final_time - time);

         Residual(state, residual);
         first = state;
         for (int i = ghost; i < ghost + cells; ++i)
            for (int c = 0; c < 3; ++c)
               first[i][c] += dt * residual[i][c];

         Residual(first, residual);
         second = state;
         for (int i = ghost; i < ghost + cells; ++i)
            for (int c = 0; c < 3; ++c)
               second[i][c] = 0.75 * state[i][c] +
                              0.25 * (first[i][c] + dt * residual[i][c]);

         Residual(second, residual);
         for (int i = ghost; i < ghost + cells; ++i)
            for (int c = 0; c < 3; ++c)
               state[i][c] = state[i][c] / 3.0 +
                             2.0 / 3.0 *
                                (second[i][c] + dt * residual[i][c]);
         time += dt;
         ++steps;
      }
      return steps;
   }

   double DensityL1Error(double time) const
   {
      if (problem != 4) { return std::numeric_limits<double>::quiet_NaN(); }
      double error = 0.0;
      for (int i = 0; i < cells; ++i)
      {
         const double x = left + (i + 0.5) * dx;
         const double exact = scale *
            (2.0 + 2.0 * std::pow(std::sin(x - time), 2));
         error += std::abs(state[ghost + i][0] - exact) * dx;
      }
      return error;
   }

   void Write(const std::string &path) const
   {
      std::ofstream output(path);
      if (!output) { throw std::runtime_error("cannot open reference output"); }
      output << "x,density,velocity,pressure\n";
      output << std::setprecision(17);
      for (int i = 0; i < cells; ++i)
      {
         const State &value = state[ghost + i];
         output << left + (i + 0.5) * dx << ',' << value[0] << ','
                << value[1] / value[0] << ',' << Pressure(value) << '\n';
      }
   }
};

int IntegerArgument(char **argv, int argc, const std::string &name, int value)
{
   for (int i = 1; i + 1 < argc; ++i)
      if (argv[i] == name) { return std::stoi(argv[i + 1]); }
   return value;
}

double RealArgument(char **argv, int argc, const std::string &name, double value)
{
   for (int i = 1; i + 1 < argc; ++i)
      if (argv[i] == name) { return std::stod(argv[i + 1]); }
   return value;
}

std::string StringArgument(char **argv, int argc, const std::string &name,
                           std::string value)
{
   for (int i = 1; i + 1 < argc; ++i)
      if (argv[i] == name) { return argv[i + 1]; }
   return value;
}

} // namespace

int main(int argc, char **argv)
{
   try
   {
      const int problem = IntegerArgument(argv, argc, "--problem", 7);
      const int cells = IntegerArgument(argv, argc, "--cells", 4096);
      const double scale = RealArgument(argv, argc, "--scale", 1.0);
      const double final_time = RealArgument(
         argv, argc, "--final-time",
         problem == 4 ? 1.1 :
         (problem == 5 ? 0.038 :
          (problem == 6 ? 0.2 : (problem == 7 ? 1.3 : 1.8))));
      const double cfl = RealArgument(argv, argc, "--cfl", 0.4);
      const std::string output = StringArgument(
         argv, argc, "--output", "measurements/reference.csv");
      ReferenceSolver solver(problem, cells, scale);
      const int steps = solver.Run(final_time, cfl);
      solver.Write(output);
      std::cout << "reference_problem=" << problem << " cells=" << cells
                << " steps=" << steps << " time=" << final_time
                << " density_l1_error=" << solver.DensityL1Error(final_time)
                << " output=" << output << '\n';
   }
   catch (const std::exception &error)
   {
      std::cerr << "Reference solver failed: " << error.what() << '\n';
      return 1;
   }
   return 0;
}
