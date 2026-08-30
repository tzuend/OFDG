#include "mfem.hpp"

#include "../src/kxrcf.hpp"

#include <iostream>

using namespace mfem;

int main(int argc, char *argv[])
{
   int resolution = 96;
   int order = 3;
   int components = 4;
   int repetitions = 100;

   OptionsParser args(argc, argv);
   args.AddOption(&resolution, "-n", "--resolution",
                  "Elements along each coordinate direction.");
   args.AddOption(&order, "-o", "--order", "DG polynomial degree.");
   args.AddOption(&components, "-c", "--components",
                  "Number of pooled state components.");
   args.AddOption(&repetitions, "-r", "--repetitions",
                  "Number of timed indicator evaluations.");
   args.ParseCheck();

   Mesh mesh = Mesh::MakeCartesian2D(resolution, resolution,
                                      Element::QUADRILATERAL,
                                      true, 1.0, 1.0);
   DG_FECollection collection(order, 2);
   FiniteElementSpace space(&mesh, &collection, components,
                            Ordering::byNODES);
   GridFunction state(&space);
   VectorFunctionCoefficient initial(components,
      [components](const Vector &x, Vector &value)
   {
      value.SetSize(components);
      for (int c = 0; c < components; ++c)
      {
         value(c) = 1.0 + 0.05 * c +
                    (x(0) + 0.1 * c * x(1) < 0.5 ? 0.0 : 0.4);
      }
   });
   state.ProjectCoefficient(initial);

   VectorFunctionCoefficient velocity(2, [](const Vector &, Vector &value)
   {
      value.SetSize(2);
      value(0) = 1.0;
      value(1) = 0.25;
   });
   KXRCFIndicator indicator(&space, &velocity);
   Array<bool> active;

   for (int warmup = 0; warmup < 5; ++warmup)
   {
      indicator.Compute(state, active);
   }

   indicator.ResetInternalTimings();
   tic_toc.Clear();
   tic_toc.Start();
   for (int repetition = 0; repetition < repetitions; ++repetition)
   {
      indicator.Compute(state, active);
   }
   tic_toc.Stop();

   int active_count = 0;
   for (int e = 0; e < active.Size(); ++e)
   {
      active_count += active[e] ? 1 : 0;
   }

   std::cout.precision(10);
   std::cout << "elements=" << mesh.GetNE()
             << " order=" << order
             << " components=" << components
             << " repetitions=" << repetitions
             << " active=" << active_count
             << " seconds=" << tic_toc.RealTime()
             << " seconds_per_call=" << tic_toc.RealTime() / repetitions
             << '\n';
   indicator.PrintInternalTimings();
   return 0;
}
