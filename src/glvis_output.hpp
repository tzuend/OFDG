#ifndef DOFDG_GLVIS_OUTPUT_HPP
#define DOFDG_GLVIS_OUTPUT_HPP

#include "mfem.hpp"

#include <iostream>
#include <string>

// Shared wrapper for optional live GLVis output. Keeping the socket protocol
// here avoids duplicating MPI-specific visualization code in every example.
class GLVisOutput
{
private:
   mfem::socketstream stream_;
   int process_count_ = 1;
   int rank_ = 0;
   bool enabled_ = false;

public:
   GLVisOutput(const mfem::ParMesh &mesh, bool requested, int precision = 8,
               const char *host = "localhost", int port = 19916)
   {
      MPI_Comm_size(mesh.GetComm(), &process_count_);
      MPI_Comm_rank(mesh.GetComm(), &rank_);
      if (!requested) { return; }

      stream_.open(host, port);
      const int locally_connected = stream_ ? 1 : 0;
      int globally_connected = 0;
      MPI_Allreduce(&locally_connected, &globally_connected, 1, MPI_INT,
                    MPI_MIN, mesh.GetComm());
      enabled_ = globally_connected != 0;
      if (!enabled_)
      {
         if (rank_ == 0)
         {
            std::cout << "Unable to connect every MPI rank to the GLVis server at "
                      << host << ':' << port << ".\n"
                      << "GLVis visualization disabled.\n";
         }
         return;
      }
      stream_.precision(precision);
   }

   bool Enabled() const { return enabled_; }

   void Send(const mfem::ParMesh &mesh,
             const mfem::ParGridFunction &solution,
             const std::string &window_title,
             const std::string &initial_commands = std::string())
   {
      if (!enabled_) { return; }
      stream_ << "parallel " << process_count_ << ' ' << rank_ << '\n';
      stream_ << "solution\n" << mesh << solution;
      if (!window_title.empty())
      {
         stream_ << "window_title '" << window_title << "'\n";
      }
      if (!initial_commands.empty()) { stream_ << initial_commands; }
      stream_ << std::flush;
   }
};

#endif // DOFDG_GLVIS_OUTPUT_HPP
