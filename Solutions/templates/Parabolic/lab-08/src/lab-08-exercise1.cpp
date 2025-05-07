#include <deal.II/base/convergence_table.h>

#include <fstream>
#include <iostream>
#include <vector>
#include "Heat.hpp"

// Main function.
/*
source /u/sw/etc/bash.bashrc
module load gcc-glibc
module load dealii

  -nel main inserisci i parametri corretti, degree, theta ...
  -cambia funzioni in .hpp e dimensione/mesh
  -cambia dirichlet neumann in .cpp
  -controlla che la weak formulation è coerente con quella della toeria
  (in particolare cambiare se ho -bugrad(v)) per le neumann

  */
int
main(int argc, char *argv[])
{
  // MPI init
  Utilities::MPI::MPI_InitFinalize mpi_init(argc, argv);
  const unsigned int mpi_rank = Utilities::MPI::this_mpi_process(MPI_COMM_WORLD);

  dealii::ConvergenceTable table;

  // Problem parameters
  const unsigned int degree = 1;
  const double       T      = 1.0;
  const double       theta  = 0.5;

  // Vector of time steps to test
  const std::vector<double> deltat_vector = {0.05, 0.025, 0.0125, 0.00625, 0.003125};
  //const std::vector<double> deltat_vector = {0.05};

  // For each \Delta t, solve PDE & store errors
  for (const double deltat : deltat_vector)
    {
      Heat problem("../mesh/mesh-square-h0.100000.msh", degree, T, deltat, theta);

      problem.setup();
      problem.solve();

      const double eL2 = problem.compute_error(dealii::VectorTools::L2_norm);
      const double eH1 = problem.compute_error(dealii::VectorTools::H1_norm);

      table.add_value("dt",   deltat);
      table.add_value("e_L2", eL2);
      table.add_value("e_H1", eH1);
    }

  // Evaluate slopes in log base 2:
  table.evaluate_convergence_rates("e_L2", "dt",
                                   dealii::ConvergenceTable::reduction_rate_log2,
                                   dealii::ConvergenceTable::reduction_rate_log2);

  table.evaluate_convergence_rates("e_H1", "dt",
                                   dealii::ConvergenceTable::reduction_rate_log2,
                                   dealii::ConvergenceTable::reduction_rate_log2);

  // Set formatting (precision, scientific notation, etc.)
  table.set_precision("dt",   4);
  table.set_scientific("dt",  true);

  table.set_precision("e_L2", 4);
  table.set_scientific("e_L2", true);

  table.set_precision("e_H1", 4);
  table.set_scientific("e_H1", true);

  if (mpi_rank == 0)
    {
      std::cout << "\nConvergence Table:\n";
      // Just call write_text(...) with one or two arguments:
      table.write_text(std::cout); // or use: (std::cout, TableHandler::org_mode_table);

      std::ofstream out("convergence_table.txt");
      table.write_text(out);       // or use: (out, TableHandler::org_mode_table);
    }

  return 0;
}