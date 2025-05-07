#include <deal.II/base/convergence_table.h>

#include <fstream>
#include <iostream>
#include <vector>
#include "Heat.hpp"

// Main function.
int
main(int argc, char *argv[])
{
  ConvergenceTable table;
  Utilities::MPI::MPI_InitFinalize mpi_init(argc, argv);
  const unsigned int               mpi_rank =
  Utilities::MPI::this_mpi_process(MPI_COMM_WORLD);

  const unsigned int degree = 1;
  const double T     = 1.0;
  const double theta = 1;
  const std::vector<unsigned int> N_values = {19};

  const std::vector<double> deltat_vector = {
    0.05};
  std::vector<double> errors_L2;
  std::vector<double> errors_H1;
  
  std::ofstream convergence_file("convergence.csv");
  convergence_file << "h,eL2,eH1" << std::endl;

  for( const double &deltat : deltat_vector){
  for (const unsigned int &N : N_values)
    {
      Heat problem(N, degree, T, deltat, theta);

      problem.setup();
      problem.solve();

      const double h        = 1.0 / (N + 1.0);
      const double error_L2 = problem.compute_error(VectorTools::L2_norm);
      const double error_H1 = problem.compute_error(VectorTools::H1_norm);

      table.add_value("h", h);
      table.add_value("dt", deltat);
      table.add_value("L2", error_L2);
      table.add_value("H1", error_H1);

      convergence_file << h << "," << error_L2 << "," << error_H1 << std::endl;
    }

  table.evaluate_all_convergence_rates(ConvergenceTable::reduction_rate_log2);

  table.set_scientific("L2", true);
  table.set_scientific("H1", true);

  table.write_text(std::cout);
  }
  return 0;
}