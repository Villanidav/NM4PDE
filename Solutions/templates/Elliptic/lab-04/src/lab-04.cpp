#include <deal.II/base/convergence_table.h>
#include <fstream>
#include <iostream>
#include <vector>

#include "Poisson2D.hpp"

// Main function.
/*
source /u/sw/etc/bash.bashrc
module load gcc-glibc
module load dealii
*/
int
main(int /*argc*/, char * /*argv*/[])
{
  ConvergenceTable table;

  const std::vector<std::string> meshes = {"../mesh/mesh-square-h0.100000.msh"};
  const std::vector<double>      h_vals = {0.1};
  const unsigned int             degree = 2;
  const unsigned int N = 0;

  // Only for Exercise 1:
  std::ofstream convergence_file("convergence.csv");
  convergence_file << "h,eL2,eH1" << std::endl;

  for (unsigned int i = 0; i < meshes.size(); ++i)
    {
      Poisson2D problem(meshes[i], degree, N);

      problem.setup();
      problem.assemble();
      problem.solve();
      problem.output();

      // Only for Exercise 1:
      const double error_L2 = problem.compute_error(VectorTools::L2_norm);
      const double error_H1 = problem.compute_error(VectorTools::H1_norm);

      table.add_value("h", h_vals[i]);
      table.add_value("L2", error_L2);
      table.add_value("H1", error_H1);

      convergence_file << h_vals[i] << "," << error_L2 << "," << error_H1
                       << std::endl;
    }

  // Only for Exercise 1:
  table.evaluate_all_convergence_rates(ConvergenceTable::reduction_rate_log2);
  table.set_scientific("L2", true);
  table.set_scientific("H1", true);
  table.write_text(std::cout);

  return 0;
}
