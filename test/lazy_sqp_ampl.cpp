#include "AmplTNLP.hpp"
#include "IpIpoptApplication.hpp"
#include "IpTNLP.hpp"
#include "restartsqp/LazySqpSolver.hpp"
#include "restartsqp/SqpIpoptTNlp.hpp"
#include "restartsqp/Utils.hpp"
#include <cstddef>
#include <iostream>
#include <memory>
#include <stdio.h>
#include <string.h>
#include <string>
#include <sys/stat.h>

using namespace Ipopt;
using namespace RestartSqp;
using namespace std;

inline bool exist(const std::string& name)
{
  struct stat buffer;
  return (stat(name.c_str(), &buffer) == 0);
}

class Table_Writer
{

public:
  Table_Writer(const std::string& filename)
  {

    if (!exist(filename)) {
      file = fopen(filename.c_str(), "a"); // append at the end
      fprintf(file, "%10s   %10s    %10s    %10s    %10s    %10s    ", "name",
              "nVar", "nConstr", "iter", "QP_iter", "exitflag");

      fprintf(file, "%23s    %23s    %23s    %23s    %23s\n", "objective",
              "primal_violation", "dual_violation", "stationarity_violation",
              "compl_violation");

    } else {
      file = fopen(filename.c_str(), "a"); // append at the end
    }
  }

  ~Table_Writer()
  {
    fclose(file);
  }

  /**
   * @brief Write a brief summary for each problem being solved,
   */
  void write_in_brief(const std::string& pname, /**<name of the problem*/
                      CrossoverSqpSolver& alg)
  {

    std::size_t found = pname.find_last_of("/\\");
    KktError kkt_error = alg.get_kkt_error();
    shared_ptr<Statistics> stats;
    stats = alg.get_stats();

    fprintf(file, "%12s  %6d  %6d   %6d   %6d   %6d   ",
            pname.substr(found + 1).c_str(), alg.get_num_var(),
            alg.get_num_constr(), stats->num_sqp_iterations_,
            stats->num_qp_iterations_, alg.get_exit_flag());
#if 1
    fprintf(file, "%23.16e  %9.2e  %9.2e  %9.2e  %9.2e\n",
            alg.get_final_objective(), kkt_error.primal_infeasibility,
            kkt_error.dual_infeasibility, kkt_error.complementarity_violation,
            kkt_error.working_set_error);
#endif
  }

private:
  FILE* file;
};

int main(int argc, char** args)
{
  // Create the SQP algorithm object
  LazySqpSolver alg;

  // Create an Ipopt::AmplTNLP object to handle the AMPL model
  SmartPtr<OptionsList> dummy_options = new OptionsList();
  SmartPtr<TNLP> ampl_tnlp =
      new AmplTNLP(ConstPtr(alg.get_jnlst()), dummy_options, args);
  string nlp_name(args[1]);
  shared_ptr<SqpIpoptTNlp> sqp_nlp =
      make_shared<SqpIpoptTNlp>(ampl_tnlp, nlp_name);

  // Get the number of constraints
  int num_variables;
  int num_constraints;
  int num_nonzeros_jacobian;
  int num_nonzeros_hessian;
  sqp_nlp->get_nlp_info(num_variables, num_constraints, num_nonzeros_jacobian,
                        num_nonzeros_hessian, nlp_name);

  assert(num_constraints > 0);

  // Memory for the initially chosen constraints
  int num_initial_constraints; // num_constraints;
  int* constraint_indices = new int[num_constraints];

#if 1
  // First solve the problem correctly
  CrossoverSqpSolver crossover;
  crossover.crossover_solve(sqp_nlp);

  // Get the activity status for the constraints and
  const ActivityStatus* constraint_activity_status =
      crossover.get_constraints_working_set();
  num_initial_constraints = 0;
  for (int i = 0; i < num_constraints; ++i) {
    // printf("activity %d = %d\n", i, constraint_activity_status[i]);
    if (constraint_activity_status[i] != INACTIVE) {
      constraint_indices[num_initial_constraints] = i + 1;
      // printf("num_initial_constraints = %d\n", num_initial_constraints);
      num_initial_constraints++;
    }
  }

  // Take the last X constraints off
  num_initial_constraints -= 2;
  assert(num_initial_constraints >= 0);

#else
  // Set the initial number of constraints
  num_initial_constraints = 1; // num_constraints;
  for (int i = 0; i < num_initial_constraints; ++i) {
    constraint_indices[i] = i + 1;
  }
#endif

  // Solve the AMPL model
  // alg.initialize(sqp_nlp, args[1]);
  alg.optimize_nlp(sqp_nlp, num_initial_constraints, constraint_indices);

#if 0
  shared_ptr<Table_Writer> writer =
      make_shared<Table_Writer>("result_table.txt");
  writer->write_in_brief(args[1], alg);
#endif

  delete[] constraint_indices;

  return 0;
}
