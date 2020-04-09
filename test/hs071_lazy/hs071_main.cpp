// Testing lazy constraints

#include "restartsqp/LazySqpSolver.hpp"
#include "restartsqp/SqpIpoptTNlp.hpp"
#include "hs071_nlp.hpp"

using namespace Ipopt;
using namespace std;
using namespace RestartSqp;

int main(int argv, char* argc[])
{
  // Create a TNLP -- this is where you create your own TNLP
  SmartPtr<TNLP> ipopt_nlp = new HS071_NLP();

  // Create a wrapper that makes the Ipopt TNLP look like an SqpNlp
  string nlp_name = "ProblemName";
  shared_ptr<SqpIpoptNlp> sqp_nlp =
      make_shared<SqpIpoptNlp>(ipopt_nlp, nlp_name);

  // Create an array with the indices of the constraints that should be considered
  // in the beginning.  Constraints are counted starting from 0.

  // In this case, the hs071 has two constraints.  The first one
  // is an inequality constraints, and the second one is an equality constraint.
  // Here, we pose the problem first without the inequality constraint
  int num_initial_constraints = 1;
  int* constraint_indices = new int[num_initial_constraints]; // allocating size 1 just for demonstration
  constraint_indices[0] = 1;  // include second constraint

  // Now get an SQP solver object that can do the lazy constraint generation
  LazySqpSolver lazy_sqp_solver;

  // Call the optimization method, providing the list of initially considered
  // constraints.
  lazy_sqp_solver.optimize_nlp(sqp_nlp, num_initial_constraints, constraint_indices);

  // At the end, the finalize_solution method for the original Ipopt NLP should be called.
  // If not, we still need to implement this.

  delete [] constraint_indices;

  // Get the exit status of the last SQP solve
  SqpSolverExitStatus exit_flag = lazy_sqp_solver.get_exit_flag();

  return (int) exit_flag;
}
