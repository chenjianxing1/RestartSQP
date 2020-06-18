/* Copyright (C) 2019
 * All Rights Reserved.
 *
 * Authors: Xinyi Luo
 * Date:2019-06
 */

#include "restartsqp/LazySqpSolver.hpp"
#include "restartsqp/MessageHandling.hpp"

#include <algorithm>
#include <fstream>

using namespace std;
using namespace Ipopt;

namespace RestartSqp {

LazySqpSolver::LazySqpSolver()
 : considered_constraints_indices_(nullptr)
 , orig_constraint_lower_bounds_(nullptr)
 , orig_constraint_upper_bounds_(nullptr)
{
  //  register_options_(reg_options_, skip_ipopt_options);

  // Create the cross-over NLP solver
  crossover_sqp_solver_ = make_shared<CrossoverSqpSolver>();

  // Get the journalist
  jnlst_ = crossover_sqp_solver_->get_jnlst();

  // Get the options TOOO
  num_max_lazy_nlp_solves_ = 20;
  max_add_lazy_constraints_ = 50;
}

/**
 * Destructor
 */
LazySqpSolver::~LazySqpSolver()
{
  delete[] considered_constraints_indices_;
  delete[] orig_constraint_lower_bounds_;
  delete[] orig_constraint_upper_bounds_;
}

void LazySqpSolver::optimize_nlp(shared_ptr<SqpTNlp> orig_sqp_tnlp,
                                 int num_initial_constraints,
                                 int* initial_constraint_indices,
                                 const string& options_file_name)
{
  // Initailize the exit_flag
  exit_flag_ = UNKNOWN_EXIT_STATUS;

  // Set the handle for the original NLP
  orig_sqp_tnlp_ = orig_sqp_tnlp;

  // Create the wrapper SqpTNlp object that can switch off constraints
  lazy_sqp_tnlp_ = make_shared<LazySqpTNlp>(orig_sqp_tnlp_);

  // Get information about the size of the original problem
  int num_nonzeros_jacobian;
  int num_nonzeros_hessian;
  string nlp_name;
  bool retval;
  retval = orig_sqp_tnlp_->get_nlp_info(
      num_orig_variables_, num_orig_constraints_, num_nonzeros_jacobian,
      num_nonzeros_hessian, nlp_name);
  if (!retval) {
    jnlst_->Printf(J_ERROR, J_MAIN, "Error calling get_nlp_info.\n");
    exit_flag_ = INVALID_NLP;
    return;
  }

  // Get the constraint bounds (so that later we can check which constraints are
  // violated. */
  double* variable_lower_bounds = new double[num_orig_variables_];
  double* variable_upper_bounds = new double[num_orig_variables_];
  orig_constraint_lower_bounds_ = new double[num_orig_constraints_];
  orig_constraint_upper_bounds_ = new double[num_orig_constraints_];

  retval = orig_sqp_tnlp_->get_bounds_info(
      num_orig_variables_, variable_lower_bounds, variable_upper_bounds,
      num_orig_constraints_, orig_constraint_lower_bounds_,
      orig_constraint_upper_bounds_);
  delete[] variable_lower_bounds;
  delete[] variable_upper_bounds;
  if (!retval) {
    jnlst_->Printf(J_ERROR, J_MAIN, "Error calling get_bounds_info.\n");
    exit_flag_ = INVALID_NLP;
    return;
  }

  // Initialize the set of constraints included in the lazy NLP with those
  // given
  // by the caller
  delete[] considered_constraints_indices_;
  considered_constraints_indices_ = nullptr;
  considered_constraints_indices_ = new int[num_orig_constraints_];
  num_considered_constraints_ = num_initial_constraints;
  for (int i = 0; i < num_considered_constraints_; ++i) {
    considered_constraints_indices_[i] = initial_constraint_indices[i];
  }

  // Give the set of constraints to the NLP object
  lazy_sqp_tnlp_->set_considered_constraints(num_considered_constraints_,
                                             considered_constraints_indices_);

  // Solve the first NLP with the initial set of constraints
  jnlst_->Printf(J_SUMMARY, J_MAIN,
                 "\n===== LazyNlpSolver: Solve the initial NLP\n\n");
  crossover_sqp_solver_->crossover_solve(lazy_sqp_tnlp_, options_file_name);

  // Check the exit flag
  exit_flag_ = crossover_sqp_solver_->get_exit_flag();
  if (exit_flag_ != OPTIMAL) {
    jnlst_->Printf(J_ERROR, J_MAIN, "Error after initial NLP solve: %d\n",
                   exit_flag_);
    return;
  }

  // Start the loop that solves the sequence of NLPs for which constraints are
  // added
  int iter_counter = 0;

  while (true) {
    // Increase iteration counter and check if the maximum number of iterations
    // is exceeded
    iter_counter++;
    if (iter_counter > num_max_lazy_nlp_solves_) {
      jnlst_->Printf(
          J_ERROR, J_MAIN,
          "Number of maximum number of lazy NLP solves (%d) exceeded.\n",
          num_max_lazy_nlp_solves_);
      exit_flag_ = EXCEED_MAX_LAZY_NLP_SOLVES;
      break;
    }

    jnlst_->Printf(J_SUMMARY, J_MAIN,
                   "\n===== Starting iteration %d in the lazy NLP solver.\n\n",
                   iter_counter);

    // Call method to update the considered constraints
    int num_added_constraints = update_considered_constraints_();

    // Check if new constraints need to be added
    if (num_added_constraints < 0) {
      jnlst_->Printf(J_ERROR, J_MAIN, "Adding constraints failed.\n");
      exit_flag_ = ERROR_IN_LAZY_NLP_UPDATE;
      return;
    } else if (num_added_constraints == 0) {
      jnlst_->Printf(
          J_SUMMARY, J_MAIN,
          "No lazy constraints added, terminating lazy NLP solver loop.\n");
      break;
    }

    // Solve the augmented NLP
    crossover_sqp_solver_->solve(lazy_sqp_tnlp_);
  }
}

typedef pair<int, double> const_viol;

static bool compare_const_viol(const const_viol& rhs, const const_viol& lhs)
{
  return (rhs.second < lhs.second);
}

int LazySqpSolver::update_considered_constraints_()
{

  // Get the solution of the NLP
  shared_ptr<const Vector> primal_variables =
      crossover_sqp_solver_->get_current_primal_variables();
  const double* primal_variables_values = primal_variables->get_values();

  // Compute the values of the full constraint set
  double* orig_constraint_values = new double[num_orig_constraints_];
  bool new_primal_variables = true;
  bool retval = orig_sqp_tnlp_->eval_constraint_values(
      num_orig_variables_, primal_variables_values, new_primal_variables,
      num_orig_constraints_, orig_constraint_values);
  assert(retval && "original constraints could not be evaluated");

  // Compute the constraint violations and store them in a list
  double constraint_violation_tol = 1e-6; // TODO make option
  list<const_viol> constraint_violations;
  for (int i = 0; i < num_orig_constraints_; ++i) {
    // Compute constraint violation
    double violation =
        max(0., orig_constraint_lower_bounds_[i] - orig_constraint_values[i]);
    violation = max(violation, orig_constraint_values[i] -
                                   orig_constraint_upper_bounds_[i]);
    // printf("i=%d lower=%e body=%e upper=%e\n", i,
    // orig_constraint_lower_bounds_[i],
    //       orig_constraint_values[i], orig_constraint_upper_bounds_[i]);

    // If there is a constraint violation, we add the new violation to the list.
    // We need to correct the counting so that it starts from 1.
    if (violation > constraint_violation_tol) {
      const_viol this_entry(i, violation);
      constraint_violations.push_back(this_entry);
    }
  }

  delete[] orig_constraint_values;

  // Determine the number of violated constraints
  int num_violated_constraints = constraint_violations.size();

  jnlst_->Printf(J_SUMMARY, J_MAIN, "\nNumber of violated constraints: %d\n",
                 num_violated_constraints);

  int num_constraints_added;

  // Select the top max_add_lazy_constraints constraints to add
  if (num_violated_constraints > 0) {
    // Sort the violated constraints by violation
    constraint_violations.sort(compare_const_viol);

    // Reserve memory for the top constraints
    num_constraints_added =
        min(num_violated_constraints, max_add_lazy_constraints_);
    int* constraints_to_add = new int[num_constraints_added];

    int counter = 0;
    for (auto it = constraint_violations.begin();
         it != constraint_violations.end(); ++it) {

      // Add the constraint index
      constraints_to_add[counter] = (*it).first;

      // Stop the loop if enough constraints have been added
      counter++;
      if (counter == num_constraints_added) {
        break;
      }
    }

    // Add the constraints to the NLP object
    retval = lazy_sqp_tnlp_->add_new_constraints(num_constraints_added,
                                                 constraints_to_add);

    jnlst_->Printf(J_SUMMARY, J_MAIN,
                   "\n%d constraints added:", num_constraints_added);
    for (int i = 0; i < num_constraints_added; ++i) {
      jnlst_->Printf(J_SUMMARY, J_MAIN, " %d", constraints_to_add[i]);
    }
    jnlst_->Printf(J_SUMMARY, J_MAIN, "\n\n");

    delete[] constraints_to_add;

    if (!retval) {
      num_constraints_added = -1;
    }
  } else {
    jnlst_->Printf(J_SUMMARY, J_MAIN, "\nAll constraints satisfied.\n\n");
    num_constraints_added = 0;
  }

  return num_constraints_added;
}

} // namespace RestartSqp
