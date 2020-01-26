#include <math.h>
#include <stdio.h>
#include <string.h>

extern "C" {
#include <qpsolver.h>
}

#include "sqphot/MessageHandling.hpp"
#include "sqphot/QoreInterface.hpp"
#include "sqphot/QpOasesInterface.hpp"
#include "sqphot/SparseHbMatrix.hpp"
#include "sqphot/SqpAlgorithm.hpp"
#include "sqphot/Vector.hpp"
#include <memory>

using namespace RestartSqp;
using namespace std;

shared_ptr<SparseHbMatrix>
convert_csr_to_csc(shared_ptr<const SparseHbMatrix> csr_matrix)
{

  double* dense_matrix;
  dense_matrix =
      new double[csr_matrix->get_num_rows() *
                 csr_matrix->get_num_columns()](); // get the dense matrix out

  csr_matrix->get_dense_matrix(dense_matrix);

  bool is_compressed_row = false;
  bool row_oriented = true;
  bool is_symmetric = csr_matrix->is_symmetric();
  shared_ptr<SparseHbMatrix> csc_matrix = make_shared<SparseHbMatrix>(
      csr_matrix->get_num_rows(), csr_matrix->get_num_columns(),
      is_compressed_row, is_symmetric);
  csc_matrix->copy_from_dense_matrix(dense_matrix, csr_matrix->get_num_rows(),
                                     csr_matrix->get_num_columns(),
                                     row_oriented, is_compressed_row);
  delete[] dense_matrix;

  return csc_matrix;
}

shared_ptr<Vector> read_vector(FILE* file, string name)
{
  char string_buf[128];
  int size;
  fscanf(file, "%*s %s %*s %d %*s", string_buf, &size);
  assert(strncmp(string_buf, name.c_str(), 128) == 0);

  shared_ptr<Vector> retval = make_shared<Vector>(size);
  double* vec_vals = retval->get_non_const_values();
  for (int i = 0; i < size; ++i) {
    fscanf(file, "%le", &vec_vals[i]);
  }

  return retval;
}

shared_ptr<SparseTripletMatrix> read_triplet_matrix(FILE* file, int nRows,
                                                    int nCols,
                                                    bool is_symmetric,
                                                    string name)
{
  char string_buf[128];
  int nnz;
  fscanf(file, "%*s %s %*s %d %*s %*s", string_buf, &nnz);
  assert(strncmp(string_buf, name.c_str(), 128) == 0);

  shared_ptr<SparseTripletMatrix> retval =
      make_shared<SparseTripletMatrix>(nnz, nRows, nCols, is_symmetric);
  int* row_indices = retval->get_nonconst_row_indices();
  int* columns_indices = retval->get_nonconst_column_indices();
  double* values = retval->get_nonconst_values();

  for (int i = 0; i < nnz; ++i) {
    fscanf(file, "%d %d %le", &row_indices[i], &columns_indices[i], &values[i]);
  }

  return retval;
}

int main(int argc, char* argv[])
{

  FILE* file = fopen(argv[1], "r");
  char buffer[100];
  qp_int tmp_int;
  double tmp_double;
  int nVar, nCon, Hnnz,
      Annz; // number of variables, number of constraints, number of nonzeros
  string exitflag_qore, exitflag_qpOASES;
  int i;
  // for H and A;
  if (argc < 2) {
    printf("Missing Filename\n");
    return 1;
  } else if (file == NULL)
    perror("Error opening file");

  shared_ptr<Vector> g;

  shared_ptr<SparseTripletMatrix> A_triplet;
  shared_ptr<SparseTripletMatrix> H_triplet;

  shared_ptr<Vector> lb;
  shared_ptr<Vector> ub;
  shared_ptr<Vector> lbA;
  shared_ptr<Vector> ubA;

  // Find out whether this is a file in the "new" format where the matrix is in
  // triplet format
  bool read_triplet_format = false;
  if (fgets(buffer, 100, file) != NULL) {
    if (strncmp(buffer, "Vector", 6) == 0) {
      read_triplet_format = true;
    }
  }
  // Rewind the file pointer
  fseek(file, 0L, SEEK_SET);

  if (read_triplet_format) {
    g = read_vector(file, "g");
    nVar = g->get_dim();
    lb = read_vector(file, "lb");
    assert(lb->get_dim() == nVar);
    ub = read_vector(file, "ub");
    assert(ub->get_dim() == nVar);

    lbA = read_vector(file, "lbA");
    nCon = lbA->get_dim();
    ubA = read_vector(file, "ubA");
    assert(ubA->get_dim() == nCon);

    H_triplet = read_triplet_matrix(file, nVar, nVar, true, "H");
    A_triplet = read_triplet_matrix(file, nCon, nVar, false, "A");
  } else {
    // begin Reading Data
    if (fgets(buffer, 100, file) != NULL) {
      sscanf(buffer, "%d", &nVar);
    }
    if (fgets(buffer, 100, file) != NULL) {
      sscanf(buffer, "%d", &nCon);
    }
    if (fgets(buffer, 100, file) != NULL) {
      sscanf(buffer, "%d", &Annz);
    }
    if (fgets(buffer, 100, file) != NULL) {
      sscanf(buffer, "%d", &Hnnz);
    }

    g = make_shared<Vector>(nVar);
    shared_ptr<Vector> all_lb = make_shared<Vector>(nCon + nVar);
    shared_ptr<Vector> all_ub = make_shared<Vector>(nCon + nVar);
    shared_ptr<SparseHbMatrix> A_csr =
        make_shared<SparseHbMatrix>(Annz, nCon, nVar, true, false);
    shared_ptr<SparseHbMatrix> H_csr =
        make_shared<SparseHbMatrix>(Hnnz, nVar, nVar, true, true);

    for (i = 0; i < nCon + nVar; i++) {
      if (fgets(buffer, 100, file) != NULL) {
        sscanf(buffer, "%lf", &tmp_double);
        all_lb->set_value(i, tmp_double);
      }
    }

    for (i = 0; i < nCon + nVar; i++) {
      if (fgets(buffer, 100, file) != NULL) {
        sscanf(buffer, "%lf", &tmp_double);
        all_ub->set_value(i, tmp_double);
      }
    }

    for (i = 0; i < nVar; i++) {
      if (fgets(buffer, 100, file) != NULL) {
        sscanf(buffer, "%lf", &tmp_double);
        g->set_value(i, tmp_double);
      }
    }

    for (i = 0; i < nCon + 1; i++) {
      if (fgets(buffer, 100, file) != NULL) {
        sscanf(buffer, "%d", &tmp_int);
        A_csr->set_row_index_at_entry(i, tmp_int);
      }
    }

    for (i = 0; i < Annz; i++) {
      if (fgets(buffer, 100, file) != NULL) {
        sscanf(buffer, "%d", &tmp_int);
        A_csr->set_column_index_at_entry(i, tmp_int);
      }
    }

    for (i = 0; i < Annz; i++) {
      if (fgets(buffer, 100, file) != NULL) {
        sscanf(buffer, "%lf", &tmp_double);
        A_csr->set_value_at_entry(i, tmp_double);
      }
    }

    for (i = 0; i < nVar + 1; i++) {
      if (fgets(buffer, 100, file) != NULL) {
        sscanf(buffer, "%d", &tmp_int);
        H_csr->set_row_index_at_entry(i, tmp_int);
      }
    }

    for (i = 0; i < Hnnz; i++) {
      if (fgets(buffer, 100, file) != NULL) {
        sscanf(buffer, "%d", &tmp_int);
        H_csr->set_column_index_at_entry(i, tmp_int);
      }
    }

    for (i = 0; i < Hnnz; i++) {
      if (fgets(buffer, 100, file) != NULL) {
        sscanf(buffer, "%lf", &tmp_double);
        H_csr->set_value_at_entry(i, tmp_double);
      }
    }

    A_csr->print("A_csr");
    H_csr->print("H_csr");

    A_triplet = make_shared<SparseTripletMatrix>(A_csr);
    H_triplet = make_shared<SparseTripletMatrix>(H_csr);

    lb = make_shared<Vector>(nVar, all_lb->get_values());
    ub = make_shared<Vector>(nVar, all_ub->get_values());
    lbA = make_shared<Vector>(nCon, all_lb->get_values() + nVar);
    ubA = make_shared<Vector>(nCon, all_ub->get_values() + nVar);
  }

  lb->print("lb");
  ub->print("ub");
  lbA->print("lbA");
  ubA->print("ubA");
  g->print("g");
  A_triplet->print("A_triplet");
  H_triplet->print("H_triplet");

  IdentityMatrixPositions identity_positions;

  /////////////////////////////////////////////////////
  //  Create a dummy Algorithm object to get options //
  /////////////////////////////////////////////////////

  shared_ptr<SqpAlgorithm> dummy_algorithm = make_shared<SqpAlgorithm>();
  Ipopt::SmartPtr<Ipopt::Journalist> jnlst = dummy_algorithm->get_jnlst();
  Ipopt::SmartPtr<Ipopt::OptionsList> options =
      dummy_algorithm->get_options_list();

  // read matrix data
  ///////////////////////////////////////////////////////////
  //                     QORE                              //
  ///////////////////////////////////////////////////////////

  //    NLPInfo nlp_info;
  //    nlp_info.nCon = nCon;
  //    nlp_info.nVar = nVar-2*nCon;
  //    nlp_info.nnz_jac_g = Annz;
  //    nlp_info.nnz_h_lag = Hnnz;

  shared_ptr<Statistics> stats_qore = make_shared<Statistics>();

  shared_ptr<QoreInterface> qore_interface =
      make_shared<QoreInterface>(nVar, nCon, QP_TYPE_QP, options, jnlst);
  qore_interface->get_lower_variable_bounds_nonconst()->copy_vector(lb);
  qore_interface->get_upper_variable_bounds_nonconst()->copy_vector(ub);
  qore_interface->get_lower_constraint_bounds_nonconst()->copy_vector(lbA);
  qore_interface->get_upper_constraint_bounds_nonconst()->copy_vector(ubA);
  qore_interface->get_linear_objective_coefficients_nonconst()->copy_vector(g);
  qore_interface->set_constraint_jacobian(A_triplet, identity_positions);
  qore_interface->set_objective_hessian(H_triplet);

  try {
    qore_interface->optimize(stats_qore);
  } catch (...) {
  }

  shared_ptr<Vector> x_qore;
  shared_ptr<Vector> y_qore_constr;
  shared_ptr<Vector> y_qore_bounds;
  double obj_qore;
  KktError kkt_error_qore;

  QpSolverExitStatus qore_exit_status = qore_interface->get_solver_status();
  if (qore_exit_status == QPEXIT_OPTIMAL) {
    x_qore = make_shared<Vector>(
        nVar, qore_interface->get_primal_solution()->get_values());
    y_qore_constr = make_shared<Vector>(
        nCon, qore_interface->get_constraint_multipliers()->get_values());
    y_qore_bounds = make_shared<Vector>(
        nVar, qore_interface->get_bounds_multipliers()->get_values());
    obj_qore = qore_interface->get_optimal_objective_value();
    kkt_error_qore = qore_interface->calc_kkt_error();
  } else {
    x_qore = make_shared<Vector>(nVar);
    y_qore_constr = make_shared<Vector>(nCon);
    y_qore_bounds = make_shared<Vector>(nVar);
    x_qore->set_to_zero();
    y_qore_constr->set_to_zero();
    y_qore_bounds->set_to_zero();
    obj_qore = 0.;
  }

  switch (qore_exit_status) {
    case QPEXIT_OPTIMAL:
      exitflag_qore = "OPTIMAL";
      break;
    case QPEXIT_INFEASIBLE:
      exitflag_qore = "INFEASIBLE";
      break;
    case QPEXIT_UNBOUNDED:
      exitflag_qore = "UNBOUNDED";
      break;
    case QPEXIT_INTERNAL_ERROR:
      exitflag_qore = "INTERNAL_ERROR";
      break;
    default:
      exitflag_qore = "UNKNOWN";
  }

  ///////////////////////////////////////////////////////////
  //                     QPOASES                           //
  ///////////////////////////////////////////////////////////
  shared_ptr<Statistics> stats_qpOASES = make_shared<Statistics>();

  shared_ptr<QpOasesInterface> qpoases_interface =
      make_shared<QpOasesInterface>(nVar, nCon, QP_TYPE_QP, options, jnlst);
  qpoases_interface->get_lower_variable_bounds_nonconst()->copy_vector(lb);
  qpoases_interface->get_upper_variable_bounds_nonconst()->copy_vector(ub);
  qpoases_interface->get_lower_constraint_bounds_nonconst()->copy_vector(lbA);
  qpoases_interface->get_upper_constraint_bounds_nonconst()->copy_vector(ubA);
  qpoases_interface->get_linear_objective_coefficients_nonconst()->copy_vector(
      g);
  qpoases_interface->set_constraint_jacobian(A_triplet, identity_positions);
  qpoases_interface->set_objective_hessian(H_triplet);

  try {
    qpoases_interface->optimize(stats_qpOASES);
  } catch (...) {
  }

#if 0
  shared_ptr<Vector> x_qpOASES = make_shared<Vector>(
      nVar, qpoases_interface->get_primal_solution()->get_values());
  shared_ptr<Vector> y_qpOASES_constr = make_shared<Vector>(
      nCon, qpoases_interface->get_constraint_multipliers()->get_values());
  shared_ptr<Vector> y_qpOASES_bounds = make_shared<Vector>(
      nVar, qpoases_interface->get_bounds_multipliers()->get_values());
#endif
  shared_ptr<Vector> x_qpOASES;
  shared_ptr<Vector> y_qpOASES_constr;
  shared_ptr<Vector> y_qpOASES_bounds;
  double obj_qpOASES;
  KktError kkt_error_qpoases;

  QpSolverExitStatus qpOASES_exit_status =
      qpoases_interface->get_solver_status();
  if (qpOASES_exit_status == QPEXIT_OPTIMAL) {
    x_qpOASES = make_shared<Vector>(
        nVar, qpoases_interface->get_primal_solution()->get_values());
    y_qpOASES_constr = make_shared<Vector>(
        nCon, qpoases_interface->get_constraint_multipliers()->get_values());
    y_qpOASES_bounds = make_shared<Vector>(
        nVar, qpoases_interface->get_bounds_multipliers()->get_values());
    obj_qpOASES = qpoases_interface->get_optimal_objective_value();
    kkt_error_qpoases = qpoases_interface->calc_kkt_error();
  } else {
    x_qpOASES = make_shared<Vector>(nVar);
    y_qpOASES_constr = make_shared<Vector>(nCon);
    y_qpOASES_bounds = make_shared<Vector>(nVar);
    x_qpOASES->set_to_zero();
    y_qpOASES_constr->set_to_zero();
    y_qpOASES_bounds->set_to_zero();
    obj_qpOASES = 0.;
  }

  //  x_qpOASES->copy_vector(qpoases_interface->get_primal_solution());
  //  y_qpOASES_constr->copy_vector(
  //      qpoases_interface->get_constraint_multipliers());
  //  y_qpOASES_bounds->copy_vector(qpoases_interface->get_bounds_multipliers());

  switch (qpoases_interface->get_solver_status()) {
    case QPEXIT_OPTIMAL:
      exitflag_qpOASES = "OPTIMAL";
      break;
    case QPEXIT_INFEASIBLE:
      exitflag_qpOASES = "INFEASIBLE";
      break;
    case QPEXIT_UNBOUNDED:
      exitflag_qpOASES = "UNBOUNDED";
      break;
    case QPEXIT_INTERNAL_ERROR:
      exitflag_qpOASES = "INTERNAL_ERROR";
      break;
#if 0
    case QPERROR_EXCEED_MAX_ITER:
      exitflag_qpOASES = "EXCEED_ITER_LIMIT";
      break;
    case QPERROR_PREPARINGAUXILIARYQP:
      exitflag_qpOASES = "PREPARINGAUXILIARYQP";
      break;
    case QPERROR_AUXILIARYQPSOLVED:
      exitflag_qpOASES = "AUXILIARYQP_SOLVED";
      break;
    case QPERROR_PERFORMINGHOMOTOPY:
      exitflag_qpOASES = "PERFORMING_HOMOTOPY";
      break;
    case QPERROR_HOMOTOPYQPSOLVED:
      exitflag_qpOASES = "HOMOTOPYQPSOLVED";
      break;
#endif
    default:
      exitflag_qpOASES = "UNKNOWN";
  }

  ///////////////////////////////////////////////////////////
  //                     PRINT OUT RESULTS                 //
  ///////////////////////////////////////////////////////////
  string* pname = new string(argv[1]);
  std::size_t found = pname->find_last_of("/\\");
  printf(DOUBLE_LONG_DIVIDER);
  printf("                         Solving QP Problem stored in %10s\n",
         pname->substr(found + 1).c_str());
  printf(DOUBLE_LONG_DIVIDER);
  printf("%30s    %23s    %23s\n", "", "QORE", "QPOASES");
  printf("%30s    %23s    %23s\n", "Exitflag", exitflag_qore.c_str(),
         exitflag_qpOASES.c_str());
  printf("%30s    %23d    %23d\n", "Iteration", stats_qore->num_qp_iterations_,
         stats_qpOASES->num_qp_iterations_);
  printf("%30s    %23.16e    %23.16e\n", "Objective", obj_qore, obj_qpOASES);
  printf("%30s    %23.16e    %23.16e\n", "||x||", x_qore->calc_one_norm(),
         x_qpOASES->calc_one_norm());
  printf("%30s    %23.16e    %23.16e\n", "||y_b||",
         y_qore_bounds->calc_inf_norm(), y_qpOASES_bounds->calc_inf_norm());
  printf("%30s    %23.16e    %23.16e\n", "||y_c||",
         y_qore_constr->calc_inf_norm(), y_qpOASES_constr->calc_inf_norm());
  printf("%30s    %23.16e    %23.16e\n", "KKT Error",
         kkt_error_qore.worst_violation, kkt_error_qpoases.worst_violation);
  printf("%30s    %23.16e    %23.16e\n", "Primal Feasibility Violation",
         kkt_error_qore.primal_infeasibility,
         kkt_error_qpoases.primal_infeasibility);
  printf("%30s    %23.16e    %23.16e\n", "Dual Feasibility  Violation",
         kkt_error_qore.dual_infeasibility,
         kkt_error_qpoases.dual_infeasibility);
  printf("%30s    %23.16e    %23.16e\n", "Complemtarity Violation",
         kkt_error_qore.complementarity_violation,
         kkt_error_qpoases.complementarity_violation);
  printf("%30s    %23.16e    %23.16e\n", "Working Set Violation",
         kkt_error_qore.working_set_error, kkt_error_qpoases.working_set_error);
  printf(DOUBLE_LONG_DIVIDER);

  delete pname;
  fclose(file);
  return 0;
}
