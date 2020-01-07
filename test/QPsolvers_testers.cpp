#include <math.h>
#include <stdio.h>
#include <string.h>

extern "C" {
#include <qpsolver.h>
}

#include "sqphot/MessageHandling.hpp"
#include "sqphot/QOREInterface.hpp"
#include "sqphot/SparseHbMatrix.hpp"
#include "sqphot/SqpAlgorithm.hpp"
#include "sqphot/Vector.hpp"
#include "sqphot/qpOASESInterface.hpp"
#include <memory>

using namespace SQPhotstart;
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
  shared_ptr<SparseHbMatrix> csc_matrix = make_shared<SparseHbMatrix>(
      csr_matrix->get_num_rows(), csr_matrix->get_num_columns(),
      is_compressed_row);
  csc_matrix->copy_from_dense_matrix(dense_matrix, csr_matrix->get_num_rows(),
                                     csr_matrix->get_num_columns(),
                                     row_oriented, is_compressed_row);
  delete[] dense_matrix;

  return csc_matrix;
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
  else {
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
  }

  shared_ptr<Vector> x_qore = make_shared<Vector>(nCon + nVar);
  shared_ptr<Vector> y_qore_constr = make_shared<Vector>(nCon);
  shared_ptr<Vector> y_qore_bounds = make_shared<Vector>(nVar);

  shared_ptr<Vector> x_qpOASES = make_shared<Vector>(nVar);
  shared_ptr<Vector> y_qpOASES_constr = make_shared<Vector>(nCon);
  shared_ptr<Vector> y_qpOASES_bounds = make_shared<Vector>(nVar);

  shared_ptr<Vector> lb_qore = make_shared<Vector>(nCon + nVar);
  shared_ptr<Vector> ub_qore = make_shared<Vector>(nCon + nVar);
  shared_ptr<Vector> g = make_shared<Vector>(nVar);

  shared_ptr<SparseHbMatrix> A_qore =
      make_shared<SparseHbMatrix>(Annz, nCon, nVar, true);
  shared_ptr<SparseHbMatrix> H_qore =
      make_shared<SparseHbMatrix>(Hnnz, nVar, nVar, true);

  qp_int A_ir_qore[nCon + 1]; // row index
  qp_int A_jc_qore[Annz];     // column index
  double A_val_qore[Annz];    // matrix value

  qp_int H_ir_qore[nVar + 1]; // row index
  qp_int H_jc_qore[Hnnz];     // column index
  double H_val_qore[Hnnz];    // matrix value

  for (i = 0; i < nCon + nVar; i++) {
    if (fgets(buffer, 100, file) != NULL) {
      sscanf(buffer, "%lf", &tmp_double);
      lb_qore->set_value(i, tmp_double);
    }
  }

  for (i = 0; i < nCon + nVar; i++) {
    if (fgets(buffer, 100, file) != NULL) {
      sscanf(buffer, "%lf", &tmp_double);
      ub_qore->set_value(i, tmp_double);
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
      A_qore->set_row_index_at_entry(i, tmp_int);
    }
  }

  for (i = 0; i < Annz; i++) {
    if (fgets(buffer, 100, file) != NULL) {
      sscanf(buffer, "%d", &tmp_int);
      A_qore->set_column_index_at_entry(i, tmp_int);
    }
  }

  for (i = 0; i < Annz; i++) {
    if (fgets(buffer, 100, file) != NULL) {
      sscanf(buffer, "%lf", &tmp_double);
      A_qore->set_value_at_entry(i, tmp_double);
    }
  }

  for (i = 0; i < nVar + 1; i++) {
    if (fgets(buffer, 100, file) != NULL) {
      sscanf(buffer, "%d", &tmp_int);
      H_qore->set_row_index_at_entry(i, tmp_int);
    }
  }

  for (i = 0; i < Hnnz; i++) {
    if (fgets(buffer, 100, file) != NULL) {
      sscanf(buffer, "%d", &tmp_int);
      H_qore->set_column_index_at_entry(i, tmp_int);
    }
  }

  for (i = 0; i < Hnnz; i++) {
    if (fgets(buffer, 100, file) != NULL) {
      sscanf(buffer, "%lf", &tmp_double);
      H_qore->set_value_at_entry(i, tmp_double);
    }
  }

  A_qore->print("A");
  H_qore->print_full("H");

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

  shared_ptr<QOREInterface> qore_inferface =
      make_shared<QOREInterface>(H_qore, A_qore, g, lb_qore, ub_qore, options);
  try {
    qore_inferface->optimize_qp(stats_qore);
  } catch (...) {
  }
  switch (qore_inferface->get_status()) {
    case QPERROR_EXCEED_MAX_ITER:
      exitflag_qore = "EXCEED_ITER_LIMIT";
      break;
    case QPSOLVER_OPTIMAL:
      exitflag_qore = "OPTIMAL";
      break;
    case QPERROR_INFEASIBLE:
      exitflag_qore = "INFEASIBLE";
      break;
    case QPERROR_UNBOUNDED:
      exitflag_qore = "UNBOUNDED";
      break;
    default:
      exitflag_qore = "UNKNOWN";
  }
  qore_inferface->test_optimality();
  x_qore->copy_vector(qore_inferface->get_primal_solution());
  y_qore_constr->copy_vector(qore_inferface->get_constraints_multipliers());
  y_qore_bounds->copy_vector(qore_inferface->get_bounds_multipliers());
  double obj_qore = qore_inferface->get_obj_value();

  ///////////////////////////////////////////////////////////
  //                     QPOASES                           //
  ///////////////////////////////////////////////////////////
  shared_ptr<Statistics> stats_qpOASES = make_shared<Statistics>();
  auto A_qpOASES = convert_csr_to_csc(A_qore);
  auto H_qpOASES = convert_csr_to_csc(H_qore);

  A_qpOASES->print("A_qpOASES");

  auto lb_qpOASES = make_shared<Vector>(nVar);
  auto ub_qpOASES = make_shared<Vector>(nVar);
  auto lbA_qpOASES = make_shared<Vector>(nCon);
  auto ubA_qpOASES = make_shared<Vector>(nCon);

  lb_qpOASES->copy_values(lb_qore->get_values());
  ub_qpOASES->copy_values(ub_qore->get_values());
  lbA_qpOASES->copy_values(lb_qore->get_values() + nVar);
  ubA_qpOASES->copy_values(ub_qore->get_values() + nVar);

  shared_ptr<qpOASESInterface> qpoases_interface =
      make_shared<qpOASESInterface>(H_qpOASES, A_qpOASES, g, lb_qore, ub_qore,
                                    lbA_qpOASES, ubA_qpOASES, options);
  try {
    qpoases_interface->optimize_qp(stats_qpOASES);
  } catch (...) {
  }
  x_qpOASES->copy_vector(qpoases_interface->get_primal_solution());
  y_qpOASES_constr->copy_vector(
      qpoases_interface->get_constraints_multipliers());
  y_qpOASES_bounds->copy_vector(qpoases_interface->get_bounds_multipliers());

  qpoases_interface->test_optimality();
  double obj_qpOASES = qpoases_interface->get_obj_value();
  switch (qpoases_interface->get_status()) {
    case QPERROR_EXCEED_MAX_ITER:
      exitflag_qpOASES = "EXCEED_ITER_LIMIT";
      break;
    case QPSOLVER_OPTIMAL:
      exitflag_qpOASES = "OPTIMAL";
      break;
    case QPERROR_INFEASIBLE:
      exitflag_qpOASES = "INFEASIBLE";
      break;
    case QPERROR_UNBOUNDED:
      exitflag_qpOASES = "UNBOUNDED";
      break;
    case QPERROR_NOTINITIALISED:
      exitflag_qpOASES = "NOTINITIALISED";
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
    default:
      exitflag_qpOASES = "UNKNOWN";
  }

  //////////////////////////////////////////////////////////
  //                     PRINT OUT RESULTS                 /"
  ///////////////////////////////////////////////////////////e
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
         qore_inferface->get_optimality_status().KKT_error,
         qpoases_interface->get_optimality_status().KKT_error);
  printf("%30s    %23.16e    %23.16e\n", "Primal Feasibility Violation",
         qore_inferface->get_optimality_status().primal_violation,
         qpoases_interface->get_optimality_status().primal_violation);
  printf("%30s    %23.16e    %23.16e\n", "Dual Feasibility  Violation",
         qore_inferface->get_optimality_status().dual_violation,
         qpoases_interface->get_optimality_status().dual_violation);
  printf("%30s    %23.16e    %23.16e\n", "Stationarity Violation",
         qore_inferface->get_optimality_status().stationarity_violation,
         qpoases_interface->get_optimality_status().stationarity_violation);
  printf("%30s    %23.16e    %23.16e\n", "Complemtarity Violation",
         qore_inferface->get_optimality_status().compl_violation,
         qpoases_interface->get_optimality_status().compl_violation);
  printf(DOUBLE_LONG_DIVIDER);

  delete pname;
  fclose(file);
  return 0;
}
