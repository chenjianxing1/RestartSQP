#include <math.h>
#include <stdio.h>
#include <string.h>

extern "C" {
#include <qpsolver.h>
}

#include "sqphot/MessageHandling.hpp"
#include "sqphot/Options.hpp"
#include "sqphot/QOREInterface.hpp"
#include "sqphot/SpHbMat.hpp"
#include "sqphot/Vector.hpp"
#include "sqphot/qpOASESInterface.hpp"
#include <memory>

using namespace SQPhotstart;
using namespace std;

shared_ptr<SpHbMat> convert_csr_to_csc(shared_ptr<const SpHbMat> csr_matrix)
{

  double* dense_matrix;
  dense_matrix = new double[csr_matrix->RowNum() *
                            csr_matrix->ColNum()](); // get the dense matrix out

  csr_matrix->get_dense_matrix(dense_matrix);

  shared_ptr<SpHbMat> csc_matrix = make_shared<SpHbMat>(
      dense_matrix, csr_matrix->RowNum(), csr_matrix->ColNum(), true, false);
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

  shared_ptr<SpHbMat> A_qore = make_shared<SpHbMat>(Annz, nCon, nVar, true);
  shared_ptr<SpHbMat> H_qore = make_shared<SpHbMat>(Hnnz, nVar, nVar, true);

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
      A_qore->setRowIndexAt(i, tmp_int);
    }
  }

  for (i = 0; i < Annz; i++) {
    if (fgets(buffer, 100, file) != NULL) {
      sscanf(buffer, "%d", &tmp_int);
      A_qore->setColIndexAt(i, tmp_int);
    }
  }

  for (i = 0; i < Annz; i++) {
    if (fgets(buffer, 100, file) != NULL) {
      sscanf(buffer, "%lf", &tmp_double);
      A_qore->setMatValAt(i, tmp_double);
    }
  }

  for (i = 0; i < nVar + 1; i++) {
    if (fgets(buffer, 100, file) != NULL) {
      sscanf(buffer, "%d", &tmp_int);
      H_qore->setRowIndexAt(i, tmp_int);
    }
  }

  for (i = 0; i < Hnnz; i++) {
    if (fgets(buffer, 100, file) != NULL) {
      sscanf(buffer, "%d", &tmp_int);
      H_qore->setColIndexAt(i, tmp_int);
    }
  }

  for (i = 0; i < Hnnz; i++) {
    if (fgets(buffer, 100, file) != NULL) {
      sscanf(buffer, "%lf", &tmp_double);
      H_qore->setMatValAt(i, tmp_double);
    }
  }

  //    H_qore->print_full("H");

  // read matrix data
  ///////////////////////////////////////////////////////////
  //                     QORE                              //
  ///////////////////////////////////////////////////////////

  //    NLPInfo nlp_info;
  //    nlp_info.nCon = nCon;
  //    nlp_info.nVar = nVar-2*nCon;
  //    nlp_info.nnz_jac_g = Annz;
  //    nlp_info.nnz_h_lag = Hnnz;

  shared_ptr<SQPhotstart::Stats> stats_qore = make_shared<SQPhotstart::Stats>();
  shared_ptr<SQPhotstart::Options> options =
      make_shared<SQPhotstart::Options>();
  options->qpPrintLevel = 0;
  Ipopt::SmartPtr<Ipopt::Journalist> jnlst = new Ipopt::Journalist();

  shared_ptr<QOREInterface> qore_inferface =
      make_shared<QOREInterface>(H_qore, A_qore, g, lb_qore, ub_qore, options);
  try {
    qore_inferface->optimizeQP(stats_qore);
  } catch (...) {
  }
  switch (qore_inferface->get_status()) {
    case QPERROR_EXCEED_MAX_ITER:
      exitflag_qore = "EXCEED_ITER_LIMIT";
      break;
    case QP_OPTIMAL:
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
  x_qore->copy_vector(qore_inferface->get_optimal_solution());
  y_qore_constr->copy_vector(qore_inferface->get_constraints_multipliers());
  y_qore_bounds->copy_vector(qore_inferface->get_bounds_multipliers());
  double obj_qore = qore_inferface->get_obj_value();

  ///////////////////////////////////////////////////////////
  //                     QPOASES                           //
  ///////////////////////////////////////////////////////////
  shared_ptr<SQPhotstart::Stats> stats_qpOASES =
      make_shared<SQPhotstart::Stats>();
  auto A_qpOASES = convert_csr_to_csc(A_qore);
  auto H_qpOASES = convert_csr_to_csc(H_qore);

  auto lb_qpOASES = make_shared<Vector>(nVar);
  auto ub_qpOASES = make_shared<Vector>(nVar);
  auto lbA_qpOASES = make_shared<Vector>(nCon);
  auto ubA_qpOASES = make_shared<Vector>(nCon);

  lb_qpOASES->copy_values(lb_qore->values());
  ub_qpOASES->copy_values(ub_qore->values());
  lbA_qpOASES->copy_values(lb_qore->values() + nVar);
  ubA_qpOASES->copy_values(ub_qore->values() + nVar);

  shared_ptr<qpOASESInterface> qpoases_interface =
      make_shared<qpOASESInterface>(H_qpOASES, A_qpOASES, g, lb_qore, ub_qore,
                                    lbA_qpOASES, ubA_qpOASES, options);
  try {
    qpoases_interface->optimizeQP(stats_qpOASES);
  } catch (...) {
  }
  x_qpOASES->copy_vector(qpoases_interface->get_optimal_solution());
  y_qpOASES_constr->copy_vector(
      qpoases_interface->get_constraints_multipliers());
  y_qpOASES_bounds->copy_vector(qpoases_interface->get_bounds_multipliers());

  qpoases_interface->test_optimality();
  double obj_qpOASES = qpoases_interface->get_obj_value();
  switch (qpoases_interface->get_status()) {
    case QPERROR_EXCEED_MAX_ITER:
      exitflag_qpOASES = "EXCEED_ITER_LIMIT";
      break;
    case QP_OPTIMAL:
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
  printf("%30s    %23d    %23d\n", "Iteration", stats_qore->qp_iter,
         stats_qpOASES->qp_iter);
  printf("%30s    %23.16e    %23.16e\n", "Objective", obj_qore, obj_qpOASES);
  printf("%30s    %23.16e    %23.16e\n", "||x||", x_qore->one_norm(),
         x_qpOASES->one_norm());
  printf("%30s    %23.16e    %23.16e\n", "||y_b||", y_qore_bounds->inf_norm(),
         y_qpOASES_bounds->inf_norm());
  printf("%30s    %23.16e    %23.16e\n", "||y_c||", y_qore_constr->inf_norm(),
         y_qpOASES_constr->inf_norm());
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
