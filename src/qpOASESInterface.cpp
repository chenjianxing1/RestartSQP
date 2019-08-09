
#include <sqphot/qpOASESInterface.hpp>


namespace SQPhotstart {
/**
 * @brief Allocate memory for the class members
 * @param nlp_index_info  the struct that stores simple nlp dimension info
 * @param qptype is the problem to be solved QP or LP or SOC?
 * @return
 */
void qpOASESInterface::allocate(Index_info nlp_index_info, QPType qptype) {

    int nVar_QP = 2 * nlp_index_info.nCon + nlp_index_info.nVar;
    int nCon_QP = nlp_index_info.nCon;

    lbA_ = make_shared<Vector>(nCon_QP);
    ubA_ = make_shared<Vector>(nCon_QP);
    lb_ = make_shared<Vector>(nVar_QP);
    ub_ = make_shared<Vector>(nVar_QP);
    g_ = make_shared<Vector>(nVar_QP);
    A_ = make_shared<qpOASESSparseMat>(
             nlp_index_info.nnz_jac_g + 2 * nlp_index_info.nCon, nCon_QP, nVar_QP);

    if (qptype != LP) {
        H_ = make_shared<qpOASESSparseMat>(nVar_QP, nVar_QP, true);
    }

    //FIXME: the qpOASES does not accept any extra input
    solver_ = std::make_shared<qpOASES::SQProblem>((qpOASES::int_t) nVar_QP,
              (qpOASES::int_t) nCon_QP);
}


/**
 * @brief Constructor which also initializes the qpOASES SQProblem objects
 * @param nlp_index_info the struct that stores simple nlp dimension info
 * @param qptype  is the problem to be solved QP or LP or SOC?
 */
qpOASESInterface::qpOASESInterface(Index_info nlp_index_info, QPType qptype) :
    status_(UNSOLVED) {

    allocate(nlp_index_info, qptype);
}


/**Default destructor*/
qpOASESInterface::~qpOASESInterface() = default;


/**
 * @brief This function solves the QP problem specified in the data, with given
 * options.
 * After the QP being solved, it updates the stats, adding the iteration
 * number used to solve the QP to the qp_iter in object stats
 */
void
qpOASESInterface::optimizeQP(shared_ptr<Stats> stats, shared_ptr<Options> options) {

    qpOASES::int_t nWSR = options->qp_maxiter;

    if (!firstQPsolved_) {//if haven't solve any QP before then initialize the first QP
        if (DEBUG) {
            if (CHECK_QP_INFEASIBILITY) {
                A_qpOASES_->print("A_");
                H_qpOASES_->print("H");
                g_->print("g");
                lb_->print("lb_");
                ub_->print("ub_");
                lbA_->print("lbA_");
                ubA_->print("ubA_");
            }
        }

        setQP_options(options);

        solver_->init(H_qpOASES_.get(), g_->values(), A_qpOASES_.get(), lb_->values(),
                      ub_->values(), lbA_->values(), ubA_->values(), nWSR);

        if (solver_->isSolved())
            firstQPsolved_ = true;
        else {
            obtain_status();
            handler_error(QP, stats, options);
        }
    }
    else {
        get_Matrix_change_status();
        if (new_QP_matrix_status_ == UNDEFINED) {
            assert(old_QP_matrix_status_ != UNDEFINED);
            if (old_QP_matrix_status_ == FIXED)
                solver_->hotstart(g_->values(), lb_->values(), ub_->values(),
                                  lbA_->values(),
                                  ubA_->values(), nWSR);
            else {
                solver_->hotstart(H_qpOASES_.get(), g_->values(), A_qpOASES_.get(),
                                  lb_->values(), ub_->values(), lbA_->values(),
                                  ubA_->values(), nWSR);

            }
        }
        else {
            if (new_QP_matrix_status_ == FIXED && old_QP_matrix_status_ == FIXED) {
                solver_->hotstart(g_->values(), lb_->values(), ub_->values(),
                                  lbA_->values(),
                                  ubA_->values(), nWSR);
            }
            else if (new_QP_matrix_status_ == VARIED &&
                     old_QP_matrix_status_ == VARIED) {
                solver_->hotstart(H_qpOASES_.get(), g_->values(), A_qpOASES_.get(),
                                  lb_->values(), ub_->values(), lbA_->values(),
                                  ubA_->values(), nWSR);
            }
            else if (new_QP_matrix_status_ != old_QP_matrix_status_) {
                solver_->init(H_qpOASES_.get(), g_->values(), A_qpOASES_.get(),
                              lb_->values(),
                              ub_->values(), lbA_->values(), ubA_->values(), nWSR);
                new_QP_matrix_status_ = old_QP_matrix_status_ = UNDEFINED;
            }
        }
    }

    reset_flags();

    stats->qp_iter_addValue((int) nWSR);

    if (!solver_->isSolved()) {
        obtain_status();
        handler_error(QP, stats, options);
    }

}


void qpOASESInterface::optimizeLP(shared_ptr<Stats> stats, shared_ptr<Options>
                                  options) {

    qpOASES::int_t nWSR = options->lp_maxiter;//TODO modify it
    if (!firstQPsolved_) {
        setQP_options(options);
        solver_->init(0, g_->values(), A_qpOASES_.get(), lb_->values(),
                      ub_->values(), lbA_->values(), ubA_->values(), nWSR);
        if (solver_->isSolved()) {
            firstQPsolved_ = true;
        }
        else {
            obtain_status();
            handler_error(LP, stats, options);
        }
    }
    else {
        get_Matrix_change_status();
        if (new_QP_matrix_status_ == UNDEFINED) {
            if (old_QP_matrix_status_ == FIXED)
                solver_->hotstart(g_->values(), lb_->values(), ub_->values(),
                                  lbA_->values(),
                                  ubA_->values(), nWSR);
            else {
                solver_->hotstart(0, g_->values(), A_qpOASES_.get(),
                                  lb_->values(), ub_->values(), lbA_->values(),
                                  ubA_->values(), nWSR);
            }
        }

        else {
            if (new_QP_matrix_status_ == FIXED && old_QP_matrix_status_ == FIXED) {
                solver_->hotstart(g_->values(), lb_->values(), ub_->values(),
                                  lbA_->values(),
                                  ubA_->values(), nWSR);
            }
            else if (new_QP_matrix_status_ == VARIED &&
                     old_QP_matrix_status_ == VARIED) {
                solver_->hotstart(0, g_->values(), A_qpOASES_.get(),
                                  lb_->values(), ub_->values(), lbA_->values(),
                                  ubA_->values(), nWSR);
            }
            else if (new_QP_matrix_status_ != old_QP_matrix_status_) {
                solver_->init(0, g_->values(), A_qpOASES_.get(), lb_->values(),
                              ub_->values(), lbA_->values(), ubA_->values(), nWSR);
                new_QP_matrix_status_ = old_QP_matrix_status_ = UNDEFINED;
            }
        }
        reset_flags();
        if (!solver_->isSolved()) {
            obtain_status();
            handler_error(LP, stats, options);
        }
        stats->qp_iter_addValue((int) nWSR);
    }
}


/**
* @brief copy the multipliers of the QP to the input pointer
*
* @param y_k   a pointer to an array with allocated memory equals to
* sizeof(double)*(num_variable+num_constraint)
*/
inline void qpOASESInterface::get_multipliers(double* y_k) {

    solver_->getDualSolution(y_k);
}


/**
 * @brief copy the optimal solution of the QP to the input pointer
 *
 * @param x_optimal a pointer to an empty array with allocated memory equals to
 * sizeof(double)*number_variables
 *
 */
inline void qpOASESInterface::get_optimal_solution(double* p_k) {

    solver_->getPrimalSolution(p_k);
}


/**
 *@brief get the objective value from the QP solvers
 *
 * @return the objective function value of the QP problem
 */


inline double qpOASESInterface::get_obj_value() {

    return (double) (solver_->getObjVal());
}


void qpOASESInterface::obtain_status() {

    qpOASES::QProblemStatus finalStatus = solver_->getStatus();

    if (solver_->isInfeasible()) {
        status_ = QP_INFEASIBLE;
    }
    else if (solver_->isUnbounded()) {
        status_ = QP_UNBOUNDED;
    }
    else
        switch (finalStatus) {
        case qpOASES::QPS_NOTINITIALISED:
            status_ = QP_NOTINITIALISED;
            break;
        case qpOASES::QPS_PREPARINGAUXILIARYQP:
            status_ = QP_PREPARINGAUXILIARYQP;
            break;
        case qpOASES::QPS_AUXILIARYQPSOLVED:
            status_ = QP_AUXILIARYQPSOLVED;
            break;
        case qpOASES::QPS_PERFORMINGHOMOTOPY:
            status_ = QP_PERFORMINGHOMOTOPY;
            break;
        case qpOASES::QPS_HOMOTOPYQPSOLVED:
            status_ = QP_HOMOTOPYQPSOLVED;
            break;
        }
}


QPReturnType qpOASESInterface::get_status() {

    return status_;
}


//@}
void qpOASESInterface::set_lb(int location, double value) {

    if (firstQPsolved_ && !data_change_flags_.Update_bounds)
        data_change_flags_.Update_bounds = true;
    lb_->setValueAt(location, value);
}


void qpOASESInterface::set_ub(int location, double value) {

    if (firstQPsolved_ && !data_change_flags_.Update_bounds)
        data_change_flags_.Update_bounds = true;
    ub_->setValueAt(location, value);
}


void qpOASESInterface::set_lbA(int location, double value) {

    if (firstQPsolved_ && !data_change_flags_.Update_bounds)
        data_change_flags_.Update_bounds = true;
    lbA_->setValueAt(location, value);
}


void qpOASESInterface::set_ubA(int location, double value) {

    if (firstQPsolved_ && !data_change_flags_.Update_bounds)
        data_change_flags_.Update_bounds = true;
    ubA_->setValueAt(location, value);
}


void qpOASESInterface::set_g(int location, double value) {

    if (firstQPsolved_ && !data_change_flags_.Update_g)
        data_change_flags_.Update_g = true;
    g_->setValueAt(location, value);
}


void qpOASESInterface::set_H_structure(shared_ptr<const SpTripletMat> rhs) {

    H_->setStructure(rhs);//TODO: move to somewhere else?
    H_qpOASES_ = std::make_shared<qpOASES::SymSparseMat>(H_->RowNum(),
                 H_->ColNum(),
                 H_->RowIndex(),
                 H_->ColIndex(),
                 H_->MatVal());
}


void qpOASESInterface::set_H_values(shared_ptr<const SpTripletMat> rhs) {

    if (firstQPsolved_ && !data_change_flags_.Update_H) {
        data_change_flags_.Update_H = true;
    }
    H_->setMatVal(rhs);
    H_qpOASES_->setVal(H_->MatVal());
    H_qpOASES_->createDiagInfo();
}


void qpOASESInterface::set_A_structure(shared_ptr<const SpTripletMat> rhs,
                                       Identity2Info I_info) {

    A_->setStructure(rhs, I_info);
    A_qpOASES_ = std::make_shared<qpOASES::SparseMatrix>(A_->RowNum(),
                 A_->ColNum(),
                 A_->RowIndex(),
                 A_->ColIndex(),
                 A_->MatVal());
}


void qpOASESInterface::set_A_values(
    shared_ptr<const SQPhotstart::SpTripletMat> rhs, Identity2Info I_info) {

    if (firstQPsolved_ && !data_change_flags_.Update_A) {
        data_change_flags_.Update_A = true;
    }
    A_->setMatVal(rhs, I_info);
    A_qpOASES_->setVal(A_->MatVal());
}


void qpOASESInterface::set_ub(shared_ptr<const Vector> rhs) {

    if (firstQPsolved_ && !data_change_flags_.Update_bounds)
        data_change_flags_.Update_bounds = true;
    ub_->copy_vector(rhs->values());
}


void qpOASESInterface::set_lb(shared_ptr<const Vector> rhs) {

    if (firstQPsolved_ && !data_change_flags_.Update_bounds)
        data_change_flags_.Update_bounds = true;
    lb_->copy_vector(rhs->values());

}


void qpOASESInterface::set_lbA(shared_ptr<const Vector> rhs) {

    if (firstQPsolved_ && !data_change_flags_.Update_bounds)
        data_change_flags_.Update_bounds = true;
    lbA_->copy_vector(rhs->values());
}


void qpOASESInterface::set_ubA(shared_ptr<const Vector> rhs) {

    if (firstQPsolved_ && !data_change_flags_.Update_bounds)
        data_change_flags_.Update_bounds = true;
    ubA_->copy_vector(rhs->values());

}


void qpOASESInterface::set_g(shared_ptr<const Vector> rhs) {

    if (firstQPsolved_ && !data_change_flags_.Update_g)
        data_change_flags_.Update_g = true;
    g_->copy_vector(rhs->values());
}


#if DEBUG
#if GET_QPOASES_MEMBERS
const shared_ptr<Vector>& qpOASESInterface::getLb() const {
    return lb_;
}

const shared_ptr<Vector>& qpOASESInterface::getUb() const {
    return ub_;
}

const shared_ptr<Vector>& qpOASESInterface::getLbA() const {
    return lbA_;
}

const shared_ptr<Vector>& qpOASESInterface::getUbA() const {
    return ubA_;
}

const shared_ptr<Vector>& qpOASESInterface::getG() const {
    return g_;
}
#endif
#endif


void qpOASESInterface::reset_flags() {

    data_change_flags_.Update_A = false;
    data_change_flags_.Update_delta = false;
    data_change_flags_.Update_H = false;
    data_change_flags_.Update_penalty = false;
    data_change_flags_.Update_g = false;
    data_change_flags_.Update_bounds = false;
}


void qpOASESInterface::handler_error(QPType qptype, shared_ptr<Stats> stats,
                                     shared_ptr<Options> options) {

    if (qptype == LP) {
        if (true) {
            qpOASES::int_t nWSR = options->lp_maxiter;//TODO modify it
            solver_->init(0, g_->values(), A_qpOASES_.get(),
                          lb_->values(), ub_->values(), lbA_->values(),
                          ubA_->values(), nWSR);
            obtain_status();
            old_QP_matrix_status_ = new_QP_matrix_status_ = UNDEFINED;
            stats->qp_iter_addValue((int) nWSR);

            if (!solver_->isSolved()) {
                if (DEBUG) {
                    if (PRINT_OUT_QP_WITH_ERROR) {
                        WriteQPDataToFile("test");
                    }

                    THROW_EXCEPTION(LP_NOT_OPTIMAL,
                                    "the QP problem didn't solved to optimality\n")
                }
            }
            else {
                if (PRINT_OUT_QP_WITH_ERROR) {
                    WriteQPDataToFile("test");
                }
                THROW_EXCEPTION(LP_NOT_OPTIMAL,
                                "the QP problem didn't solved to optimality\n")
            }

        }
    }
    else {
        if (true) {
            qpOASES::int_t nWSR = options->lp_maxiter;//TODO modify it
            solver_->init(H_qpOASES_.get(), g_->values(), A_qpOASES_.get(),
                          lb_->values(), ub_->values(), lbA_->values(),
                          ubA_->values(), nWSR);
            obtain_status();
            old_QP_matrix_status_ = new_QP_matrix_status_ = UNDEFINED;
            stats->qp_iter_addValue((int) nWSR);
            if (!solver_->isSolved()) {
                if (DEBUG) {
                    if (PRINT_OUT_QP_WITH_ERROR) {
                        WriteQPDataToFile("test");
                    }
                }
                THROW_EXCEPTION(QP_NOT_OPTIMAL,
                                "the QP problem didn't solved to optimality\n")

            }
        }
        else {
            if (DEBUG)
                if (PRINT_OUT_QP_WITH_ERROR) {
                    WriteQPDataToFile("test");
                }
            THROW_EXCEPTION(QP_NOT_OPTIMAL,
                            "the QP problem didn't solved to optimality\n")
        }
    }
}


void qpOASESInterface::setQP_options(shared_ptr<Options> options) {

    qpOASES::Options qp_options;
    qp_options.setToReliable();
    //setup the printlevel of q
    switch (options->qpPrintLevel) {
    case 0:
        qp_options.printLevel = qpOASES::PL_NONE;
        break;
    case 1:
        qp_options.printLevel = qpOASES::PL_TABULAR;
        break;
    case 2:
        qp_options.printLevel = qpOASES::PL_LOW;
        break;
    case 3:
        qp_options.printLevel = qpOASES::PL_MEDIUM;
        break;
    case 4:
        qp_options.printLevel = qpOASES::PL_HIGH;
        break;
    case -2:
        qp_options.printLevel = qpOASES::PL_DEBUG_ITER;
        break;
    }

    solver_->setOptions(qp_options);
}


void qpOASESInterface::WriteQPDataToFile(const char* const filename) {

#if DEBUG
    FILE* file;
    file = fopen(filename,"w");
    lb_->write_to_file(file,"lb");
    ub_->write_to_file(file,"ub");
    lbA_->write_to_file(file,"lbA");
    ubA_->write_to_file(file,"ubA");
    g_->write_to_file(file,"g");
    A_->write_to_file(file,"A");
    H_->write_to_file(file,"H");
#endif

}


void qpOASESInterface::get_Matrix_change_status() {

    if (old_QP_matrix_status_ == UNDEFINED) {
        old_QP_matrix_status_ = data_change_flags_.Update_A ||
                                data_change_flags_.Update_H ? VARIED
                                : FIXED;
    }
    else {
        if (new_QP_matrix_status_ != UNDEFINED)
            old_QP_matrix_status_ = new_QP_matrix_status_;
        new_QP_matrix_status_ = data_change_flags_.Update_A ||
                                data_change_flags_.Update_H ? VARIED
                                : FIXED;
    }


}

}//SQPHOTSTART

