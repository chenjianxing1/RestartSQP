#include <sqphot/QPsolverInterface.hpp>

namespace SQPhotstart {
/**
 * @brief Allocate memory for the class members
 * @param nlp_index_info  the struct that stores simple nlp dimension info
 * @param qptype is the problem to be solved QP or LP or SOC?
 * @return
 */
bool qpOASESInterface::allocate(Index_info nlp_index_info, QPType qptype) {
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
    qp_ = std::make_shared<qpOASES::SQProblem>((qpOASES::int_t) nVar_QP,
            (qpOASES::int_t) nCon_QP);
    return true;
}

/**
 * @brief Constructor which also initializes the qpOASES SQProblem objects
 * @param nlp_index_info the struct that stores simple nlp dimension info
 * @param qptype  is the problem to be solved QP or LP or SOC?
 */
qpOASESInterface::qpOASESInterface(Index_info nlp_index_info, QPType qptype):
    status_(UNSOLVED)
{
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
bool
qpOASESInterface::optimizeQP(shared_ptr<Stats> stats, shared_ptr<Options> options) {
    H_qpOASES_ = std::make_shared<qpOASES::SymSparseMat>(H_->RowNum(), H_->ColNum(),
                 H_->RowIndex(),
                 H_->ColIndex(),
                 H_->MatVal());
    A_qpOASES_ = std::make_shared<qpOASES::SparseMatrix>(A_->RowNum(), A_->ColNum(),
                 A_->RowIndex(),
                 A_->ColIndex(),
                 A_->MatVal());

    H_qpOASES_->createDiagInfo();
    if(DEBUG) {
        if(CHECK_QP_INFEASIBILITY) {
            A_qpOASES_->print("A_");
            lb_->print("lb_");
            ub_->print("ub_");
            lbA_->print("lbA_");
            ubA_->print("ubA_");
        }
    }
    qpOASES::int_t nWSR = options->qp_maxiter;

    if (!firstQPsolved) {//if haven't solve any QP before then initialize the first QP
        qpOASES::Options qp_options;

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


        qp_->setOptions(qp_options);

        qp_->init(H_qpOASES_.get(), g_->values(), A_qpOASES_.get(), lb_->values(),
                  ub_->values(), lbA_->values(), ubA_->values(), nWSR, 0);
        if (qp_->isSolved())
            firstQPsolved = true;
        else {
            get_status();

            //  THROW_EXCEPTION(QP_NOT_OPTIMAL,
            //                  "the QP problem didn't solved to optimality\n")

        }
    } else
        qp_->hotstart(H_qpOASES_.get(), g_->values(), A_qpOASES_.get(),
                      lb_->values(), ub_->values(), lbA_->values(), ubA_->values(),
                      nWSR, 0);
    if (!qp_->isSolved()) {
        get_status();
        //  THROW_EXCEPTION(QP_NOT_OPTIMAL,
        //                  "the LP problem didn't solved to optimality\n")
    }
    stats->qp_iter_addValue((int) nWSR);
    return true;
}

bool qpOASESInterface::optimizeLP(shared_ptr<Stats> stats, shared_ptr<Options>
                                  options) {

    A_qpOASES_ = std::make_shared<qpOASES::SparseMatrix>(A_->RowNum(), A_->ColNum(),
                 A_->RowIndex(),
                 A_->ColIndex(),
                 A_->MatVal());

    qpOASES::int_t nWSR = options->lp_maxiter;//TODO modify it

    if (!firstQPsolved) {//if haven't solve any LP before then initialize the
        //  first qp
        qpOASES::Options qp_options;
        //setup the printlevel of qpOASES
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

        qp_->setOptions(qp_options);
        qp_->init(0, g_->values(), A_qpOASES_.get(), lb_->values(),
                  ub_->values(), lbA_->values(), ubA_->values(), nWSR, 0);
        if (qp_->isSolved())
            firstQPsolved = true;
        else {
            get_status();
            //    THROW_EXCEPTION(QP_NOT_OPTIMAL,
            //                    "the LP problem didn't solved to optimality\n")
        }

    } else
        qp_->hotstart(0, g_->values(), A_qpOASES_.get(),
                      lb_->values(), ub_->values(), lbA_->values(), ubA_->values(),
                      nWSR, 0);
    if (!qp_->isSolved()) {
        get_status();
        //THROW_EXCEPTION(QP_NOT_OPTIMAL,
        //                 "the LP problem didn't solved to optimality\n")
    }
    stats->qp_iter_addValue((int) nWSR);

    return true;

}

/**
* @brief copy the multipliers of the QP to the input pointer
*
* @param y_k   a pointer to an array with allocated memory equals to
* sizeof(double)*(num_variable+num_constraint)
*/
inline bool qpOASESInterface::get_multipliers(double* y_k) {
    qp_->getDualSolution(y_k);
    return true;
}

/**
 * @brief copy the optimal solution of the QP to the input pointer
 *
 * @param x_optimal a pointer to an empty array with allocated memory equals to
 * sizeof(double)*number_variables
 *
 */
inline bool qpOASESInterface::get_optimal_solution(double* p_k) {
    qp_->getPrimalSolution(p_k);
    return true;
}


/**
 *@brief get the objective value from the QP solvers
 *
 * @return the objective function value of the QP problem
 */


inline double qpOASESInterface::get_obj_value() {
    return (double) (qp_->getObjVal());
}

/** Getters, extract private member information*/
//@{
shared_ptr<Vector>& qpOASESInterface::getLb() {
    return lb_;
}

shared_ptr<Vector>& qpOASESInterface::getUb() {
    return ub_;
}

shared_ptr<Vector>& qpOASESInterface::getLbA() {
    return lbA_;
}

shared_ptr<Vector>& qpOASESInterface::getUbA() {
    return ubA_;
}

shared_ptr<Vector>& qpOASESInterface::getG() {
    return g_;
}

shared_ptr<qpOASESSparseMat>& qpOASESInterface::getH() {
    return H_;
}

shared_ptr<qpOASESSparseMat>& qpOASESInterface::getA() {
    return A_;
}

bool qpOASESInterface::get_status() {
    qpOASES::QProblemStatus finalStatus = qp_->getStatus();

    //switch(finalStatus) {
    //case qpOASES::QPS_NOTINITIALISED:
    //    printf("QP_NOTINITIALISED\n");
    //    break;
    //case qpOASES::QPS_PREPARINGAUXILIARYQP:
    //    printf("QPS_PREPARINGAUXILIARYQP\n");
    //    break;
    //case qpOASES::QPS_AUXILIARYQPSOLVED:
    //    printf( "QPS_AUXILIARYQPSOLVED\n");
    //    break;
    //case qpOASES::QPS_PERFORMINGHOMOTOPY:
    //    printf( "QPS_PERFORMINGHOMOTOPY\n");
    //    break;
    //case qpOASES::QPS_HOMOTOPYQPSOLVED:
    //    printf("QPS_HOMOTOPYQPSOLVED\n");
    //    break;

    //}

    if (finalStatus == qpOASES::QPS_NOTINITIALISED)
        status_ = QP_NOTINITIALISED;
    else if (qp_->isInfeasible()) {
        status_= QP_INFEASIBLE;
    }
    else if (qp_->isUnbounded()) {
        status_ = QP_UNBOUNDED;
    }
    else
        status_ = QP_UNKNOWN_ERROR;
    return true;
}


QPReturnType qpOASESInterface::status() {
    return status_;
}
//@}


}//SQPHOTSTART

