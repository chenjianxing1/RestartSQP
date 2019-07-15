#include <sqphot/QPsolverInterface.hpp>

namespace SQPhotstart {
    
    bool qpOASESInterface::allocate(Index_info nlp_index_info, QPType qptype) {
        int nVar_QP = 2 * nlp_index_info.nCon + nlp_index_info.nVar;
        int nCon_QP = nlp_index_info.nCon;
        lbA_ = make_shared<Vector>(nCon_QP);
        ubA_ = make_shared<Vector>(nCon_QP);
        lb_ = make_shared<Vector>(nVar_QP);
        ub_ = make_shared<Vector>(nVar_QP);
        g_ = make_shared<Vector>(nVar_QP);
        A_ = make_shared<qpOASESSparseMat>(nlp_index_info.nnz_jac_g + 2 * nlp_index_info.nCon, nCon_QP, nVar_QP);
        if(qptype!=LP){
            H_ = make_shared<qpOASESSparseMat>(nlp_index_info.nnz_h_lag, nVar_QP, nVar_QP);
        }//FIXME: the qpOASES does not accept any extra input
        qp_ = std::make_shared<qpOASES::SQProblem>((qpOASES::int_t) nVar_QP, (qpOASES::int_t) nCon_QP);
        return true;
    }
    
    /**
     * @name Constructor which also initializes the qpOASES SQProblem objects
     * @param nlp_index_info the number of variables in QP problem
     * @param nCon_QP the number of constraints in QP problem (the number of rows of A)
     */
    qpOASESInterface::qpOASESInterface(Index_info nlp_index_info, QPType qptype) {
        allocate(nlp_index_info, qptype);
    }
    
    /**Default destructor*/
    qpOASESInterface::~qpOASESInterface() = default;
    
    /**
     * @name This function solves the QP problem specified in the data, with given options. After the QP
     * is solved, it updates the stats, adding the iteration number used to solve the QP to the qp_iter
     * in object stats
     */
    bool qpOASESInterface::optimizeQP(shared_ptr<Stats> stats, shared_ptr<Options> options) {
        H_qpOASES_ = std::make_shared<qpOASES::SymSparseMat>(H_->RowNum(), H_->ColNum(),H_->RowIndex(), H_->ColIndex(), H_->MatVal());
        A_qpOASES_ = std::make_shared<qpOASES::SparseMatrix>(A_->RowNum(), A_->ColNum(),A_->RowIndex(), A_->ColIndex(), A_->MatVal());
        H_qpOASES_->createDiagInfo();
        
        qpOASES::int_t nWSR = options->qp_maxiter;//TODO modify it
        if (firstQPsolved==false) {//if haven't solve any QP before then initialize the first qp
            qpOASES::Options qp_options;
            if (options->qpPrintLevel == 0)//else use the default print level in qpOASES
                qp_options.printLevel = qpOASES::PL_NONE;
            qp_->setOptions(qp_options);
            qp_->init(H_qpOASES_.get(), g_->values(), A_qpOASES_.get(), lb_->values(), ub_->values(),
                      lbA_->values(), ubA_->values(), nWSR, 0);
            if(qp_->isSolved())
                firstQPsolved = true;
        } else
            qp_->hotstart(H_qpOASES_.get(), g_->values(), A_qpOASES_.get(), lb_->values(), ub_->values(),
                          lbA_->values(), ubA_->values(), nWSR, 0);
        stats->qp_iter_addValue((int) nWSR);
        return true;
    }
    
    bool qpOASESInterface::get_multipliers(double* y_k) {
        qp_->getDualSolution(y_k);
        return true;
    }
    
    bool qpOASESInterface::get_optimal_solution(double* p_k) {
        qp_->getPrimalSolution(p_k);
        return true;
    }
    
    double qpOASESInterface::get_obj_value() {
        return (double) (qp_->getObjVal());
    }
    
    
    shared_ptr<Vector> &qpOASESInterface::getLb() {
        return lb_;
    }
    
    shared_ptr<Vector> &qpOASESInterface::getUb() {
        return ub_;
    }
    
    shared_ptr<Vector> &qpOASESInterface::getLbA() {
        return lbA_;
    }
    
    shared_ptr<Vector> &qpOASESInterface::getUbA() {
        return ubA_;
    }
    
    shared_ptr<Vector> &qpOASESInterface::getG() {
        return g_;
    }
    
    shared_ptr<qpOASESSparseMat> &qpOASESInterface::getH(){
        return H_;
    }
    
    shared_ptr<qpOASESSparseMat> &qpOASESInterface::getA(){
        return A_;
    }
    
}//SQPHOTSTART

