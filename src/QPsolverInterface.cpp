#include <sqphot/QPsolverInterface.hpp>

namespace SQPhotstart {

    bool qpOASESInterface::allocate(Index_info nlp_index_info) {
        int nVar_QP = 2 * nlp_index_info.nCon + nlp_index_info.nVar;
        int nCon_QP = nlp_index_info.nCon;
        lbA_ = make_shared<Vector>(nCon_QP);
        ubA_ = make_shared<Vector>(nCon_QP);
        lb_ = make_shared<Vector>(nVar_QP);
        ub_ = make_shared<Vector>(nVar_QP);
        g_ = make_shared<Vector>(nVar_QP);
        A_ = make_shared<qpOASESSparseMat>(nlp_index_info.nnz_jac_g + 2 * nlp_index_info.nCon, nCon_QP, nVar_QP);
        H_ = make_shared<qpOASESSparseMat>(nlp_index_info.nnz_h_lag, nVar_QP, nVar_QP);
        qp_ = std::make_shared<qpOASES::SQProblem>((qpOASES::int_t) nVar_QP, (qpOASES::int_t) nCon_QP);
        return true;
    }

    /**
     * @name Constructor which also initializes the qpOASES SQProblem objects
     * @param nlp_index_info the number of variables in QP problem
     * @param nCon_QP the number of constraints in QP problem (the number of rows of A)
     */
    qpOASESInterface::qpOASESInterface(Index_info nlp_index_info) {
        //FIXME: the qpOASES does not accept any extra input
        allocate(nlp_index_info);


    }

    /**Default destructor*/
    qpOASESInterface::~qpOASESInterface() = default;

/**
 * @name This function solves the QP problem specified in the data, with given options. After the QP
 * is solved, it updates the stats, adding the iteration number used to solve the QP to the qp_iter
 * in object stats
 */
    bool qpOASESInterface::optimizeQP(shared_ptr<Stats> stats, shared_ptr<Options> options) {
//        qpOASESMatrixAdapter(H, H_qpOASES_);
//        qpOASESMatrixAdapter(A, A_qpOASES_);
//        H_qpOASES_->createDiagInfo();
//
//        qpOASES::int_t nWSR = options->qp_maxiter;//TODO modify it
//        if (stats->qp_iter == 0) {//if haven't solve any QP before then initialize the first qp
//            qpOASES::Options qp_options;
//            if (options->qpPrintLevel == 0)//else use the default print level in qpOASES
//                qp_options.printLevel = qpOASES::PL_NONE;
//            qp_->setOptions(qp_options);
//            qp_->init(H_qpOASES_.get(), g->vector(), A_qpOASES_.get(), lb->vector(), ub->vector(), lbA->vector(), ubA->vector(), nWSR,
//                      0);
//        } else
//            qp_->hotstart(H_qpOASES_.get(), g->vector(), A_qpOASES_.get(), lb->vector(), ub->vector(), lbA->vector(), ubA->vector(),
//                          nWSR, 0);
//        stats->qp_iter_addValue((int) nWSR);
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

//
///**
// * This function transforms the representation form of the sparse matrix in triplet form
// * to adapt the format required by the qpOASES(Harwell-Boeing Sparse Matrix).
// *
// * @param M_in_triplet Matrix objects contains data in triplet form
// * @param M_result     Matrix object prepared as the input for qpOASES
// */
//    bool qpOASESInterface::qpOASESMatrixAdapter(shared_ptr<Matrix> M_in_triplet,
//                                                shared_ptr<qpOASES::SparseMatrix> M_result) {
//
//        if (!A_.isinitialized) {
//            std::vector<std::tuple<int, int, Number, int>> tmp;
//            for (int i = 0; i < M_in_triplet->EntryNum(); i++) {
//                tmp.push_back(
//                        make_tuple(M_in_triplet->RowIndex_at(i), M_in_triplet->ColIndex_at(i),
//                                   M_in_triplet->MatVal_at(i),
//                                   M_in_triplet->order_at(i)));
//            }
//            for (int i = 0; i < M_in_triplet->num_I(); i++)
//                for (int j = 0; j < M_in_triplet->size_I_at(i); j++)
//                    tmp.push_back(make_tuple(M_in_triplet->RowIndex_I_at(i) + j, M_in_triplet->ColIndex_I_at(i) + j,
//                                             M_in_triplet->sign_I_at(i), tmp.size()));
//            std::sort(tmp.begin(), tmp.end(), tuple_sort_rule);
//            //copy the order information back
//            for (int i = 0; i < M_in_triplet->EntryNum(); i++) {
//                M_in_triplet->get_order(i, get<3>(tmp[i]));
//            }
//            delete_repetitive_entry(tmp);
//            for(int i  = 0; i<tmp.size(); i++){
//                cout<<get<0>(tmp[i])<<"    ";
//                cout<<get<1>(tmp[i])<<"    ";
//                cout<<get<2>(tmp[i])<<"    ";
//                cout<<get<3>(tmp[i])<<endl;
//            }
//            initialize_qpOASES_input(tmp, M_result, M_in_triplet->RowNum(), M_in_triplet->ColNum(), true);
//            A_.isinitialized = true;
//            tmp.clear();
//        } else {
//            update_qpOASES_input(M_in_triplet, M_result);
//        }
//        return true;
//    }
//
//    bool qpOASESInterface::qpOASESMatrixAdapter(shared_ptr<Matrix> M_in_triplet,
//                                                shared_ptr<qpOASES::SymSparseMat> M_result) {
//        if (!H_.isinitialized) {
//            std::vector<std::tuple<int, int, Number, int>> tmp;
//            for (int i = 0; i < M_in_triplet->EntryNum(); i++) {
//                tmp.push_back(
//                        make_tuple(M_in_triplet->RowIndex_at(i), M_in_triplet->ColIndex_at(i),
//                                   M_in_triplet->MatVal_at(i),
//                                   M_in_triplet->order_at(i)));
//            }
//            std::sort(tmp.begin(), tmp.end(), tuple_sort_rule);
//            //copy the order information back
//            for (int i = 0; i < M_in_triplet->EntryNum(); i++) {
//                M_in_triplet->get_order(i, get<3>(tmp[i]));
//            }
//
//            delete_repetitive_entry(tmp);
//            for(int i  = 0; i<tmp.size(); i++){
//                cout<<get<0>(tmp[i])<<"    ";
//                cout<<get<1>(tmp[i])<<"    ";
//                cout<<get<2>(tmp[i])<<"    ";
//                cout<<get<3>(tmp[i])<<endl;
//            }
//            initialize_qpOASES_input(tmp, M_result, M_in_triplet->RowNum(), M_in_triplet->ColNum(), false);
//            H_.isinitialized = true;
//            tmp.clear();
//        } else {
//            update_qpOASES_input(M_in_triplet, M_result);
//        }
//
//        return true;
//    }
//
///**
// * This is part of qpOASESMatrixAdapter
// *@name process the input data. Transform sorted data in tuple form to the format required by qpOASES QProblem class.
// * It will initialize the matrix object required by the qpOASES
// *
// * @tparam T either be SparseMatrix or SymSparseMat
// * @param input The input data in tuple form, it is the form used to store data in Matrix class
// * @param results The Matrix to be initialize
// * @param RowNum The number of rows of the matrix
// * @param ColNum The number of columns of the matrix
// */
//    template<typename T>
//    bool qpOASESInterface::initialize_qpOASES_input(vector<tuple<int, int, Number, int>> input,
//                                                    shared_ptr<T> results, Index RowNum, Index ColNum,
//                                                    bool isA) {
//        if (isA) {
//            A_.RowInd_ = new qpOASES::sparse_int_t[(int) input.size()];
//            A_.ColInd_ = new qpOASES::sparse_int_t[ColNum + 1];
//            A_.MatVal_ = new qpOASES::real_t[(int) input.size()];
//
//            int j = get<1>(input[0]);
//            for (int i = 0; i < (int) input.size(); i++) {
//                A_.RowInd_[i] = (qpOASES::sparse_int_t) get<0>(input[i]) - 1;
//                A_.ColInd_[i] = (qpOASES::real_t) get<2>(input[i]);
//                if (get<1>(input[i]) == j)A_.ColInd_[j]++;
//                else {
//                    A_.ColInd_[get<1>(input[i])] = (qpOASES::sparse_int_t) A_.ColInd_[get<1>(input[i - 1])] + 1;
//                    j = get<1>(input[i]);
//                }
//
//            }
//            for (int i = get<1>(input[input.size() - 1]); i <= ColNum; i++) {
//                A_.ColInd_[i] = A_.ColInd_[get<1>(input[input.size() - 1])];
//            }
//            for(int i =0; i<input.size();i++) cout<<A_.RowInd_[i]<<"   ";
//            cout<<" "<<endl;
//            for(int i =0; i<=ColNum;i++) cout<<A_.ColInd_[i]<<"   ";
//            cout<<" "<<endl;
//            for(int i =0; i<input.size();i++) cout<<A_.MatVal_at(i)<<"   ";
//            cout<<" "<<endl;
//            results = std::make_shared<T>(RowNum, ColNum, A_.RowInd_, A_.ColInd_, A_.MatVal_);
//        } else {
//            H_.RowInd_ = new qpOASES::sparse_int_t[(int) input.size()]();
//            H_.ColInd_ = new qpOASES::sparse_int_t[ColNum + 1]();
//            H_.MatVal_ = new qpOASES::real_t[(int) input.size()]();
//
//            int j = get<1>(input[0]);
//            for (int i = 0; i < (int) input.size(); i++) {
//                H_.RowInd_[i] = (qpOASES::sparse_int_t) get<0>(input[i]) - 1;
//                H_.MatVal_[i] = (qpOASES::real_t) get<2>(input[i]);
//                if (get<1>(input[i]) == j)H_.ColInd_[j]++;
//                else {
//                    H_.ColInd_[get<1>(input[i])] = (qpOASES::sparse_int_t) H_.ColInd_[get<1>(input[i - 1])] + 1;
//                    j = get<1>(input[i]);
//                }
//
////                for(int k =0; k<=ColNum;k++) cout<<H_tmp_.ColInd_[k];
////                cout<<"  "<<endl;
//            }
//            for (int i = get<1>(input[input.size() - 1]); i <= ColNum; i++) {
//                H_.ColInd_[i] = H_.ColInd_[get<1>(input[input.size() - 1])];
//            }
////            for(int i =0; i<input.size();i++) cout<<H_tmp_.RowInd_[i]<<endl;
////            for(int i =0; i<=ColNum;i++) cout<<H_tmp_.ColInd_[i]<<endl;
////            for(int i =0; i<input.size();i++) cout<<H_tmp_.MatVal_at(i)<<endl;
//            results = std::make_shared<T>(RowNum, ColNum, H_.RowInd_, H_.ColInd_, H_.MatVal_);
//
//        }
//        return true;
//    }
//
//
}//SQPHOTSTART

