#ifndef SQPHOTSTART_QPSOLVER_INTERFACE_HPP
#define SQPHOTSTART_QPSOLVER_INTERFACE_HPP

#include <memory>
#include <vector>
#include <sqphot/Vector.hpp>
#include <sqphot/Matrix.hpp>
#include <qpOASES.hpp>
#include <sqphot/Stats.hpp>
#include <sqphot/Options.hpp>


using namespace std;
namespace SQPhotstart {

    /**
     * Base class for all standard QP solvers that use standard triplet matrix form and dense vectors.
     *
     * It can optimize QP problem in the following format
     *
     *  minimize 1/2 x^T H x + g^T x
     *  subject  lb_A<=Ax<=ub_A
     *              lb<=x<=ub
     */
    class QPSolverInterface {

    public:
        /** Default constructor*/
        QPSolverInterface() = default;

        /** Default destructor*/
        virtual ~QPSolverInterface() = default;

        /**
         * @name Solve a regular QP with given data and options.
         *
         * overload this method to optimize a QP with the data specified, update the stats by adding the iteration
         * number used to solve this QP to stats.qp_iter
         */
        virtual bool optimizeQP(shared_ptr<Matrix> H, shared_ptr<Vector> g, shared_ptr<Matrix> A,
                                shared_ptr<Vector> lbA, shared_ptr<Vector> ubA, shared_ptr<Vector> lb,
                                shared_ptr<Vector> ub, shared_ptr<Stats> stats, shared_ptr<Options> options) = 0;


        /**
         * @name copy the optimal solution of the QP to the input pointer
         *
         * @param x_optimal a pointer to an empty array with allocated memory euqal to sizeof(double)*number_variables
         */
        virtual bool get_optimal_solution(double* x_optimal) = 0;

        /**
         *@name get the objective value from the QP solvers
         *
         * @return the objective function value of the QP problem
         */
        virtual double get_obj_value() = 0;


        /**
         * @name copy the multipliers of the QP to the input pointer
         *
         * @param y_k   a pointer to an array with allocated memory
         */
        virtual bool get_multipliers(double* y_optimal) = 0;
    };

    /**
     *This is a derived class of QPsolverInterface. It uses qpOASES as the QP solver which favors the hotstart
     * option. It is used as the default QP solver for SQPhostart.
     *
     */
    class qpOASESInterface : public QPSolverInterface {
        typedef struct {
            qpOASES::sparse_int_t* RowInd_ = NULL;
            qpOASES::sparse_int_t* ColInd_ = NULL;
            qpOASES::real_t* MatVal_ = NULL;
            bool isinitialized = false;
        } qpOASESSparseMat;
    public:
        /**Default constructor*/
        qpOASESInterface();

        /**
         * @name Constructor which also initializes the qpOASES SQProblem objects
         * @param nVar_QP the number of variables in QP problem
         * @param nCon_QP the number of constraints in QP problem (the number of rows of A)
         */

        qpOASESInterface(const int nVar_QP,    //number of variables in the QP problem
                         const int nCon_QP);    //number of constraints in the QP problem

        /**Default destructor*/
        virtual ~qpOASESInterface();


        /**
         * @name This function solves the QP problem specified in the data, with given options. After the QP
         * is solved, it updates the stats, adding the iteration number used to solve the QP to the qp_iter
         * in object stats
         */
        virtual bool optimizeQP(shared_ptr<Matrix> H, shared_ptr<Vector> g, shared_ptr<Matrix> A,
                                shared_ptr<Vector> lbA, shared_ptr<Vector> ubA, shared_ptr<Vector> lb,
                                shared_ptr<Vector> ub, shared_ptr<Stats> stats, shared_ptr<Options> options);

        /**
         * @name copy the optimal solution of the QP to the input pointer
         *
         * @param x_optimal a pointer to an empty array with allocated memory equals to sizeof(double)*number_variables
         */
        inline bool get_optimal_solution(double* p_k) {
            _qp->getPrimalSolution(p_k);
            return true;
        }

        /**
         * @name copy the multipliers of the QP to the input pointer
         *
         * @param y_k   a pointer to an array with allocated memory equals to sizeof(double)*(num_variable+num_constraint)
         */
        inline bool get_multipliers(double* y_k) {
            _qp->getDualSolution(y_k);
            return true;
        }

        /**
         *@name get the objective value from the QP solvers
         *
         * @return the objective function value of the QP problem
         */
        inline double get_obj_value() {
            return (double) (_qp->getObjVal());
        }

        /**
         * @name This function transforms the representation form of the sparse matrix in triplet form
         * to adapt the format required by the qpOASES(Harwell-Boeing Sparse Matrix)
         *
         * @param M_in_triplet Matrix objects contains data in triplet form
         * @param M_result     Matrix object prepared as the input for qpOASES
         */
        virtual bool
        qpOASESMatrixAdapter(shared_ptr<Matrix> M_in_triplet, shared_ptr<qpOASES::SparseMatrix> M_result);

        virtual bool qpOASESMatrixAdapter(shared_ptr<Matrix> M_in_triplet, shared_ptr<qpOASES::SymSparseMat> M_result);

    public:
        shared_ptr<qpOASES::SQProblem> _qp;// the qpOASES object used for solving a qp

    private:
        shared_ptr<qpOASES::SymSparseMat> H_;
        shared_ptr<qpOASES::SparseMatrix> A_;
        qpOASESSparseMat A_tmp_;
        qpOASESSparseMat H_tmp_;

        /**
         * This is part of qpOASESMatrixAdapter
         *@name process the input data. Transform sorted data in tuple form to the format required by qpOASES QProblem class.
         * It will initialize the matrix object required by the qpOASES
         *
         * @tparam T either be SparseMatrix or SymSparseMat
         * @param input The input data in tuple form, it is the form used to store data in Matrix class
         * @param results The Matrix to be initialize
         * @param RowNum The number of rows of the matrix
         * @param ColNum The number of columns of the matrix
         */
        template<typename T>
        bool initialize_qpOASES_input(vector<tuple<int, int, Number, int>> input,
                                      shared_ptr<T> results, Index RowNum, Index ColNum, bool isA);

        template<typename T>
        bool update_qpOASES_input(shared_ptr<Matrix> input, shared_ptr<T> results) { return false; };


        /**
         *This is part of qpOASESMatrixAdapter
         * @name it will process the input Matrix data in tuple form by detecting if there is any repeated entry.
         * If there are any, it will add the MatVal with the same row index and column index together and keep only one
         * in the vector of tuples
         * @param mat the data in the Matrix object
         */
        bool delete_repetitive_entry(vector<tuple<int, int, Number, int>> mat) { return false; };

        /**
         * This is part of qpOASESMatrixAdapter
         * @name This is the sorted rule that used to sort data, first based on column index then based on row index
         */
        static bool
        tuple_sort_rule(const tuple<int, int, Number, int> left, const tuple<int, int, Number, int> right) {
            if (get<1>(left) < get<1>(right)) return true;
            if (get<1>(left) > get<1>(right)) return false;
            if (get<1>(left) == get<1>(right)) {
                if (get<0>(left) < get<0>(right))return true;
                else return false;
            }
        }

    };
}//SQPHOTSTART
#endif

