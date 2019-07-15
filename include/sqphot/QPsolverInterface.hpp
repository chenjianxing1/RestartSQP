#ifndef SQPHOTSTART_QPSOLVER_INTERFACE_HPP
#define SQPHOTSTART_QPSOLVER_INTERFACE_HPP

#include <memory>
#include <vector>
#include <sqphot/Vector.hpp>
#include <sqphot/Matrix.hpp>
#include <qpOASES.hpp>
#include <sqphot/Stats.hpp>
#include <sqphot/Options.hpp>
//#include <sqphot/SQPDebug.hpp>


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
        QPSolverInterface() {}
        
        /** Default destructor*/
        virtual ~QPSolverInterface() {}
        
        /**
         * @name Solve a regular QP with given data and options.
         *
         * overload this method to optimize a QP with the data specified, update the stats by adding the iteration
         * number used to solve this QP to stats.qp_iter
         */
        virtual bool optimizeQP(shared_ptr<Stats> stats, shared_ptr<Options> options) = 0;
        
        
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
        
        
        virtual  shared_ptr<Vector> &getLb() = 0;
        
        virtual  shared_ptr<Vector> &getUb() = 0;
        
        virtual  shared_ptr<Vector> &getLbA() = 0;
        
        virtual  shared_ptr<Vector> &getUbA()  = 0;
        
        virtual  shared_ptr<Vector> &getG() = 0;
        
        //        virtual  shared_ptr<SpMatrix> &getH() =0;
        //
        //        virtual  shared_ptr<SpMatrix> &getA() = 0;
        
        
    private:
        
        /** Copy Constructor */
        QPSolverInterface(const QPSolverInterface &);
        
        /** Overloaded Equals Operator */
        void operator=(const QPSolverInterface &);
    };
    
    
    /**
     *This is a derived class of QPsolverInterface. It uses qpOASES as the QP solver which favors the hotstart
     * option. It is used as the default QP solver for SQPhostart.
     */
    class qpOASESInterface : public QPSolverInterface {
    public:
        
        virtual ~qpOASESInterface();
        
        /**
         * @name Constructor which also initializes the qpOASES SQProblem objects
         * @param nlp_index_info the number of variables in QP problem
         * @param nCon_QP the number of constraints in QP problem (the number of rows of A)
         */
        qpOASESInterface(Index_info nlp_index_info, QPType qptype);    //number of constraints in the QP problem
        
        virtual bool optimizeQP(shared_ptr<Stats> stats, shared_ptr<Options> options);
        
        /**
         * @name copy the optimal solution of the QP to the input pointer
         *
         * @param x_optimal a pointer to an empty array with allocated memory equals to sizeof(double)*number_variables
         *
         */
        
        inline bool get_optimal_solution(double* p_k);
        
        /**
         * @name copy the multipliers of the QP to the input pointer
         *
         * @param y_k   a pointer to an array with allocated memory equals to sizeof(double)*(num_variable+num_constraint)
         */
        inline bool get_multipliers(double* y_k);
        
        /**
         *@name get the objective value from the QP solvers
         *
         * @return the objective function value of the QP problem
         */
        inline double get_obj_value();
        
        
        shared_ptr<Vector> &getLb() ;
        
        shared_ptr<Vector> &getUb() ;
        
        shared_ptr<qpOASESSparseMat> &getH();
        
        shared_ptr<qpOASESSparseMat> &getA();
        
        shared_ptr<Vector> &getLbA();
        
        shared_ptr<Vector> &getUbA();
        
        shared_ptr<Vector> &getG()  ;
        
        
    public:
        shared_ptr<qpOASES::SQProblem> qp_;// the qpOASES object used for solving a qp
        
        
    private:
        shared_ptr<qpOASES::SymSparseMat> H_qpOASES_;
        shared_ptr<qpOASES::SparseMatrix> A_qpOASES_;
        shared_ptr<Vector> lb_;  // lower bounds of x
        shared_ptr<Vector> ub_;  // upper bounds of x
        shared_ptr<Vector> lbA_; //lower bounds of Ax
        shared_ptr<Vector> ubA_; //upper bounds of Ax
        shared_ptr<Vector> g_;
        shared_ptr<qpOASESSparseMat> H_;
        shared_ptr<qpOASESSparseMat> A_;
        bool firstQPsolved = false;
        
    private:
        
        /** default constructor*/
        qpOASESInterface();
        
        bool allocate(Index_info nlp_index_info, QPType qptype);
        
        /** Copy Constructor */
        qpOASESInterface(const qpOASESInterface &);
        
        /** Overloaded Equals Operator */
        void operator=(const qpOASESInterface &);
        
    };
}//SQPHOTSTART
#endif


