#ifndef SQPHOTSTART_MATRIX_HPP_
#define SQPHOTSTART_MATRIX_HPP_

#include <sqphot/Types.hpp>
#include <sqphot/Vector.hpp>
#include <algorithm>
#include <vector>
#include <memory>
#include <tuple>
#include <cassert>
#include <qpOASES.hpp>

namespace SQPhotstart {
    //TODO: have a virtual base class for SparseMatrix, DenseMatrix,....
    class Matrix {
        
    public:
        /** Default constructor*/
        Matrix() = default;
        
        /** Default destructor*/
        virtual ~Matrix() = default;
        
        virtual void print() = 0;
    };
    
    class qpOASESSparseMat;
    
    typedef struct{
        int irow1 = 0;
        int jcol1 = 0;
        int size = 0;
        double scaling_factor1 = 1.0;
        int irow2 = 0;
        int jcol2 = 0;
        double scaling_factor2 = -1.0;
    }Identity2Info;
    
    /**
     *
     * This is the class for SparseMatrix, it stores the sparse matrix data in triplet,
     * in class member @_vector. It contains the methods that can copy a Matrix,
     * allocate data to the class member, and perform a matrix vector multiplication.
     *
     */
    
    class SpMatrix : public Matrix {
    public:
        
        /** Default constructor */
        SpMatrix():
        RowIndex_(NULL),
        ColIndex_(NULL),
        MatVal_(NULL),
        order_(NULL)
        {};
        
        /** Constructor for an empty Sparse Matrix with N non-zero entries*/
        SpMatrix(int nnz, int RowNum, int ColNum);
        
        /** Default destructor */
        virtual ~SpMatrix();
        
        /**
         * @name allocate the data to the class members
         *
         * @param RowIndex the row index of a entry in a matrix, starting from 1
         * @param ColIndex the column index of a entry in a matrix, starting from 1
         * @param MatVal   the entry value corresponding to (RowIndex,ColIndex)
         *
         */
        
        bool assignMatrix(Index* RowIndex, Index* ColIndex, Number* MatVal);
        
        /**
         *@name print the sparse matrix in triplet form
         */
        void print();
        
        
        
        /**
         * @name Times a matrix with a vector p, the pointer to the matrix-vector product will be
         * stored in the class member of another Vector class object called "result"
         * */
        virtual bool
        times(double* p,
              std::shared_ptr<SpMatrix> result  // the class object that will store the pointer of the matrix-vector product
        ) {
            return false;
        }
        
        
        
        //        virtual bool getStructure(shared_ptr<const SpMatrix> rhs){
        //            getStructure_sub(rhs);
        //            return true;};
        //
        //        virtual bool addIdentityMat(int irow, int jcol, int size, double scaling_factor){
        //            addIdentityMat_sub(irow,jcol,size,scaling_factor);
        //        };
        //
        //        virtual bool setMatVal(const double* MatVal){
        //            getMatVal_sub(MatVal);
        //            return true;
        //        };
        
        double onenorm() { return 0; };
        
        double infnorm() { return 0; };
        
        /**Extract Matrix info*/
        inline int ColNum() { return ColNum_; }
        
        inline int RowNum() { return RowNum_; }
        
        inline int EntryNum() { return EntryNum_; }
        
        inline int EntryNum() const { return EntryNum_; }
        
        inline int* RowIndex() { return RowIndex_; }
        
        inline int* ColIndex() { return ColIndex_; }
        
        inline double* MatVal() { return MatVal_; }
        
        inline int* order() { return order_; }
        
        inline const int* RowIndex() const { return RowIndex_; }
        
        inline const int* ColIndex() const { return ColIndex_; }
        
        inline const double* MatVal() const { return MatVal_; }
        
        inline const int* order() const { return order_; }
        
        inline int RowIndex_at(int i) { return RowIndex_[i]; }
        
        inline int ColIndex_at(int i) { return ColIndex_[i]; }
        
        inline double MatVal_at(int i) { return MatVal_[i]; }
        
        inline int order_at(int i) { return order_[i]; }
        
        inline bool get_order(int location, int order_to_assign);
        
        inline bool get_MatVal_at(int location, int value_to_assign);
        
        
        /** Private Method */
    private:
        
        /** free all memory*/
        bool freeMemory();
        
        /** Copy Constructor */
        SpMatrix(const SpMatrix &);
        
        /** Overloaded Equals Operator */
        void operator=(const SpMatrix &);
        /** Private Class Members*/
        
        
    private:
        double* MatVal_;  // the entry data of a matrix
        int* order_;    // the corresponding original position of a matrix entry
        int* RowIndex_;// the row number of a matrix entry
        int* ColIndex_;// the column number of a matrix entry
        int ColNum_;    // the number columns of a matrix
        int RowNum_;    // the number of rows of a matrix
        int EntryNum_;  // number of non-zero entries in  matrix
    };
    
    /**
     *
     *
     */
    class qpOASESSparseMat : public SpMatrix {
        
    public:
        
        qpOASESSparseMat(int nnz, int RowNum, int ColNum);
        
        ~qpOASESSparseMat();
        
        virtual bool setMatVal(const double* MatVal);
        
        /**Method that convert automatically convert data in triplet form
         * (stored in spMatrix) to Harwell-Boeing Sparse Matrix
         *@param rhs
         *@param numI number of Identity matrix
         */
        virtual bool setStructure(std::shared_ptr<const SpMatrix> rhs, Identity2Info I_info);
        
        virtual bool addIdentityMat(int irow, int jcol, int size, double scaling_factor);
        
        bool updateMatVal(const qpOASES::real_t) { return false; }
        
        
        void print() override;
        
        inline const qpOASES::sparse_int_t* RowIndex() const { return RowIndex_; }
        
        inline const qpOASES::sparse_int_t* ColIndex() const { return ColIndex_; }
        
        inline const qpOASES::real_t* MatVal() const { return MatVal_; }
        
        inline const int* order() const { return order_; }
        
        /**Private methods*/
    private:
        
        /**Default constructor*/
        qpOASESSparseMat();
        
        
        /** free all memory*/
        bool freeMemory();
        
        /** Copy Constructor */
        qpOASESSparseMat(const qpOASESSparseMat &);
        
        /** Overloaded Equals Operator */
        void operator=(const qpOASESSparseMat &);
        
        /**
         * This is part of qpOASESMatrixAdapter
         * @name This is the sorted rule that used to sort data, first based on column index then based on row index
         */
        static bool
        tuple_sort_rule(const std::tuple<int, int, int> left, const std::tuple<int, int, int> right) {
            if (std::get<1>(left) < std::get<1>(right)) return true;
            else if (std::get<1>(left) > std::get<1>(right)) return false;
            else {
                if (std::get<0>(left) < std::get<0>(right))return true;
                else return false;
            }
        }
        
        /** Private members*/
    private:
        qpOASES::sparse_int_t* RowIndex_;
        qpOASES::sparse_int_t* ColIndex_;
        qpOASES::real_t* MatVal_;
        int* order_;
        int RowNum_;
        int ColNum_;
        
        int EntryNum_;
        bool isinitialized;
        
    };
    
    
    //do nothing unless any data is to be assigned
    
    //   /**
    //    * @MatrixWithIdentity this is a derived class of SparseMatrix. In addition to the data stored in the class
    //    * members stored in the class SparseMatrix, it contains additional pointer to arrays stores the data for
    //    * identity matrix.
    //    *
    //    */
    //   class MatrixWithIdentity : public SparseMatrix {
    
    //   public:
    //       /** Constructor which assigns the number of nonzero entry and number of identity matrix to
    //        * class members and allocate memory to the class member to be used*/
    //       MatrixWithIdentity(int num_nontrivial_entry, int num_of_I);
    
    //       /** Default destructor, will free all memory*/
    //       virtual ~MatrixWithIdentity();
    
    
    //       /**
    //        * This method adds the data of a identity matrix to be stored in the class members
    //        *
    //        * @param row_index the starting row index of the identity matrix
    //        * @param col_index the starting column index of the identity matrix
    //        * @param size      the size of the identity matrix
    //        * @param isPositive is the entry positive or negative? If isPositive==true, then the sign_I
    //        *                      will record 1; otherwise, it will record -1
    //        */
    //       virtual bool
    //       add_I(const SQPhotstart::Index row_index, const SQPhotstart::Index col_index, const SQPhotstart::Index size,
    //             bool isPositive);
    
    //       /**
    //        * Get class member info
    //        */
    //       inline int num_I() { return num_I_; }
    
    //       inline int current_I_num() { return num_I_assigned_; }
    
    //       inline int* RowIndex_I() { return RowIndex_I_; }
    
    //       inline int* ColIndex_I() { return ColIndex_I_; }
    
    //       inline int* sizeI() { return size_I_; }
    
    //       inline int* signI() { return sign_I_; }
    
    //       inline int RowIndex_I_at(int i) { return RowIndex_I_[i]; }
    
    //       inline int ColIndex_I_at(int i) { return ColIndex_I_[i]; }
    
    //       inline int size_I_at(int i) { return size_I_[i]; }
    
    //       inline int sign_I_at(int i) { return sign_I_[i]; }
    
    //   private:
    //       int num_I_;             //number of identity submatrix in total
    //       int* ColIndex_I_;// each entry will record the starting column of the identity matrix
    //       int* RowIndex_I_;//each entry will record the starting row of the identity matrix
    //       int* size_I_;     //each entry will record the size of a identity matrix
    //       int* sign_I_;    //each entry will record the sign of identity matrix, either -1 or 1
    //       int num_I_assigned_;    //number of identity submatrix already assigned to this class
    
    //   private:
    //       /** Free all the memory ヽ( ^∀^)ﾉ*/
    //       bool freeMemory();
    
    //       /** Default constructor*/
    //       MatrixWithIdentity();
    //   };
    
    
}
#endif //SQPHOTSTART_MATRIX_HPP_

