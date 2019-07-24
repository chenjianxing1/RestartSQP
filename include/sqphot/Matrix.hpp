/* Copyright (C) 2019
 * All Rights Reserved.
 *
 * Authors: Xinyi Luo
 * Date:2019-07
 */
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

    class Matrix {

    public:
        /** Default constructor*/
        Matrix() {};

        /** Default destructor*/
        virtual ~Matrix() {};

        virtual void print() = 0;

    };

    class qpOASESSparseMat;

    typedef struct {
        int irow1 = 0;
        int jcol1 = 0;
        int size = 0;
        double scaling_factor1 = 1.0;
        int irow2 = 0;
        int jcol2 = 0;
        double scaling_factor2 = -1.0;
    } Identity2Info;

    /**
     *
     * This is the class for SparseMatrix, it stores the sparse matrix data in triplet,
     * in class member @_vector. It contains the methods that can copy a Matrix,
     * allocate data to the class member, and perform a matrix vector multiplication.
     *
     */

    class SpTripletMat : public Matrix {
    public:

        /** constructor/destructor */
        //@{

        /** Constructor for an empty Sparse Matrix with N non-zero entries*/
        SpTripletMat(int nnz, int RowNum, int ColNum, bool isSymmetric);

        SpTripletMat(int nnz, int RowNum, int ColNum);


        /** Default destructor */
        ~SpTripletMat() override ;

        //@}
        /**
         * @brief allocate the data to the class members
         *
         * @param RowIndex the row index of a entry in a matrix, starting from 1
         * @param ColIndex the column index of a entry in a matrix, starting from 1
         * @param MatVal   the entry value corresponding to (RowIndex,ColIndex)
         *
         */

        bool setMatrix(Index* RowIndex, Index* ColIndex, Number* MatVal){return false;}

        /**
         *@brief print the sparse matrix in triplet form
         */
        void print() const;

        void print() override;

        /**
         * @brief Times a matrix with a vector p, the pointer to the matrix-vector
         * product  will be stored in the class member of another Vector class object
         * called "result"
         * */
        virtual bool times(std::shared_ptr<const Vector> p,
                           std::shared_ptr<Vector> result);


        /**
         * @brief calculate the one norm of the matrix
         *
         * @return the calculated one-norm
         */
        double onenorm();

        /**
         * @brief calculate the infinity norm of the matrix
         *
         * @return the calculated inf-norm
         */
        double infnorm();;

        /**
         *@brief make a deep copy of a matrix information
         */
        virtual bool copy(std::shared_ptr<const SpTripletMat> rhs);;


        /**Extract Matrix info*/
        inline int ColNum() const { return ColNum_; }

        inline int RowNum() const { return RowNum_; }

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

        bool isSymmetric() const;

        inline bool setOrderAt(int location, int order_to_assign);

        inline bool setMatValAt(int location, int value_to_assign);


        /** Private Method */
    private:

        /** Default constructor */
        SpTripletMat();


    private:

        /** free all memory*/
        bool freeMemory();

        /** Copy Constructor */
        SpTripletMat(const SpTripletMat&);

        /** Overloaded Equals Operator */
        void operator=(const SpTripletMat&);
        /** Private Class Members*/


    private:
        double* MatVal_;  /**< the entry data of a matrix */
        int* order_;    /**< the corresponding original position of a matrix entry */
        int* RowIndex_;/**< the row number of a matrix entry */
        int* ColIndex_;/**< the column number of a matrix entry */
        int ColNum_;    /**< the number columns of a matrix */
        int RowNum_;    /**< the number of rows of a matrix */
        int EntryNum_;  /**< number of non-zero entries in  matrix */
        bool isSymmetric_;

    };

    /**
     *
     *
     */
    class qpOASESSparseMat : public Matrix {

    public:
        /** constructor/destructor */
        /**Default constructor*/
        qpOASESSparseMat(int RowNum, int ColNum, bool isSymmetric);

        //@{
        /**
         * @brief A constructor
         * @param nnz the number of nonzero entries
         * @param RowNum the number of rows
         * @param ColNum the number of columns
         */
        qpOASESSparseMat(int nnz, int RowNum, int ColNum);

        bool isSymmetric() const;

        /**
         * @brief Default destructor
         */
        ~qpOASESSparseMat() override;
        //@}


        /**
         * @brief set the Matrix values to the matrix, convert from triplet format to
         * Harwell-Boeing Matrix format.
         * @param MatVal entry values(orders are not yet under permutation)
         * @param I_info the 2 identity matrices information
         */
        virtual bool setMatVal(const double* MatVal, Identity2Info I_info);

        virtual bool setMatVal(std::shared_ptr<const SpTripletMat> rhs);

        /**
         * @brief setup the structure of the sparse matrix for solver qpOASES(should
         * be called only for once).
         *
         * This method will convert the strucutre information from the triplet form from a
         * SpMatrix object to the format required by the QPsolver qpOASES.
         *
         * @param rhs a SpMatrix object whose content will be copied to the class members
         * (in a different sparse matrix representations)
         * @param I_info the information of 2 identity sub matrices.
         *
         */
        bool setStructure(std::shared_ptr<const SpTripletMat> rhs, Identity2Info I_info);

        bool setStructure(std::shared_ptr<const SpTripletMat> rhs);


        bool updateMatVal(const double* MatVal) { return false; }

        /**
         * @brief print the matrix information
         */
        void print() override;


        /**
         * @brief make a deep copy of a matrix information
         */

        virtual bool copy(std::shared_ptr<const qpOASESSparseMat> rhs);

        /** Extract class member information*/
        //@{

        inline int EntryNum() const { return EntryNum_; }

        inline int ColNum() const { return ColNum_; }

        inline int RowNum() const { return RowNum_; }

        inline qpOASES::sparse_int_t* RowIndex() { return RowIndex_; }

        inline qpOASES::sparse_int_t* ColIndex() { return ColIndex_; }

        inline qpOASES::real_t* MatVal() { return MatVal_; }

        inline int* order() { return order_; }

        inline const qpOASES::sparse_int_t* RowIndex() const { return RowIndex_; }

        inline const qpOASES::sparse_int_t* ColIndex() const { return ColIndex_; }

        inline const qpOASES::real_t* MatVal() const { return MatVal_; }

        inline const int* order() const { return order_; }

        bool isIsinitialized() const;
        //@}

        /**Private methods*/
    private:



        /** free all memory*/
        bool freeMemory();

        /** Copy Constructor */
        qpOASESSparseMat(const qpOASESSparseMat&);

        /** Overloaded Equals Operator */
        void operator=(const qpOASESSparseMat&);

        /**
         * This is part of qpOASESMatrixAdapter
         * @brief This is the sorted rule that used to sort data, first based on column
         * index then based on row index
         */
        static bool
        tuple_sort_rule(const std::tuple<int, int, int> left,
                        const std::tuple<int, int, int> right) {
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
        bool isInitialised_;
        bool isSymmetric_;
    };


}
#endif //SQPHOTSTART_MATRIX_HPP_

