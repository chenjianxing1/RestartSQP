/* Copyright (C) 2019
 * All Rights Reserved.
 *
 * Authors: Xinyi Luo
 * Date:2019-07
 */
#ifndef SQPHOTSTART_MATRIX_HPP_
#define SQPHOTSTART_MATRIX_HPP_

#include <algorithm>
#include <cassert>
#include <cmath>
#include <memory>
#include <IpJournalist.hpp>
#include <qpOASES.hpp>
#include <IpTNLP.hpp>
#include <sqphot/Utils.hpp>
#include <sqphot/Types.hpp>
#include <sqphot/Vector.hpp>
#include <tuple>
#include <vector>


namespace SQPhotstart {

typedef struct {
    int irow1 = 0;
    int jcol1 = 0;
    int size = 0;
    int irow2 = 0;
    int jcol2 = 0;
} Identity2Info;


/**
 *@brief
 * This is a virtual base class...
 *
 */
class Matrix {

///////////////////////////////////////////////////////////
//                     PUBLIC  METHODS                   //
///////////////////////////////////////////////////////////

public:

    /** Default constructor*/
    Matrix() = default;


    /** Default destructor*/
    virtual ~Matrix() = default;

    /**
     * @brief Print Matrix. If matrix is sparse, then print it in sparse form.
     */
    virtual void
    print(const char* name = nullptr, Ipopt::SmartPtr<Ipopt::Journalist> jnlst = nullptr,
          Ipopt::EJournalLevel level = Ipopt::J_ALL, Ipopt::EJournalCategory category
          =Ipopt::J_DBG)
    const = 0;


    virtual void
    print_full(const char* name = nullptr, Ipopt::SmartPtr<Ipopt::Journalist> jnlst = nullptr,
               Ipopt::EJournalLevel level=Ipopt::J_ALL, Ipopt::EJournalCategory
               category=Ipopt:: J_DBG)
    const
        = 0;

    virtual bool isCompressedRow() = 0;

    virtual int EntryNum() const = 0;

    virtual bool isSymmetric() const =0;

    virtual double MatVal(int i) =0;

    virtual int ColIndex(int i) = 0;

    virtual int RowIndex(int i) = 0;

    virtual int order(int i) = 0;

    virtual double* MatVal() =0;

    virtual int* ColIndex() = 0;

    virtual int* RowIndex() = 0;

    virtual int* order() = 0;
///////////////////////////////////////////////////////////
//                     PRIVATE  METHODS                  //
///////////////////////////////////////////////////////////

private:
    /** Copy Constructor */
    Matrix(const Matrix &);


    /** Overloaded Equals Operator */
    void operator=(const Matrix &);

};



/**
 *
 * This is the class for SparseMatrix, it stores the sparse matrix data in triplet,
 * in class member @values_. It contains the methods that can copy a Matrix,
 * allocate data to the class member, and perform a matrix vector multiplication.
 *
 */

class SpTripletMat :
    public Matrix {
public:

    /** constructor/destructor */
    //@{
    /** Constructor for an empty Sparse Matrix with N non-zero entries*/
    SpTripletMat(int nnz, int RowNum, int ColNum, bool isSymmetric=false,
                 bool allocate=true);



    /**
     *@brief
     *
     */
    SpTripletMat(const double* data, int RowNum, int ColNum, bool row_oriented);


    /** Default destructor*/
    ~SpTripletMat() override;
    //@}

    /**
     *@brief print the sparse matrix in triplet form
     */
    void print(const char* name, Ipopt::SmartPtr<Ipopt::Journalist> jnlst = nullptr,
               Ipopt::EJournalLevel level =Ipopt::J_ALL,
               Ipopt::EJournalCategory category=Ipopt::J_DBG) const override;


    /**
     * @brief print the sparse matrix in the sense form
     */
    void print_full(const char* name, Ipopt::SmartPtr<Ipopt::Journalist> jnlst = nullptr,
                    Ipopt::EJournalLevel level =Ipopt::J_ALL,
                    Ipopt::EJournalCategory category=Ipopt::J_DBG) const override;


    /**
     * @brief Times a matrix with a vector p, the pointer to the matrix-vector
     * product  will be stored in the class member of another Vector class object
     * called "result"
     * */
    virtual void times(std::shared_ptr<const Vector> p,
                       std::shared_ptr<Vector> result) const;


    /**
     * @brief Times the matrix transpose with a vector p, the pointer to the matrix-vector
     * product  will be stored in the class member of another Vector class object
     * called "result"
     * */
    virtual void transposed_times(std::shared_ptr<const Vector> p,
                                  std::shared_ptr<Vector> result) const;


    //@}
    /**
     * @brief get the dense matrix corresponding to the matrix data stored in class
     * members. The stored matrix will be column oriented.
     */
    void get_dense_matrix(double* dense_matrix, bool row_oriented = true) const;


    /**
     * @brief calculate the one norm of the matrix
     *
     * @return the calculated one-norm
     */
    double OneNorm();


    bool isCompressedRow() override {
        return false;
    };


    bool isSymmetric() const override {
        return isSymmetric_;
    }


    /**
     * @brief calculate the infinity norm of the matrix
     *
     * @return the calculated inf-norm
     */
    double InfNorm();


    /**
     *@brief make a deep copy of a matrix information
     */
    virtual void copy(std::shared_ptr<const SpTripletMat> rhs, bool deep_copy = true);


    /**@name Extract Matrix info*/
    //@{
    inline int ColNum() const {

        return ColNum_;
    }


    inline int RowNum() const {

        return RowNum_;
    }


    inline int EntryNum() {

        return EntryNum_;
    }


    inline int EntryNum() const override {

        return EntryNum_;
    }


    inline int* RowIndex() override {

        return RowIndex_;
    }


    inline int* ColIndex() override {

        return ColIndex_;
    }


    inline double* MatVal() override {

        return MatVal_;
    }


    inline int* order() override {

        return order_;
    }

    inline int RowIndex(int i) override {

        return RowIndex_[i];
    }


    inline int ColIndex(int i ) override {

        return ColIndex_[i];
    }


    inline double MatVal(int i) override {

        return MatVal_[i];
    }


    inline int order(int i) override {

        return order_[i];
    }


    inline const int RowIndex(int i) const {

        return RowIndex_[i];
    }


    inline const int ColIndex(int i ) const {

        return ColIndex_[i];
    }


    inline const double MatVal(int i) const {

        return MatVal_[i];
    }


    inline const int order(int i) const {

        return order_[i];
    }

    inline const int* RowIndex() const {

        return RowIndex_;
    }


    inline const int* ColIndex() const {

        return ColIndex_;
    }


    inline const double* MatVal() const {

        return MatVal_;
    }


    inline const int* order() const {

        return order_;
    }


    //@}

    inline void setMatValAt(int location, double value_to_assign) {

        MatVal_[location] = value_to_assign;
    }

    inline void setRowIndex(int i, int value) {
        RowIndex_[i] = value;
    }

    inline void setColIndex(int i, int value) {
        ColIndex_[i] = value;
    }

    inline void setOrderAt(int location, int order_to_assign) {
        order_[location] = order_to_assign;
    }

    void set_zero() {
        for(int i = 0; i<EntryNum_; i++) {
            order_[i] = 0;
            MatVal_[i] = 0;
            ColIndex_[i] = 0;
            RowIndex_[i] = 0;
        }
    }



///////////////////////////////////////////////////////////
//                     PRIVATE  METHODS                  //
///////////////////////////////////////////////////////////
private:
    /** Default constructor*/

    SpTripletMat();


    /** free all memory*/
    void freeMemory();


    /** Copy Constructor */
    SpTripletMat(const SpTripletMat &);


    /** Overloaded Equals Operator */
    void operator=(const SpTripletMat &);


///////////////////////////////////////////////////////////
//                     PRIVATE  MEMBERS                  //
///////////////////////////////////////////////////////////

private:
    bool isAllocated_;
    bool isSymmetric_;/**< is the matrix symmetric, if yes, the non-diagonal
                                *data will only be stored for once*/
    double* MatVal_;  /**< the entry data of a matrix */
    int ColNum_;    /**< the number columns of a matrix */
    int EntryNum_;  /**< number of non-zero entries in  matrix */
    int RowNum_;    /**< the number of rows of a matrix */
    int* ColIndex_;/**< the column number of a matrix entry */
    int* RowIndex_;/**< the row number of a matrix entry */
    int* order_;    /**< the corresponding original position of a matrix entry */

};


/**
 *@brief This is a derived class of Matrix.
 * It strored matrix in Harwell-Boeing format which is required by qpOASES and QORE.
 * It contains method to transform matrix format from Triplet form to
 * Harwell-Boeing format and then stored to its class members
 */

class SpHbMat :
    public Matrix {

public:
    /** constructor/destructor */
    //@{

    /**Default constructor*/
    SpHbMat(int RowNum, int ColNum, bool isCompressedRow);


    /**
     * @brief A constructor
     * @param nnz the number of nonzero entries
     * @param RowNum the number of rows
     * @param ColNum the number of columns
     */
    SpHbMat(int nnz, int RowNum, int ColNum, bool isCompressedRow = false);


    /**
     * @brief constructor which generate matrix data directly from a dense matrix
     *
     */
    SpHbMat(const double* data, int RowNum, int ColNum, bool row_oriented = true,
            bool isCompressedRow = false);

    /**
     * @brief Default destructor
     */
    ~SpHbMat() override;
    //@}


    /**
     * @brief set the Matrix values to the matrix, convert from triplet format to
     * Harwell-Boeing Matrix format.
     * @param rhs entry values(orders are not yet under permutation)
     * @param I_info the 2 identity matrices information
     */
    virtual void setMatVal(std::shared_ptr<const SpTripletMat> rhs, Identity2Info
                           I_info);



    virtual void setMatVal(std::shared_ptr<const SpTripletMat> rhs);


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
//@{
    void setStructure(std::shared_ptr<const SpTripletMat> rhs, Identity2Info I_info);


    void setStructure(std::shared_ptr<const SpTripletMat> rhs);
//@}

    void get_dense_matrix(double* dense_matrix,bool row_oriented = true);

    /**
     * @brief print the matrix information
     */
    void print(const char* name = nullptr, Ipopt::SmartPtr<Ipopt::Journalist> jnlst =
                   nullptr,
               Ipopt::EJournalLevel level=Ipopt::J_ALL, Ipopt::EJournalCategory
               category=Ipopt::J_DBG)
    const
    override;


    void print_full(const char* name = nullptr, Ipopt::SmartPtr<Ipopt::Journalist> jnlst = nullptr,
                    Ipopt::EJournalLevel level=Ipopt::J_ALL, Ipopt::EJournalCategory
                    category=Ipopt::J_DBG) const
    override;




    void times(std::shared_ptr<const Vector> p,
               std::shared_ptr<Vector> result) const;

    void transposed_times(std::shared_ptr<const Vector> p,
                          std::shared_ptr<Vector> result) const;
    /**
     * @brief make a deep copy of a matrix information
     */

    virtual  void copy(std::shared_ptr<const SpHbMat> rhs);


    /**
     * @convert the matrix data stored in the class members to a triplet matrix
     * specified by rhs */
    shared_ptr<SpTripletMat> convert_to_triplet() const ;


    /** Extract class member information*/
    //@{

    inline int EntryNum() const  override {

        return EntryNum_;
    }


    inline int ColNum() const {

        return ColNum_;
    }


    inline int RowNum() const {

        return RowNum_;
    }

    inline int RowIndex(int i) override {

        return RowIndex_[i];
    }

    inline int ColIndex(int i ) override {

        return ColIndex_[i];
    }


    inline double MatVal(int i) override {

        return MatVal_[i];
    }


    inline int order(int i) override {

        return order_[i];
    }

    inline int* RowIndex() override {

        return RowIndex_;
    }


    inline int* ColIndex() override {

        return ColIndex_;
    }


    inline double* MatVal() override {

        return MatVal_;
    }


    inline int* order() override {

        return order_;
    }


    inline const int* RowIndex() const {

        return RowIndex_;
    }


    inline const int* ColIndex() const {

        return ColIndex_;
    }


    inline const double* MatVal() const {

        return MatVal_;
    }


    inline const int* order() const {

        return order_;
    }


    inline bool isSymmetric() const override {

        return isSymmetric_;
    };


    inline bool isinitialized() {

        return isInitialised_;
    };

    inline bool isCompressedRow() override {
        return isCompressedRow_;
    }

    /**
     *@brief write data to a file
     * Only works when DEBUG is enabled.
     */
    void write_to_file(const char* name,
                       Ipopt::SmartPtr<Ipopt::Journalist> jnlst,
                       Ipopt::EJournalLevel level,
                       Ipopt::EJournalCategory category,
                       Solver solver);


///////////////////////////////////////////////////////////
//                     PRIVATE  METHODS                  //
///////////////////////////////////////////////////////////

private:

    /**Default constructor*/
    SpHbMat();


    /** free all memory*/
    void freeMemory();


    /** Copy Constructor */
    SpHbMat(const SpHbMat &);


    /** Overloaded Equals Operator */
    void operator=(const SpHbMat &);


    void set_zero();


    template <typename T>
    static void print_tuple(vector<tuple<int,int,T>> tuple) {
        for(int i = 0; i<tuple.size(); i++) {
            printf("%d %d ", std::get<0>(tuple[i]),std::get<1>(tuple[i]));
            std::cout<<std::get<2>(tuple[i])<<std::endl;
        }
    }

    /**
     * @brief This is the sorted rule that used to sort data, first based on column
     * index then based on row index
     */
    static bool
    tuple_sort_rule_compressed_column(const std::tuple<int, int, int> left,
                                      const std::tuple<int, int, int> right) {


        if (std::get<1>(left) < std::get<1>(right)) return true;
        else if (std::get<1>(left) > std::get<1>(right)) return false;
        else {
            return std::get<0>(left) < std::get<0>(right);
        }
    }

    /**
     * @brief This is the sorted rule that used to sort data, first based on row
     * index then based on column index
     */
    static bool
    tuple_sort_rule_compressed_row(const std::tuple<int, int, int> left,
                                   const std::tuple<int, int, int> right) {


        if (std::get<0>(left) < std::get<0>(right)) return true;
        else if (std::get<0>(left) > std::get<0>(right)) return false;
        else {
            return std::get<1>(left) < std::get<1>(right);
        }
    }

///////////////////////////////////////////////////////////
//                     PRIVATE  MEMBERS                  //
///////////////////////////////////////////////////////////
private:
    bool isInitialised_;
    bool isSymmetric_;
    bool isCompressedRow_;
    int ColNum_;
    int EntryNum_;
    int RowNum_;
    int* order_;
    double * MatVal_;
    int * ColIndex_;
    int *RowIndex_;
};
}

#endif //SQPHOTSTART_MATRIX_HPP_

