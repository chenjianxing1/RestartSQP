/* Copyright (C) 2019
 * All Rights Reserved.
 *
 * Authors: Xinyi Luo
 * Date:2019-10
 */

#ifndef SQPHOTSTART_SPTRIPLETMAT_HPP_
#define SQPHOTSTART_SPTRIPLETMAT_HPP_

#include "sqphot/Matrix.hpp"
#include "sqphot/Utils.hpp"
#include "sqphot/Vector.hpp"
#include <memory>

namespace SQPhotstart {

/** This is the class for SparseMatrix, it stores the sparse matrix
 *  data in triplet, in class member @values_. It contains the methods
 *  that can copy a Matrix, allocate data to the class member, and
 *  perform a matrix vector multiplication.
 */

class SpTripletMat : public Matrix
{
public:
  /** constructor/destructor */
  //@{
  /** Constructor for an empty Sparse Matrix with N non-zero entries.
   *  It sets the property of the matrix, but */
  SpTripletMat(int nnz, int RowNum, int ColNum, bool isSymmetric = false,
               bool allocate = true);

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
  void print(const char* name,
             Ipopt::SmartPtr<Ipopt::Journalist> jnlst = nullptr,
             Ipopt::EJournalLevel level = Ipopt::J_ALL,
             Ipopt::EJournalCategory category = Ipopt::J_DBG) const override;

  /**
   * @brief print the sparse matrix in the sense form
   */
  void
  print_full(const char* name,
             Ipopt::SmartPtr<Ipopt::Journalist> jnlst = nullptr,
             Ipopt::EJournalLevel level = Ipopt::J_ALL,
             Ipopt::EJournalCategory category = Ipopt::J_DBG) const override;

  /**
   * @brief Times a matrix with a vector p, the pointer to the matrix-vector
   * product  will be stored in the class member of another Vector class object
   * called "result"
   * */
  virtual void times(std::shared_ptr<const Vector> p,
                     std::shared_ptr<Vector> result) const;

  /**
   * @brief Times the matrix transpose with a vector p, the pointer to the
   * matrix-vector
   * product  will be stored in the class member of another Vector class object
   * called "result"
   * */
  virtual void transposed_times(std::shared_ptr<const Vector> p,
                                std::shared_ptr<Vector> result) const;

  //@}
  /**
   * @brief get the dense matrix corresponding to the matrix data stored in
   * class
   * members. The stored matrix will be column oriented.
   */
  void get_dense_matrix(double* dense_matrix, bool row_oriented = true) const;

  /**
   * @brief calculate the one norm of the matrix
   *
   * @return the calculated one-norm
   */
  double OneNorm();

  bool isCompressedRow() override
  {
    return false;
  };

  bool isSymmetric() const override
  {
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
  virtual void copy(std::shared_ptr<const SpTripletMat> rhs,
                    bool deep_copy = true);

  /**@name Extract Matrix info*/
  //@{
  inline int ColNum() const
  {

    return ColNum_;
  }

  inline int RowNum() const
  {

    return RowNum_;
  }

  inline int EntryNum()
  {

    return EntryNum_;
  }

  inline int EntryNum() const override
  {

    return EntryNum_;
  }

  inline int* RowIndex() override
  {

    return RowIndex_;
  }

  inline int* ColIndex() override
  {

    return ColIndex_;
  }

  inline double* MatVal() override
  {

    return MatVal_;
  }

  inline int* order() override
  {

    return order_;
  }

  inline int RowIndex(int i) override
  {

    return RowIndex_[i];
  }

  inline int ColIndex(int i) override
  {

    return ColIndex_[i];
  }

  inline double MatVal(int i) override
  {

    return MatVal_[i];
  }

  inline int order(int i) override
  {

    return order_[i];
  }

  inline const int RowIndex(int i) const
  {

    return RowIndex_[i];
  }

  inline const int ColIndex(int i) const
  {

    return ColIndex_[i];
  }

  inline const double MatVal(int i) const
  {

    return MatVal_[i];
  }

  inline const int order(int i) const
  {

    return order_[i];
  }

  inline const int* RowIndex() const
  {

    return RowIndex_;
  }

  inline const int* ColIndex() const
  {

    return ColIndex_;
  }

  inline const double* MatVal() const
  {

    return MatVal_;
  }

  inline const int* order() const
  {

    return order_;
  }

  //@}

  inline void setMatValAt(int location, double value_to_assign)
  {

    MatVal_[location] = value_to_assign;
  }

  inline void setRowIndex(int i, int value)
  {
    RowIndex_[i] = value;
  }

  inline void setColIndex(int i, int value)
  {
    ColIndex_[i] = value;
  }

  void set_zero()
  {
    for (int i = 0; i < EntryNum_; i++) {
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
  SpTripletMat(const SpTripletMat&);

  /** Overloaded Equals Operator */
  void operator=(const SpTripletMat&);

  ///////////////////////////////////////////////////////////
  //                     PRIVATE  MEMBERS                  //
  ///////////////////////////////////////////////////////////

private:
  bool isAllocated_;
  bool isSymmetric_; /**< is the matrix symmetric, if yes, the non-diagonal
                               *data will only be stored for once*/
  double* MatVal_;   /**< the entry data of a matrix */
  int ColNum_;       /**< the number columns of a matrix */
  int EntryNum_;     /**< number of non-zero entries in  matrix */
  int RowNum_;       /**< the number of rows of a matrix */
  int* ColIndex_;    /**< the column number of a matrix entry */
  int* RowIndex_;    /**< the row number of a matrix entry */
  int* order_; /**< the corresponding original position of a matrix entry */
};
}

#endif
