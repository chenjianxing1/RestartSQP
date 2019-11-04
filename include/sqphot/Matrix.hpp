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




}

#endif //SQPHOTSTART_MATRIX_HPP_

