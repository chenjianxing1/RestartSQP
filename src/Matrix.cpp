#include <sqphot/Matrix.hpp>

namespace SQPhotstart {

    /** Constructor for an empty matrix with N non-zero entries*/
    SpMatrix::SpMatrix(int EntryNum, int RowNum, int ColNum) :
            RowIndex_(NULL),
            ColIndex_(NULL),
            MatVal_(NULL),
            order_(NULL),
            EntryNum_(EntryNum),
            RowNum_(RowNum),
            ColNum_(ColNum) {
        //do nothing unless any data is to be assigned
        RowIndex_ = new int[EntryNum];
        ColIndex_ = new int[EntryNum];
        MatVal_ = new Number[EntryNum];
        order_ = new int[EntryNum];
        //initialize the order to 0:N-1
        for (int i = 0; i < EntryNum; i++) {
            order_[i] = i;
        }

    }

    /** Default destructor */
    SpMatrix::~SpMatrix() {
        freeMemory();
    }


    //    /**
//     * @name allocate the data to the class members
//     *
//     * @param RowIndex the row index of a entry in a matrix, starting from 1
//     * @param ColIndex the column index of a entry in a matrix, starting from 1
//     * @param MatVal   the entry value corresponding to (RowIndex,ColIndex)
//     *
//     */
//    bool SparseMatrix::assignMatrix(SQPhotstart::Index *RowIndex, SQPhotstart::Index *ColIndex,
//                              SQPhotstart::Number *MatVal) {
//        RowIndex_ = RowIndex;
//        ColIndex_ = ColIndex;
//        MatVal_ = MatVal;
//
//        return true;
//    }

    /**
     *@name print the sparse matrix in triplet form
     */
    void SpMatrix::print() {
        std::cout << "Row Column Entry Order" << std::endl;
        for (int i = 0; i < EntryNum_; i++) {
            std::cout << RowIndex_[i] << "    ";
            std::cout << ColIndex_[i] << "    ";
            std::cout << MatVal_[i] << "    ";
            std::cout << order_[i] << std::endl;
        }
    }

    /** free all memory*/
    bool SpMatrix::freeMemory() {
        delete[] RowIndex_;
        RowIndex_ = NULL;
        delete[] ColIndex_;
        ColIndex_ = NULL;
        delete[] MatVal_;
        MatVal_ = NULL;
        delete[] order_;
        order_ = NULL;
        return true;
    }

    bool SpMatrix::get_MatVal_at(int location, int value_to_assign) {
        MatVal_[location] = value_to_assign;
        return true;
    }

    bool SpMatrix::get_order(int location, int order_to_assign) {
        order_[location] = order_to_assign;
        return true;
    }


    /** MATRIX WITH IDENTITY*/


//    /** Constructor which assigns the number of nonzero entry and number of identity matrix to
//     * class members and allocate memory to the class member to be used*/
//    MatrixWithIdentity::MatrixWithIdentity(int num_nontrivial_entry, int num_of_I)
//            :
//            num_I_(num_of_I),
//            num_I_assigned_(0),
//            RowIndex_I_(NULL),
//            ColIndex_I_(NULL),
//            size_I_(NULL),
//            sign_I_(NULL) {
//        EntryNum_ = num_nontrivial_entry;
//        RowIndex_I_ = new int[num_of_I]();
//        ColIndex_I_ = new int[num_of_I]();
//        size_I_ = new int[num_of_I]();
//        sign_I_ = new int[num_of_I]();
//    }
//
//
//    /** Default destructor, will free all memory*/
//    MatrixWithIdentity::~MatrixWithIdentity() {
//        freeMemory();
//    }
//
//    /**
//     * This method adds the data of a identity matrix to be stored in the class members
//     *
//     * @param row_index the starting row index of the identity matrix
//     * @param col_index the starting column index of the identity matrix
//     * @param size      the size of the identity matrix
//     * @param isPositive is the entry positive or negative? If isPositive==true, then the sign_I
//     * 		     will record 1; otherwise, it will record -1
//     */
//    bool MatrixWithIdentity::add_I(const SQPhotstart::Index row_index, const SQPhotstart::Index col_index,
//                                   const SQPhotstart::Index size, bool isPositive) {
//
//
//        RowIndex_I_[num_I_assigned_] = row_index;
//        ColIndex_I_[num_I_assigned_] = col_index;
//        size_I_[num_I_assigned_] = size;
//
//        if (isPositive) sign_I_[num_I_assigned_] = 1;
//        else sign_I_[num_I_assigned_] = -1;
//        num_I_assigned_++;
//        return true;
//    }
//
//
//    /** Free all the memory ヽ( ^∀^)ﾉ*/
//    bool MatrixWithIdentity::freeMemory() {
//        delete[] RowIndex_;
//        RowIndex_ = NULL;
//        delete[] ColIndex_;
//        ColIndex_ = NULL;
//        delete[] MatVal_;
//        MatVal_ = NULL;
//        delete[] order_;
//        order_ = NULL;
//        delete[] ColIndex_I_;
//        ColIndex_I_ = NULL;
//        delete[] RowIndex_I_;
//        RowIndex_I_ = NULL;
//        delete[] size_I_;
//        size_I_ = NULL;
//        delete[] sign_I_;
//        sign_I_ = NULL;
//        return true;
//    }

    void qpOASESSparseMat::print() {}

    qpOASESSparseMat::~qpOASESSparseMat() { freeMemory(); }

    qpOASESSparseMat::qpOASESSparseMat(const SpMatrix &) {

    }

    qpOASESSparseMat::qpOASESSparseMat(int RowNum, int ColNum) :
            RowIndex_(NULL),
            ColIndex_(NULL),
            MatVal_(NULL),
            order_(NULL),
            isinitialized(false),
            RowNum_(RowNum),
            ColNum_(ColNum) {
    }
}

