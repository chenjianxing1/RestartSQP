#include <sqphot/Matrix.hpp>



namespace SQPhotstart {
/** Constructor for an empty matrix with N non-zero
 * entries*/

SpTripletMat::SpTripletMat(int nnz, int RowNum, int ColNum, bool isSymmetric,
                           bool allocate) :
    RowIndex_(nullptr),
    ColIndex_(nullptr),
    MatVal_(nullptr),
    order_(nullptr),
    isAllocated_(allocate),
    isSymmetric_(isSymmetric) {
    EntryNum_ = nnz;
    RowNum_ = RowNum;
    ColNum_ = ColNum;
    //do nothing unless any data is to be assigned
    if(allocate) {
        RowIndex_ = new int[nnz]();
        ColIndex_ = new int[nnz]();
        MatVal_ = new double[nnz]();
        order_ = new int[nnz]();
        //initialize the order to 0:N-1
        for (int i = 0; i < nnz; i++) {
            order_[i] = i;
        }
    }

}


/** Default destructor */
SpTripletMat::~SpTripletMat() {
    if(isAllocated_)
        freeMemory();
}


/**
 *@name print the sparse matrix in triplet form
 */
//@{
void
SpTripletMat::print(const char* name, Ipopt::SmartPtr<Ipopt::Journalist> jnlst,
                    Ipopt::EJournalLevel level,
                    Ipopt::EJournalCategory category) const {

    std::cout << "Row    Column    Entry    Order" << std::endl;
    for (int i = 0; i < EntryNum_; i++) {
        printf("%d       ", RowIndex_[i]);
        printf("%d       ", ColIndex_[i]);
        printf("%10e     ", MatVal_[i]);
        printf("%d     \n", order_[i]);
    }
}

/**
 *@name print the sparse matrix in dense form
 */
//@{
void
SpTripletMat::print_full(const char* name, Ipopt::SmartPtr<Ipopt::Journalist> jnlst,
                         Ipopt::EJournalLevel level,
                         Ipopt::EJournalCategory category) const {
    char mat_val[99];
    auto dense_matrix = new double[RowNum_ * ColNum_]();

    for (int i = 0; i < EntryNum_; i++) {
        dense_matrix[ColNum_ * (RowIndex_[i] - 1) + ColIndex_[i] - 1] = MatVal_[i];
        if (isSymmetric_ && RowIndex_[i] != ColIndex_[i])
            dense_matrix[ColNum_ * (ColIndex_[i] - 1) + RowIndex_[i] -
                         1] = MatVal_[i];
    }
    if(!IsNull(jnlst)) {
//    if (name != nullptr) {
//            jnlst->Printf(Ipopt::J_DBG,Ipopt::J_MATRIX,name);
//            jnlst->Printf(Ipopt::J_DBG,Ipopt::J_MATRIX," =: {\n");
//    }
//    for (int i = 0; i < RowNum_; i++) {
//        for (int j = 0; j < ColNum_; j++) {
//            sprintf(mat_val, "%f  ", dense_matrix[i * ColNum() + j]);
//               jnlst->Print(Ipopt::J_DBG,Ipopt::J_MATRIX,mat_val);
//        }
//           jnlst->Printf(Ipopt::J_DBG,Ipopt::J_MATRIX,"\n");
//    }
//       jnlst->Printf(Ipopt::J_DBG,Ipopt::J_MATRIX,"}\n\n");
    } else {
        if(name!=nullptr)
            printf("%s =:{\n", name);

        for (int i = 0; i < RowNum_; i++) {
            for (int j = 0; j < ColNum_; j++) {
                printf("%10e  ", dense_matrix[i * ColNum() + j]);
            }
            printf("\n");
        }
        printf("}\n\n");
    }
    delete[] dense_matrix;
}
//@}


/** free all memory*/
void SpTripletMat::freeMemory() {
    if(isAllocated_) {
        delete[] RowIndex_;
        RowIndex_ = nullptr;
        delete[] ColIndex_;
        ColIndex_ = nullptr;
        delete[] MatVal_;
        MatVal_ = nullptr;
        delete[] order_;
        order_ = nullptr;
    }
}


void SpTripletMat::setMatValAt(int location, int value_to_assign) {

    MatVal_[location] = value_to_assign;
}


void SpTripletMat::setOrderAt(int location, int order_to_assign) {

    order_[location] = order_to_assign;
}


/**
 * @brief Times a matrix with a vector p, the pointer to the matrix-vector
 * product  will be stored in the class member of another Vector class object
 * called "result"
 */

void SpTripletMat::times(std::shared_ptr<const Vector> p,
                         std::shared_ptr<Vector> result) const {

    assert(ColNum_ == p->Dim());
    if (isSymmetric_) {
        result->set_zeros();
        for (int i = 0; i < EntryNum_; i++) {
            result->addNumberAt(RowIndex_[i] - 1, MatVal_[i] * p->values()
                                [ColIndex_[i] - 1]);
            if (RowIndex_[i] != ColIndex_[i]) {
                result->addNumberAt(ColIndex_[i] - 1, MatVal_[i] * p->values()
                                    [RowIndex_[i] - 1]);
            }
        }
    } else {
        result->set_zeros(); //set all entries to be 0
        for (int i = 0; i < EntryNum_; i++) {
            result->addNumberAt(RowIndex_[i] - 1, MatVal_[i] * p->values()
                                [ColIndex_[i] - 1]);
        }
    }
}


double SpTripletMat::InfNorm() {
    //TODO: test it!
    std::shared_ptr<Vector> rowSums = std::make_shared<Vector>(RowNum_);
    for (int i = 0; i < EntryNum_; i++) {
        rowSums->addNumberAt(RowIndex()[i] - 1, abs(MatVal_[i]));
    }

    double InfNorm = rowSums->getInfNorm();//same as calculating the MAX of an array

    return InfNorm;
}


double SpTripletMat::OneNorm() {

    std::shared_ptr<Vector> colSums = std::make_shared<Vector>(ColNum_);
    for (int i = 0; i < EntryNum_; i++) {
        colSums->addNumberAt(ColIndex()[i] - 1, abs(MatVal_[i]));
    }

    double OneNorm = colSums->getInfNorm();//same as calculating the MAX of an array

    return OneNorm;
}


void SpTripletMat::copy(std::shared_ptr<const SpTripletMat> rhs, bool deep_copy) {

    if(!deep_copy) {
        RowIndex_= (int*)rhs->RowIndex();
        ColIndex_= (int*)rhs->ColIndex();
        MatVal_ = (double*)rhs->MatVal();
        order_= (int*)rhs->order();
    }
    for (int i = 0; i < EntryNum_; i++) {
        RowIndex_[i] = rhs->RowIndex()[i];
        ColIndex_[i] = rhs->ColIndex()[i];
        MatVal_[i] = rhs->MatVal()[i];
        order_[i] = rhs->order()[i];
    }
}


bool SpTripletMat::isSymmetric() const {

    return isSymmetric_;
}


void SpTripletMat::transposed_times(std::shared_ptr<const Vector> p,
                                    std::shared_ptr<Vector> result) const {

    if (isSymmetric_) {
        times(p, result);
    } else {
        result->set_zeros(); //set all entries to be 0
        for (int i = 0; i < EntryNum_; i++) {
            result->addNumberAt(ColIndex_[i] - 1, MatVal_[i] * p->values()
                                [RowIndex_[i] - 1]);
        }
    }

}


/**
 * qpOASESSparseMatrix
 */
}




