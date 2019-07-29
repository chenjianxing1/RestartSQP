#include <sqphot/Matrix.hpp>

namespace SQPhotstart {

/** Constructor for an empty matrix with N non-zero
 * entries*/

    SpTripletMat::SpTripletMat(int nnz, int RowNum, int ColNum, bool isSymmetric) :
            RowIndex_(NULL),
            ColIndex_(NULL),
            MatVal_(NULL),
            order_(NULL),
            isSymmetric_(isSymmetric) {
        EntryNum_ = nnz;
        RowNum_ = RowNum;
        ColNum_ = ColNum;
        //do nothing unless any data is to be assigned
        RowIndex_ = new int[nnz]();
        ColIndex_ = new int[nnz]();
        MatVal_ = new Number[nnz]();
        order_ = new int[nnz]();
        //initialize the order to 0:N-1
        for (int i = 0; i < nnz; i++) {
            order_[i] = i;
        }

    }

/** Default destructor */
    SpTripletMat::~SpTripletMat() {
        freeMemory();
    }


/**
 *@name print the sparse matrix in triplet form
 */
//@{
    void SpTripletMat::print() const {
        std::cout << "Row Column Entry Order" << std::endl;
        for (int i = 0; i < EntryNum_; i++) {
            std::cout << RowIndex_[i] << "    ";
            std::cout << ColIndex_[i] << "    ";
            std::cout << MatVal_[i] << "    ";
            std::cout << order_[i] << std::endl;
        }
    }

/**
 *@name print the sparse matrix in dense form
 */
//@{
    void SpTripletMat::print_full(const char* name) const {
        if (name != NULL) {
            std::cout << name << "is" << std::endl;
        }
        auto dense_matrix = new double[RowNum_ * ColNum_]();

        for(int i = 0; i < EntryNum_; i++){
		dense_matrix[ColNum_*(RowIndex_[i]-1)+ColIndex_[i]-1] = MatVal_[i];
		if(isSymmetric_&& RowIndex_[i]!=ColIndex_[i])
			dense_matrix[ColNum_*(ColIndex_[i]-1)+RowIndex_[i]-1] = MatVal_[i];
	}

        for (int i = 0; i < RowNum_; i++) {
            for (int j = 0; j < ColNum_; j++) {
                std::cout << dense_matrix[i * ColNum_ + j] << " ";
            }
            std::cout << std::endl;
        }
        std::cout << std::endl;
        delete[] dense_matrix;
    }
//@}


/** free all memory*/
    bool SpTripletMat::freeMemory() {
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

    bool SpTripletMat::setMatValAt(int location, int value_to_assign) {
        MatVal_[location] = value_to_assign;
        return true;
    }


    bool SpTripletMat::setOrderAt(int location, int order_to_assign) {
        order_[location] = order_to_assign;
        return true;
    }


/**
 * @brief Times a matrix with a vector p, the pointer to the matrix-vector
 * product  will be stored in the class member of another Vector class object
 * called "result"
 */

    bool SpTripletMat::times(std::shared_ptr<const Vector> p,
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
        return true;
    }


    double SpTripletMat::infnorm() {
        double InfNorm = 0;
        //FIXME: finish it!

        return InfNorm;
    }


    double SpTripletMat::onenorm() {
        double OneNorm = 0;
        return OneNorm;
    }

    bool SpTripletMat::copy(std::shared_ptr<const SpTripletMat> rhs) {
        RowNum_ = rhs->RowNum();
        ColNum_ = rhs->ColNum();
        EntryNum_ = rhs->EntryNum();

        for (int i = 0; i < EntryNum_; i++) {
            RowIndex_[i] = rhs->RowIndex()[i];
            ColIndex_[i] = rhs->ColIndex()[i];
            MatVal_[i] = rhs->MatVal()[i];
            order_[i] = rhs->order()[i];
        }
        return true;

    }

    bool SpTripletMat::isSymmetric() const {
        return isSymmetric_;
    }

    bool SpTripletMat::transposed_times(std::shared_ptr<const Vector> p,
                                        std::shared_ptr<Vector> result) const {
        if (isSymmetric_) {
            times(p, result);
        } else {
            result->set_zeros(); //set all entries to be 0
            for (int i = 0; i < EntryNum_; i++) {
                result->addNumberAt(ColIndex_[i] - 1, MatVal_[i] * p->values()
                [RowIndex_[i] - 1]);
            }
        }//TODO:test it..

        return true;
    }

/**
 * qpOASESSparseMatrix
 */
    void qpOASESSparseMat::print() const {
        std::cout << "ColIndex: ";
        for (int i = 0; i < ColNum_ + 1; i++)
            std::cout << ColIndex()[i] << " ";

        std::cout << " " << std::endl;

        std::cout << "RowIndex: ";

        for (int i = 0; i < EntryNum_; i++)
            std::cout << RowIndex()[i] << " ";
        std::cout << " " << std::endl;

        std::cout << "MatVal:   ";

        for (int i = 0; i < EntryNum_; i++)
            std::cout << MatVal()[i] << " ";
        std::cout << " " << std::endl;

        std::cout << "order:    ";
        for (int i = 0; i < EntryNum_; i++)
            std::cout << order()[i] << " ";
        std::cout << " " << std::endl;


    }

/** Default constructor*/
    qpOASESSparseMat::qpOASESSparseMat(int RowNum, int ColNum, bool
    isSymmetric) :
            RowIndex_(NULL),
            ColIndex_(NULL),
            MatVal_(NULL),
            order_(NULL),
            EntryNum_(-1),
            isInitialised_(false),
            RowNum_(RowNum),
            ColNum_(ColNum),
            isSymmetric_(isSymmetric) {
        ColIndex_ = new qpOASES::int_t[ColNum + 1]();
    }

/**
 *@brief A constructor with the number of non-zero entries, row number and column
 * number specified.
 *
 * @param: nnz number of nonzero entry
 * @param RowNum: number of rows of a matrix
 * @param ColNum: number of columns of a matrix
 */
    qpOASESSparseMat::qpOASESSparseMat(int nnz, int RowNum, int ColNum) :
            RowIndex_(NULL),
            ColIndex_(NULL),
            MatVal_(NULL),
            order_(NULL),
            isSymmetric_(false),
            isInitialised_(false) {
        EntryNum_ = nnz;
        RowNum_ = RowNum;
        ColNum_ = ColNum;
        ColIndex_ = new qpOASES::int_t[ColNum + 1]();
        RowIndex_ = new qpOASES::int_t[nnz]();
        MatVal_ = new qpOASES::real_t[nnz]();
        order_ = new int[nnz]();
        for (int i = 0; i < nnz; i++)
            order_[i] = i;

    }

/**
 *Default destructor
 */
    qpOASESSparseMat::~qpOASESSparseMat() {
        freeMemory();
    }


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
    bool qpOASESSparseMat::setStructure(std::shared_ptr<const SpTripletMat> rhs,
                                        Identity2Info I_info) {
        assert(isInitialised_ == false);
        int counter = 0; // the counter for recording the index location
        std::vector<std::tuple<int, int, int>> sorted_index_info;
        for (int i = 0; i < rhs->EntryNum(); i++) {
            sorted_index_info.push_back(std::make_tuple(rhs->RowIndex()[i],
                                                        rhs->ColIndex()[i], counter));
            counter++;
        }

        // adding 2 identity matrix to the tuple array.
        if (I_info.irow1 != 0) {
            for (int j = 0; j < I_info.size; j++) {
                sorted_index_info.push_back(
                        std::make_tuple(I_info.irow1 + j, I_info.jcol1 + j, counter));
                sorted_index_info.push_back(
                        std::make_tuple(I_info.irow2 + j, I_info.jcol2 + j,
                                        counter + 1));
                counter += 2;
            }

        }

        assert(counter == EntryNum_);

        std::sort(sorted_index_info.begin(), sorted_index_info.end(),
                  tuple_sort_rule);
        //copy the order information back

        for (int i = 0; i < EntryNum_; i++) {
            RowIndex_[i] = std::get<0>(sorted_index_info[i]) - 1;
            order_[i] = std::get<2>(sorted_index_info[i]);
            if (i < EntryNum_ - 1) {
                if (std::get<1>(sorted_index_info[i]) <
                    std::get<1>(sorted_index_info[i + 1])) {
                    ColIndex_[std::get<1>(sorted_index_info[i])] = i + 1;
                }
            }

            int j = ColNum_;
            while (j >= 0 && ColIndex_[j] == 0) {
                ColIndex_[j] = EntryNum_;
                if (ColIndex_[j - 1] == EntryNum_)
                    break;
                else
                    j--;

            }
        }

        sorted_index_info.clear();
        return true;
    }


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
    bool qpOASESSparseMat::setStructure(std::shared_ptr<const SpTripletMat> rhs) {
        assert(isInitialised_ == false);
        int counter = 0; // the counter for recording the index location
        std::vector<std::tuple<int, int, int>> sorted_index_info;
        if (isSymmetric_) {
            for (int i = 0; i < rhs->EntryNum(); i++) {

                sorted_index_info.push_back(std::make_tuple(rhs->RowIndex()[i],
                                                            rhs->ColIndex()[i],
                                                            counter));
                counter++;
                if (rhs->RowIndex()[i] != rhs->ColIndex()[i]) {
                    sorted_index_info.push_back(std::make_tuple(rhs->ColIndex()[i],
                                                                rhs->RowIndex()[i],
                                                                counter));
                    counter++;
                }

            }
//            std::cout << "EntryNum is " << EntryNum_;
            if (EntryNum_ == -1) {
                EntryNum_ = counter;
                RowIndex_ = new qpOASES::sparse_int_t[counter]();
                MatVal_ = new double[counter]();
                order_ = new int[counter]();
            }
        }
        assert(counter == EntryNum_);

        std::sort(sorted_index_info.begin(), sorted_index_info.end(),
                  tuple_sort_rule);
        //copy the order information back

        for (int i = 0; i < EntryNum_; i++) {
            RowIndex_[i] = std::get<0>(sorted_index_info[i]) - 1;
            order_[i] = std::get<2>(sorted_index_info[i]);
            if (i < EntryNum_ - 1) {
                if (std::get<1>(sorted_index_info[i]) <
                    std::get<1>(sorted_index_info[i + 1])) {
                    ColIndex_[std::get<1>(sorted_index_info[i])] = i + 1;
                }
            }
            int j = ColNum_;

            while (j >= 1 && ColIndex_[j] == 0) {
                ColIndex_[j] = EntryNum_;
                if (ColIndex_[j - 1] == EntryNum_)
                    break;
                else
                    j--;

            }
        }

        sorted_index_info.clear();
        return true;
    }


/**
 * @brief set the Matrix values to the matrix, convert from triplet format to
 * Harwell-Boeing Matrix format.
 * @param MatVal entry values(orders are not yet under permutation)
 * @param I_info the 2 identity matrices information
 */
    bool qpOASESSparseMat::setMatVal(const double* MatVal, Identity2Info I_info) {
        //adding the value to the matrix
        if (isInitialised_ == false) {
            for (int i = 0; i < I_info.size; i++) {
                MatVal_[order()[EntryNum_ - i - 1]] = -1;
                MatVal_[order()[EntryNum_ - i - I_info.size - 1]] = 1;
            }

            isInitialised_ = true;
        }

        //assign each matrix entry to the corresponding
        //position after permutation
        for (int i = 0; i < EntryNum_ - 2 * I_info.size; i++) {
            MatVal_[order()[i]] = MatVal[i];
        }
        return true;
    }

    bool qpOASESSparseMat::setMatVal(std::shared_ptr<const SpTripletMat> rhs) {
        int j = 0;
        for (int i = 0; i < rhs->EntryNum(); i++) {
            MatVal_[order()[j]] = rhs->MatVal()[i];
            j++;
            if (isSymmetric_ && (rhs->ColIndex()[i] != rhs->RowIndex()[i])) {
                MatVal_[order()[j]] = rhs->MatVal()[i];
                j++;
            }
        }
//        print();
        return true;
    }

/**
* Free all memory allocated
*/
    bool qpOASESSparseMat::freeMemory() {
        delete[] ColIndex_;
        ColIndex_ = NULL;
        delete[] RowIndex_;
        RowIndex_ = NULL;
        delete[] MatVal_;
        MatVal_ = NULL;
        delete[] order_;
        order_ = NULL;
        return true;
    }

    bool qpOASESSparseMat::copy(std::shared_ptr<const qpOASESSparseMat> rhs) {
        assert(isInitialised_ == rhs->isIsinitialized());
        assert(EntryNum_ == rhs->EntryNum());
        assert(RowNum_ == rhs->RowNum());
        assert(ColNum_ == rhs->ColNum());
        for (int i = 0; i < EntryNum_; i++) {
            RowIndex_[i] = rhs->RowIndex()[i];
            MatVal_[i] = rhs->MatVal()[i];
            order_[i] = rhs->order()[i];
        }

        for (int i = 0; i < ColNum_ + 1; i++) {
            ColIndex_[i] = rhs->ColIndex()[i];
        }
        return true;
    }


}//END_OF_NAMESPACE

