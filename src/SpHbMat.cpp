#include <sqphot/Matrix.hpp>


namespace SQPhotstart {

/** Default constructor*/
SpHbMat::SpHbMat(int RowNum, int ColNum, bool isSymmetric, bool isCompressedRow) :
    RowIndex_(NULL),
    ColIndex_(NULL),
    MatVal_(NULL),
    order_(NULL),
    EntryNum_(-1),
    isInitialised_(false),
    RowNum_(RowNum),
    ColNum_(ColNum),
    isSymmetric_(isSymmetric),
    isCompressedRow_(isCompressedRow) {

    if(isCompressedRow_) {
        RowIndex_ = new int[RowNum + 1]();
    } else
        ColIndex_ = new int[ColNum + 1]();
}


/**
 *@brief A constructor with the number of non-zero entries, row number and column
 * number specified.
 *
 * @param: nnz number of nonzero entry
 * @param RowNum: number of rows of a matrix
 * @param ColNum: number of columns of a matrix
 */
SpHbMat::SpHbMat(int nnz, int RowNum, int ColNum, bool isCompressedRow) :
    RowIndex_(NULL),
    ColIndex_(NULL),
    MatVal_(NULL),
    order_(NULL),
    EntryNum_(nnz),
    RowNum_(RowNum),
    ColNum_(ColNum),
    isSymmetric_(false),
    isInitialised_(false),
    isCompressedRow_(isCompressedRow) {

    if(isCompressedRow) {
        ColIndex_ = new int[nnz]();
        RowIndex_ = new int[RowNum + 1]();
    } else {
        ColIndex_ = new int[ColNum + 1]();
        RowIndex_ = new int[nnz]();
    }
    MatVal_ = new double[nnz]();
    order_ = new int[nnz]();
    for (int i = 0; i < nnz; i++)
        order_[i] = i;
}


SpHbMat::SpHbMat(double* data, int RowNum, int ColNum, bool row_oriented,
                 bool isCompressedRow):
    RowIndex_(NULL),
    ColIndex_(NULL),
    MatVal_(NULL),
    order_(NULL),
    EntryNum_(0),
    RowNum_(RowNum),
    ColNum_(ColNum),
    isCompressedRow_(isCompressedRow)
{

    int* RowIndex_tmp = NULL;
    int* ColIndex_tmp = NULL;
    double* MatVal_tmp = NULL ;
    //allocate memory
    if(isCompressedRow) {
        RowIndex_ = new int[RowNum + 1]();
        ColIndex_tmp = new int[RowNum*ColNum]();
        MatVal_tmp = new double[RowNum*ColNum]();
    } else {
        ColIndex_ = new int[ColNum + 1]();
        RowIndex_tmp = new int[RowNum*ColNum]();
        MatVal_tmp = new double[RowNum*ColNum]();
    }

    EntryNum_ = 0;
    if(row_oriented) {
        if(isCompressedRow) {
            for(int i = 0; i<RowNum; i++) {
                for(int j = 0; j<ColNum; j++) {
                    //identify nonzero entry
                    if(data[i]>m_eps||data[i]<-m_eps) {
                        MatVal_tmp[EntryNum_] = data[i*ColNum+j];
                        ColIndex_tmp[EntryNum_] = j;
                        EntryNum_++;
                    }
                }
                RowIndex_[i+1] = EntryNum_;
            }
        }
        else { //if it is condensed column
            for(int j = 0; j<ColNum; j++) {
                for(int i = 0; i<RowNum; i++) {
                    if(data[j+i*ColNum]>m_eps||data[j+i*ColNum]<-m_eps) {
                        MatVal_tmp[EntryNum_] = data[j+i*ColNum];
                        RowIndex_tmp[EntryNum_] = i;
                        EntryNum_++;
                    }
                }
                ColIndex_[j+1] = EntryNum_;
            }
        }
    }
    else {
        if(isCompressedRow) {
            for(int j = 0; j<RowNum; j++) {
                for(int i = 0; i<ColNum; i++) {
                    if(data[j+i*RowNum]>m_eps||data[j+i*RowNum]<-m_eps) {
                        MatVal_tmp[EntryNum_] = data[j+i*RowNum];
                        ColIndex_tmp[EntryNum_] = i;
                        EntryNum_++;
                    }
                }
                RowIndex_[j+1] = EntryNum_;
            }
        }
        else {
            for(int i = 0; i<ColNum; i++) {
                for(int j = 0; j<RowNum; j++) {
                    //identify nonzero entry
                    if(data[i]>m_eps||data[i]<-m_eps) {
                        MatVal_tmp[EntryNum_] = data[i*RowNum+j];
                        RowIndex_tmp[EntryNum_] = j;
                        EntryNum_++;
                    }
                }
                ColIndex_[i+1] = EntryNum_;
            }
        }
    }
    //allocate memory for class member
    MatVal_ = new double[EntryNum_]();
    order_ = new int[EntryNum_]();
    if(isCompressedRow) {
        ColIndex_ = new int[EntryNum_]();
        for(int i = 0; i<EntryNum_; i++) {
            ColIndex_[i] = ColIndex_tmp[i];
            MatVal_[i] = MatVal_tmp[i];
        }
        delete[] MatVal_tmp;
        MatVal_tmp = NULL;
        delete[] ColIndex_tmp;
        ColIndex_tmp = NULL;

    } else {
        RowIndex_ = new int[EntryNum_]();
        for(int i = 0; i<EntryNum_; i++) {
            RowIndex_[i] = RowIndex_tmp[i];
            MatVal_[i] = MatVal_tmp[i];
        }
        delete[] MatVal_tmp;
        MatVal_tmp = NULL;
        delete[] RowIndex_tmp;
        RowIndex_tmp = NULL;
    }
    order_ = new int[EntryNum()]();
}



/**
 *Default destructor
 */
SpHbMat::~SpHbMat() {

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
void SpHbMat::setStructure(std::shared_ptr<const SpTripletMat> rhs,
                           Identity2Info I_info) {


    set_zero();
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

    //  for(int i = 0; i <sorted_index_info.size(); i++) {
    //      printf("%i ", std::get<0>(sorted_index_info[i]));
    //      printf("%i ", std::get<1>(sorted_index_info[i]));
    //      printf("%i \n", std::get<2>(sorted_index_info[i]));
    //  }

    assert(counter == EntryNum_);

    if(isCompressedRow_) {
        std::sort(sorted_index_info.begin(), sorted_index_info.end(),
                  tuple_sort_rule_compressed_row);
        for (int i = 0; i < EntryNum_; i++) {
            ColIndex_[i] = std::get<1>(sorted_index_info[i]) - 1;
            order_[std::get<2>(sorted_index_info[i])]= i;

            for (int j = std::get<0>(sorted_index_info[i]); j < RowNum_; j++) {
                RowIndex_[j]++;
            }
        }
        RowIndex_[RowNum_] = EntryNum_;

    } else {
        std::sort(sorted_index_info.begin(), sorted_index_info.end(),
                  tuple_sort_rule_compressed_column);
        for (int i = 0; i < EntryNum_; i++) {
            RowIndex_[i] = std::get<0>(sorted_index_info[i]) - 1;
            order_[std::get<2>(sorted_index_info[i])] =i;

            for (int j = std::get<1>(sorted_index_info[i]); j < ColNum_; j++) {
                ColIndex_[j]++;
            }
        }
        ColIndex_[ColNum_] = EntryNum_;
    }
    sorted_index_info.clear();
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
 *
 */
void SpHbMat::setStructure(std::shared_ptr<const SpTripletMat> rhs) {

    set_zero();
    assert(isInitialised_ == false);
    int counter = 0; // the counter for recording the index location
    std::vector<std::tuple<int, int, int>> sorted_index_info;

    //if it is symmetric, it will calculate the
    // number of entry and allocate_memory the memory
    // of RowIndex and MatVal by going through
    // all of its entries one by one
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
        if (EntryNum_ == -1) {

            EntryNum_ = counter;
            MatVal_ = new double[counter]();
            order_ = new int[counter]();

            if (isCompressedRow_) {
                ColIndex_ = new int[counter]();
            } else {
                RowIndex_ = new int[counter]();
            }
        }
    }
    assert(counter == EntryNum_);


    if (isCompressedRow_) {
        std::sort(sorted_index_info.begin(), sorted_index_info.end(),
                  tuple_sort_rule_compressed_row);
        //copy the order information back

        for (int i = 0; i < EntryNum_; i++) {
            ColIndex_[i] = std::get<1>(sorted_index_info[i]) - 1;
            order_[std::get<2>(sorted_index_info[i])]= i;
            for (int j = std::get<0>(sorted_index_info[i]); j < RowNum_; j++) {
                RowIndex_[j]++;
            }
        }
        RowIndex_[RowNum_] = EntryNum_;
    } else {
        std::sort(sorted_index_info.begin(), sorted_index_info.end(),
                  tuple_sort_rule_compressed_column);
        //copy the order information back

        for (int i = 0; i < EntryNum_; i++) {
            RowIndex_[i] = std::get<0>(sorted_index_info[i]) - 1;
            order_[std::get<2>(sorted_index_info[i])] = i;

            for (int j = std::get<1>(sorted_index_info[i]); j < ColNum_; j++) {
                ColIndex_[j]++;
            }
        }
        ColIndex_[ColNum_] = EntryNum_;

    }
    sorted_index_info.clear();
}



/**
 * @brief set the Matrix values to the matrix, convert from triplet format to
 * Harwell-Boeing Matrix format.
 * @param rhs entry values(orders are not yet under permutation)
 * @param I_info the 2 identity matrices information
 */
void SpHbMat::setMatVal(std::shared_ptr<const SpTripletMat> rhs,
                        Identity2Info I_info) {                       //adding the value to the matrix
    if(isCompressedRow_) {
        if (isInitialised_ == false) {
            int i = rhs->EntryNum();
            while(i<EntryNum_) {
                MatVal_[order_[i]] = 1;
                MatVal_[order_[i+1]] = -1;
                i +=2;
            }
            isInitialised_ = true;
        }
    }
    else {
        if (isInitialised_ == false) {
            for (int i = 0; i < I_info.size; i++) {
                MatVal_[EntryNum_ - i - 1] = -1;
                MatVal_[EntryNum_ - i - I_info.size - 1] = 1;
            }
            isInitialised_ = true;
        }

    }
    //assign each matrix entry to the corresponding position after permutation
    for (int i = 0; i < EntryNum_ - 2 * I_info.size; i++) {
        MatVal_[order()[i]] = rhs->MatVal()[i];
    }
}


void SpHbMat::setMatVal(std::shared_ptr<const SpTripletMat> rhs) {

    int j = 0;
    for (int i = 0; i < rhs->EntryNum(); i++) {
        MatVal_[order_[j]] = rhs->MatVal()[i];
        j++;
        if (isSymmetric_ && (rhs->ColIndex()[i] != rhs->RowIndex()[i])) {
            //            MatVal_[order()[j]] = rhs->MatVal()[i];
            MatVal_[order_[j]] = rhs->MatVal()[i];
            j++;
        }
    }
}


/**
 * Free all memory allocated
 */
void SpHbMat::freeMemory() {
    delete[] ColIndex_;
    ColIndex_ = NULL;
    delete[] RowIndex_;
    RowIndex_ = NULL;
    delete[] MatVal_;
    MatVal_ = NULL;
    delete[] order_;
    order_ = NULL;
}


void SpHbMat::copy(std::shared_ptr<const SpHbMat> rhs) {

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
}


void SpHbMat::write_to_file(const char* name,
                            Ipopt::SmartPtr<Ipopt::Journalist> jnlst,
                            Ipopt::EJournalLevel level,
                            Ipopt::EJournalCategory category,
                            Solver solver) {
#if DEBUG
#if PRINT_QP_IN_CPP
    const char* var_type_int;
    const char* var_type_double;
    var_type_int = (solver == QPOASES) ? "sparse_int_t" : "qp_int";
    var_type_double = (solver == QPOASES) ? "real_t" : "double";
    jnlst->Printf(level, category, "%s %s_jc[] = \n{", var_type_int, name);
    int i;
    for (i = 0; i < ColNum_ + 1; i++) {
        if (i % 10 == 0 && i > 1)
            jnlst->Printf(level, category, "\n");
        if (i == ColNum_)
            jnlst->Printf(level, category, "%d};\n\n", ColIndex_[i]);
        else
            jnlst->Printf(level, category, "%d, ", ColIndex_[i]);
    }
    jnlst->Printf(level, category, "%s %s_ir[] = \n{", var_type_int, name);
    for (i = 0; i < EntryNum_; i++) {
        if (i % 10 == 0 && i > 1)
            jnlst->Printf(level, category, "\n");
        if (i == EntryNum_ - 1)
            jnlst->Printf(level, category, "%d};\n\n", RowIndex_[i]);
        else
            jnlst->Printf(level, category, "%d, ", RowIndex_[i]);
    }
    jnlst->Printf(level, category, "%s %s_val[] = \n{", var_type_double, name);
    for (i = 0; i < EntryNum_; i++) {
        if (i % 10 == 0 && i > 1)
            jnlst->Printf(level, category, "\n");
        if (i == EntryNum_ - 1)
            jnlst->Printf(level, category, "%23.16e};\n\n", MatVal_[i]);
        else
            jnlst->Printf(level, category, "%23.16e, ", MatVal_[i]);
    }
#else
    int i;
    for (i = 0; i < ColNum_ + 1; i++) {
        jnlst->Printf(level, category, "%d\n", ColIndex_[i]);
    }
    for (i = 0; i < EntryNum_; i++) {
        jnlst->Printf(level, category, "%d\n", RowIndex_[i]);
    }
    for (i = 0; i < EntryNum_; i++) {
        jnlst->Printf(level, category, "%23.16e\n", MatVal_[i]);
    }

#endif
#endif
}




void SpHbMat::get_dense_matrix(double* dense_matrix,bool row_oriented) {

    int row;
    if(row_oriented) {
        if(isCompressedRow_) {
            for(int i = 1; i <RowNum_+1; i++) {
                if(RowIndex_[i]>0) {
                    row = i-1;
                    break;
                }
            }
            for (int i = 0; i < EntryNum_; i++) {
                while(i==RowIndex_[row+1]) {
                    row++;
                }
                dense_matrix[ColNum_ * row + ColIndex_[i]] = MatVal_[i];
            }
        }
        else {
            int col;
            for(int i = 1; i <ColNum_+1; i++) {
                if(ColIndex_[i]>0) {
                    col = i-1;
                    break;
                }
            }
            for (int i = 0; i < EntryNum_; i++) {
                while(i==ColIndex_[col+1]) {
                    col++;
                }

                dense_matrix[ColNum_ * RowIndex_[i]+col] = MatVal_[i];
//                    printf("dense_matrix[%d] = %10e\n",ColNum_ * RowIndex_[i]+col, MatVal_[i]);
            }
        }
    }
    else {
        if(isCompressedRow_) {
            for(int i = 1; i <RowNum_+1; i++) {
                if(RowIndex_[i]>0) {
                    row = i-1;
                    break;
                }
            }
            for (int i = 0; i < EntryNum_; i++) {
                while(i==RowIndex_[row+1]) {
                    row++;
                }
                dense_matrix[RowNum_ * ColIndex_[i] + row] = MatVal_[i];
            }
        }
        else {
            int col;
            for(int i = 1; i <ColNum_+1; i++) {
                if(ColIndex_[i]>0) {
                    col = i-1;
                    break;
                }
            }
            for (int i = 0; i < EntryNum_; i++) {
                while(i==ColIndex_[col+1]) {
                    col++;
                }

                dense_matrix[RowNum_ * col+RowIndex_[i]] = MatVal_[i];
            }
        }

    }
}


void SpHbMat::transposed_times(shared_ptr<const Vector> p,
                               shared_ptr<Vector> result) const {

    result->set_zeros();

    int row;
    if(isCompressedRow_) {
        for(int i = 1; i <RowNum_+1; i++) {
            if(RowIndex_[i]>0) {
                row = i-1;
                break;
            }
        }
        for(int i = 0; i<EntryNum_; i++) {
            if(i==RowIndex_[row+1]) {
                row++;
            }
            result->addNumberAt(row, MatVal_[i]*p->values(ColIndex_[i]));
        }
    }
    else {


    }


}

void SpHbMat::times(std::shared_ptr<const Vector> p,
                    std::shared_ptr<Vector> result) const {

    result->set_zeros();

    if(isCompressedRow_) {
        int row;
        //find the row corresponding to the first nonzero entry
        for(int i = 1; i <RowNum_+1; i++) {
            if(RowIndex_[i]>0) {
                row = i-1;
                break;
            }
        }
        for(int i = 0; i<EntryNum_; i++) {
            //go to the next row
            while(i==RowIndex_[row+1]) {
                row++;
            }
            result->addNumberAt(row, MatVal_[i]*p->values(ColIndex_[i]));
        }
    }
    else {
        int col;
        //find the col corresponding to the first nonzero entry
        for(int i = 1; i<ColNum_+1; i++) {
            if(ColIndex_[i]>0) {
                col = i-1;
                break;
            }
        }
        for(int i = 0; i<EntryNum_; i++) {
            //go to the next col
            while(i == ColIndex_[col+1]) {
                col++;
            }
            result->addNumberAt(RowIndex_[i], MatVal_[i]*p->values(col));
        }
    }
}


//@{
void
SpHbMat::print_full(const char* name, Ipopt::SmartPtr<Ipopt::Journalist> jnlst,
                    Ipopt::EJournalLevel level,
                    Ipopt::EJournalCategory category) const {
    char mat_val[99];
    auto dense_matrix = new double[RowNum_ * ColNum_]();
    int row;
    if(isCompressedRow_) {
        for(int i = 1; i <RowNum_+1; i++) {
            if(RowIndex_[i]>0) {
                row = i-1;
                break;
            }
        }
        for (int i = 0; i < EntryNum_; i++) {
            if(i==RowIndex_[row+1]) {
                row++;
            }
            dense_matrix[ColNum_ * row + ColIndex_[i]] = MatVal_[i];
        }
    }
    else {
        int col;
        for(int i = 1; i <ColNum_+1; i++) {
            if(ColIndex_[i]>0) {
                col = i-1;
                break;
            }
        }
        for (int i = 0; i < EntryNum_; i++) {
            if(i==ColIndex_[col+1]) {
                col++;
            }
            dense_matrix[ColNum_ * RowIndex_[i]+col] = MatVal_[i];
        }
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
                printf("%23.16e  ", dense_matrix[i * ColNum() + j]);
            }
            printf("\n");
        }
        printf("}\n\n");
    }
    delete[] dense_matrix;
}


void SpHbMat::print(const char* name, Ipopt::SmartPtr<Ipopt::Journalist> jnlst,
                    Ipopt::EJournalLevel level,
                    Ipopt::EJournalCategory category) const {

    if(isCompressedRow_) {
        std::cout<< name <<"= "<<std::endl;
        std::cout << "ColIndex: ";
        for (int i = 0; i < EntryNum_; i++)
            std::cout << ColIndex()[i] << " ";
        std::cout << " " << std::endl;

        std::cout << "RowIndex: ";
        for (int i = 0; i < RowNum_+1; i++)
            std::cout << RowIndex()[i] << " ";
        std::cout << " " << std::endl;
    } else {
        //for compressed column format
        std::cout<< name <<"= "<<std::endl;
        std::cout << "ColIndex: ";
        for (int i = 0; i < ColNum_ + 1; i++)
            std::cout << ColIndex()[i] << " ";

        std::cout << " " << std::endl;
        std::cout << "RowIndex: ";

        for (int i = 0; i < EntryNum_; i++)
            std::cout << RowIndex()[i] << " ";
        std::cout << " " << std::endl;

    }
    std::cout << "MatVal:   ";

    for (int i = 0; i < EntryNum_; i++)
        std::cout << MatVal()[i] << " ";
    std::cout << " " << std::endl;

    std::cout << "order:    ";
    for (int i = 0; i < EntryNum_; i++)
        std::cout << order()[i] << " ";
    std::cout << " " << std::endl;
}

//@}
}//END_OF_NAMESPACE


