#include <sqphot/Matrix.hpp>

namespace SQPhotstart {
    
    /** Constructor for an empty matrix with N non-zero entries*/
    SpMatrix::SpMatrix(int nnz, int RowNum, int ColNum) :
    RowIndex_(NULL),
    ColIndex_(NULL),
    MatVal_(NULL),
    order_(NULL){
        EntryNum_ = nnz;
        RowNum_ = RowNum;
        ColNum_ = ColNum;
        //do nothing unless any data is to be assigned
        RowIndex_ = new int[nnz];
        ColIndex_ = new int[nnz];
        MatVal_ = new Number[nnz];
        order_ = new int[nnz];
        //initialize the order to 0:N-1
        for (int i = 0; i < nnz; i++) {
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
    
    
    
    /**
     * qpOASESSparseMatrix
     */
    void qpOASESSparseMat::print() {
        std::cout<<"ColIndex: ";
        for(int i=0;i<ColNum_+1; i++)
            std::cout<<ColIndex()[i]<<" ";
        
        std::cout<<" "<<std::endl;
        
        std::cout<<"RowIndex: ";
        
        for(int i =0; i<EntryNum_; i++)
            std::cout<<RowIndex()[i]<<" ";
        std::cout<<" "<<std::endl;
        
        std::cout<<"MatVal:   ";
        
        for(int i =0; i<EntryNum_; i++)
            std::cout<<MatVal()[i]<<" ";
        std::cout<<" "<<std::endl;
        
        std::cout<<"order:    ";
        for(int i =0; i<EntryNum_; i++)
            std::cout<<order()[i]<<" ";
        std::cout<<" "<<std::endl;
        
        
    }
    
    /**
     *
     */
    qpOASESSparseMat::~qpOASESSparseMat() { freeMemory(); }
    
    /**
     *
     *
     * @param nnz
     * @param RowNum
     * @param ColNum
     */
    qpOASESSparseMat::qpOASESSparseMat(int nnz, int RowNum, int ColNum) :
    RowIndex_(NULL),
    ColIndex_(NULL),
    MatVal_(NULL),
    order_(NULL){
        isinitialized = false;
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
    
    bool qpOASESSparseMat::setStructure(std::shared_ptr<const SpMatrix> rhs,
                                        Identity2Info I_info){
        assert(isinitialized==false);
        
        int counter = 0; // the counter for recording the index location
        std::vector<std::tuple<int, int, int>> sorted_index_info;
        for (int i = 0; i < rhs->EntryNum(); i++){
            sorted_index_info.push_back(std::make_tuple(rhs->RowIndex()[i],rhs->ColIndex()[i],counter));
            counter++;
        }
        
        // adding 2 identity matrix to the tuple array.
        if(I_info.irow1!=0){
            for (int j=0; j<I_info.size;j++){
                sorted_index_info.push_back(std::make_tuple(I_info.irow1+j,I_info.jcol1+j,counter));
                sorted_index_info.push_back(std::make_tuple(I_info.irow2+j,I_info.jcol2+j,counter+1));
                counter +=2;
            }
            
        }
        
        assert(counter==EntryNum_);
        
        std::sort(sorted_index_info.begin(),sorted_index_info.end(), tuple_sort_rule);
        //copy the order information back
        
        for (int i = 0; i <EntryNum_; i++) {
            RowIndex_[i] = std::get<0>(sorted_index_info[i])-1;
            order_[i] = std::get<2>(sorted_index_info[i]);
            if(i<EntryNum_-1){
                if(std::get<1>(sorted_index_info[i])<std::get<1>(sorted_index_info[i+1])){
                    ColIndex_[std::get<1>(sorted_index_info[i])] = i+1;
                }
            }
            
            int j = ColNum_;
            while(j>=0 && ColIndex_[j] ==0){
                ColIndex_[j] = EntryNum_;
                if(ColIndex_[j-1] ==EntryNum_)
                    break;
                else
                    j--;
                
            }
        }
        
        sorted_index_info.clear();
        return true;
    }
    
    
    bool qpOASESSparseMat::setMatVal(const double *MatVal, Identity2Info I_info){
        //adding the value to the matrix
        if(isinitialized == false){
            for(int i = 0; i<I_info.size; i++){
                MatVal_[order()[EntryNum_-i-1]] = -1;
                MatVal_[order()[EntryNum_-i-I_info.size-1]]=1;
            }
            
            isinitialized = true;
        }
        
        //assign each matrix entry to the corresponding position after
        //permutaion
        for(int i=0; i<EntryNum_-2*I_info.size; i++){
            MatVal_[order()[i]] = MatVal[i];
        }
        return true;
    }
    
    bool qpOASESSparseMat::addIdentityMat(int irow, int jcol, int size, double scaling_factor){
        
        return true;
    }
    
    
    bool qpOASESSparseMat::freeMemory(){
        delete[] ColIndex_;
        ColIndex_=NULL;
        delete[] RowIndex_;
        RowIndex_=NULL;
        delete[] MatVal_;
        MatVal_ =NULL;
        delete[] order_;
        order_ = NULL;
        return true;
    }
    
    
    
}

