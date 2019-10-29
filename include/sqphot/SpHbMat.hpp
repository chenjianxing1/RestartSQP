
/* Copyright (C) 2019
 * All Rights Reserved.
 *
 * Authors: Xinyi Luo
 * Date:2019-10
 */

#ifndef SQPHOTSTART_SPHBMAT_HPP_
#define SQPHOTSTART_SPHBMAT_HPP_

#include <sqphot/SpTripletMat.hpp>

namespace SQPhotstart {

//forward declaration

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
    SpHbMat(int nnz, int RowNum, int ColNum, bool isCompressedRow);


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

//@{
    inline void setMatValAt(int i, double value) {
        assert(i<EntryNum_);
        MatVal_[i] = value;
    }

    inline void setRowIndexAt(int i, int value) {
        if(isCompressedRow_) {
            assert(i<RowNum_+1);
        }
        else {
            assert(i<EntryNum_);
        }

        RowIndex_[i] = value;
    }

    inline void setColIndexAt(int i, int value) {
        if(isCompressedRow_) {
            assert(i<EntryNum_);
        }
        else {
            assert(i<ColNum_+1);
        }

        ColIndex_[i] = value;
    }

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

    void get_dense_matrix(double* dense_matrix,bool row_oriented = true) const ;

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

    void transposed_times(const double* p, double* result) const {
        for(int i = 0; i<ColNum_; i++)
            result[i] = 0.0;

        if(isCompressedRow_) {
            int row;
            for(int i = 1; i <RowNum_+1; i++) {
                if(RowIndex_[i]>0) {
                    row = i-1;
                    break;
                }
            }
            for(int i = 0; i<EntryNum_; i++) {
                while(i==RowIndex_[row+1]) {
                    row++;
                }
                result[ColIndex_[i]] +=MatVal_[i]*p[row];
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
                result[col] += MatVal_[i]*p[RowIndex_[i]];
            }
        }

    }
    /**
     * @brief make a deep copy of a matrix information
     */

    virtual  void copy(std::shared_ptr<const SpHbMat> rhs);



    const double oneNorm() const ;

    const double infNorm() const ;
    /**
     * @brief convert the matrix data stored in the class members to a triplet matrix
     * speci fied by rhs */
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

    inline const int RowIndex(int i) const {

        return RowIndex_[i];
    }

    inline const int ColIndex(int i ) const {

        return ColIndex_[i];
    }


    inline const double MatVal(int i)  const {

        return MatVal_[i];
    }


    inline const int order(int i)  const {

        return order_[i];
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


//    /** Copy Constructor */
//    SpHbMat(const SpHbMat &);


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

#endif
