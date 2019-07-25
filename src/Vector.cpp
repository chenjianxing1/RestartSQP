/* Copyright (C) 2019
 * All Rights Reserved.
 *
 * Authors: Xinyi Luo
 * Date:2019-07
 */
#include <sqphot/Vector.hpp>


namespace SQPhotstart {
    /**
     * Default constructor
     * @name initialize the size of the vector and allocate
     * the memory to the array
     */
    Vector::Vector(int vector_size)
    :
    size_(vector_size),
    values_(NULL)
    {
        values_ = new double[size_]();
    }
    
    /**
     * constructor that initializes the size of the vector
     * and initializes @vector_ to be @vector_value
     */
    
    Vector::Vector(int vector_size, const double *vector_value)
    :
    size_(vector_size),
    values_(NULL)
    {
        values_ = new double[size_]();
        std::copy(vector_value, vector_value+size_,values_);
    }
    
    /** Default destructor*/
    Vector::~Vector() {
        delete[] values_;
        values_ = NULL;
    }
    
    
    /** assign a sub-vector into the class member vector
     * without shifting elements' positions*/
    void Vector::assign(int Location, int subvector_size, const double *subvector) {
        for (int i = 0; i < subvector_size; i++) {
            values_[Location + i - 1] = subvector[i];
        }
    }
    
    
    /** assign a sub-vector with a specific size to the class member vector at a specific location.
     * The sub-vector will have all elements equal to the scaling factor*/
    void Vector::assign_n(int Location, int subvector_size, double scaling_factor) {
        for (int i = 0; i < subvector_size; i++) {
            values_[Location + i - 1] = scaling_factor;
        }
    }

    /** print the vector*/
    void Vector::print() const {
        for (int i = 0; i < size_; i++) std::cout << values_[i] << std::endl;
    }
    
    /* add all the element in the array by a number*/
    void Vector::addNumber(double increase_amount) {
        for (int i = 0; i < size_; i++) {
            values_[i] += increase_amount;
        }
    }
    
    /* add the element in a specific location by a number*/
    void Vector::addNumberAt(int Location, double increase_amount) {
        values_[Location] += increase_amount;
    }
    
    /* add all elements from the initial location to
     * the end location  specified by user by a number*/
    void Vector::add_number(int initloc, int endloc, double increase_amount) {
        for (int i = initloc - 1; i < endloc; i++) {
            values_[i] += increase_amount;
        }
    }
    
    /** add _vector with another vector, store the results in _vector*/
    void Vector::add_vector(const double *rhs) {
        for (int i = 0; i < size_; i++) {
            values_[i] += rhs[i];
        }
    }
    
    /** subtract a vector @rhs from the class member @_vector*/
    void Vector::subtract_vector(const double *rhs) {
        for (int i = 0; i < size_; i++) {
            values_[i] -= rhs[i];
        }
    }
    
    
    /**
     * subtract  a subvector with length @subvec_size from the class member @_vector
     * from the location @iloc
     *
     * @param iloc the starting location to subtract the subvector
     * @param subvec_size the size of the subvector
     */
    void Vector::subtract_subvector(int iloc, int subvec_size, const double *subvector) {
        for(int i = 0;i<subvec_size; i++ ){
            values_[i+iloc-1] -= subvector[i];
        }
    }
    
    /*copy all the entries from another vector*/
    void Vector::copy_vector(const double *rhs) {
        for (int i = 0; i < size_; i++) {
            values_[i] = (double) rhs[i];
        }
    }

    /*copy all the entries from another vector*/
    void Vector::copy_vector(std::shared_ptr<const Vector> rhs) {
        assert(size_ == rhs->Dim());
        for(int i=0; i<size_; i++){
            values_[i] = rhs->values()[i];
        }
//        print();
    }


    /** copy a subvector from member_vector from (Location) to (Location+subvector_size) to the pointer (results)*/
    void Vector::get_subVector(int Location, int subvector_size, std::shared_ptr<Vector> rhs) const {
        //TODO::test it! not sure if it will work...
        double *tmp = &(values_[Location - 1]);
        rhs->assign(1, subvector_size, tmp);
    }
    
    
    /** calculate one norm of the member _vector*/
    double Vector::getInfNorm() const {
        double infnorm = 0;
        for (int i = 0; i < size_; i++) {
            double absxk;
            if (values_[i] < 0) absxk = -values_[i];
            else absxk = values_[i];
            if (absxk > infnorm) infnorm = absxk;
        }
        return infnorm;
    }
    
    /** set all entries to be 0*/
    void Vector::set_zeros() {
        for (int i = 0; i < size_; i++) {
            values_[i] = 0;
        }
    }
    /**
     * @brief
     * @param rhs
     * @param compared_const
     * @return
     */
    bool Vector::setValue2Min(const Vector &rhs, double compared_const) {
        //        if(size_==0){
        //            size_ = rhs.size_;
        //            values_ = new double[size_];
        //        }
        //        for(int i=0; i<size_;i++)
        //            setValueAt(i,std::min(rhs->values()[i],compared_const));
        return true;
    }
    
    Vector Vector::operator+(const Vector& rhs) {
        assert(this->size_==rhs.size_);
        auto added_result = Vector(this->size_);
        for (int i = 0; i < this->size_; i++)
            added_result.setValueAt(i, this->values()[i] + rhs.values()[i]);
        return added_result;
    }
    
    const Vector Vector::operator-(const Vector &rhs) {
        assert(this->size_==rhs.size_);
        auto added_result = Vector(this->size_);
        for (int i = 0; i < this->size_; i++)
            added_result.setValueAt(i, this->values()[i] - rhs.values()[i]);
        return added_result;
    }
    
    /**
     * @brief set each entry of values_ to be the max of the same entry of rhs and a
     * constant number.
     *
     * For example, values_[i] = max(rhs.values[i],compared_const).
     *
     * @param rhs the vector to be compared with the constant
     * @param compared_const the constant to be compared with each entry of the vector
     */
    bool Vector::setValue2Max(const Vector rhs, double compared_const) {
        assert(this->size_==rhs.Dim());
        if(size_==0){
            size_ = rhs.size_;
            values_ = new double[size_]();
        }
        for(int i=0; i<size_;i++)
            setValueAt(i,std::max(rhs.values()[i],compared_const));
        return true;
    }

    double Vector::times(Vector& rhs) {
        double product = 0;
        for(int i = 0; i<size_;i++){
            product+=values_[i]*rhs.values()[i];
        }
        return product;
    }


}

