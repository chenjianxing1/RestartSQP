/* Copyright (C) 2019
 * All Rights Reserved.
 *
 * Authors: Xinyi Luo
 * Date:2019-06
 */
#ifndef SQPHOTSTART_VECTOR_HPP
#define SQPHOTSTART_VECTOR_HPP

#include <sqphot/Utils.hpp>

namespace SQPhotstart {


/**
 * This is part of SQPhotstart.
 * @Vector This is a class for handling vectors and some operations that are used in the
 * SQP algorithms
 */

class Vector {
public:

    /**
     * Default constructor
     * @name initialize the size of the vector and allocate the memory to the array
     */
    Vector(int vector_size, bool allocate = true);


    void swp(double* rhs) {
        assert(!isAllocated());
        this->values_ = rhs;
    }

    void free() {
        delete[] values_;
        values_ = NULL;
        isAllocated_ = false;
    }

    /**
     * constructor that initializes the _size and make a deep copy from vector_value to
     * class member _vector
     */
    Vector(int vector_size, const double* vector_value);


    void allocate_memory(int size = 0) {
        isAllocated_ = true;

        if (size != 0) {
            assert(size_ == size);
            values_ = new double[size]();
        } else {
            values_ = new double[size_]();
        }
    }

    /** Default destructor*/
    virtual ~Vector();


    /** assign a sub-vector into the class member vector without shifting elements' positions
     * */
    void assign(int Location, int subvector_size, const double* subvector);


    /** assign a sub-vector with a specific size to the class member vector at a specific location.
     * The sub-vector will have all elements equal to the scaling factor
     * */
    void assign_n(int Location, int subvector_size, double scaling_factor);

    /** print the vector*/
    void print(const char* name, Ipopt::SmartPtr<Ipopt::Journalist> jnlst= nullptr,
               Ipopt::EJournalLevel level=Ipopt::J_ALL, Ipopt::EJournalCategory
               category =Ipopt::J_DBG) const;


    /* add all the element in the array by a number*/
    void addNumber(double increase_amount);

    /* add the element in a specific location by a number*/
    void addNumberAt(int Location, double increase_amount);

    /** add all elements from the initial location to the end location  specified
     * by user by a number*/
    void add_number(int initloc, int endloc, double increase_amount);

    /** add _vector with another vector, store the results in _vector*/
    void add_vector(const double* rhs);

    /** subtract a vector @rhs from the class member @_vector*/
    void subtract_vector(const double* rhs);

    /** subtract class member @vector_ to a vector @rhs, modified @vector_ to be the result*/
    void subtract_vector_to(const double* rhs);

    /**
     *
     * subtract  a subvector with length @subvec_size from the class member @_vector
     * from the location @iloc
     *
     * @param begin  the starting location to subtract the subvector
     * @param subvec_size the size of the subvector
     */
    void subtract_subvector(int begin, int subvec_size, const double* subvector);

    /**copy all the entries from another vector*/
    void copy_vector(const double* rhs);

    /**copy all the entries from another vector*/
    void copy_vector(std::shared_ptr<const Vector> rhs);

    /**
     * copy a subvector from member_vector from (Location) to (Location+subvector_size)
     * to the pointer (results)
     */
    void get_subVector(int Location, int subvector_size,
                       std::shared_ptr<Vector> rhs) const;

    /** set all entries to be 0*/
    void set_zeros();

    /** calculate one norm of the member _vector*/
    double getOneNorm() const;;

    /** calculate the infinity norm of the member _vector*/
    double getInfNorm() const;

    double times(std::shared_ptr<Vector> rhs);



    void scale(double scaling_factor) {
        for(int i=0; i<size_; i++) {
            values_[i]*=scaling_factor;
        }
    }

    /**
     * get the class member
     */

    inline int Dim() {
        return size_;
    }

    inline int Dim() const {
        return size_;
    }

    inline double* values() {
        return values_;
    }



    inline const double* values() const {
        return values_;
    }

    inline const double values(int i) const {
        return values_[i];
    }
    inline void setValueAt(int location, double value) {
        values_[location] = value;
    }

//        /**
//         * @brief overload equal operator
//         * @param rhs
//         */
//        void operator=(const Vector &rhs) {
//            size_ = rhs.size_;
//            values_ = rhs.values_;
//        }
    void write_to_file(const char* name,
                       Ipopt::SmartPtr<Ipopt::Journalist> jnlst,
                       Ipopt::EJournalLevel level,
                       Ipopt::EJournalCategory category,
                       Solver qpsolver);


    bool isAllocated() {
        return isAllocated_;
    }

//    double operator[int i](const Vector&){
//        return values_[i];
//    }
private:


    /** Default Constructor*/
    Vector();

    /** Copy Constructor */
    Vector(const Vector &);

    /** Overloaded Equals Operator */
    void operator=(const Vector &);

private:
    bool isAllocated_;
    int size_; /* the size of an vector*/
    double* values_;/*the array stored vector information*/
};

}

#endif



