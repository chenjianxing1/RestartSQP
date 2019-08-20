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
        assert(isAllocated() == false);
        this->values_ = rhs;
        printf("rhs is\n");
        for (int i = 0; i < size_; i++) {
            std::cout << rhs[i] << std::endl;
        }
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
        if (size == 0) {
            values_ = new double[size_]();
        } else {
            assert(size == size_);
            values_ = new double[size]();
        }
    }

    /** Default destructor*/
    virtual ~Vector();


    /** assign a sub-vector into the class member vector without shifting elements' positions
     * */
    void assign(int Location, int subvector_size, const double* subvector);

    inline void setValueAt(int location, double value) {
        if(value>INF)
            values_[location] = INF;
        else
            values_[location] = value;

    }

    /** assign a sub-vector with a specific size to the class member vector at a specific location.
     * The sub-vector will have all elements equal to the scaling factor
     * */
    void assign_n(int Location, int subvector_size, double scaling_factor);

    /** print the vector*/
    void print(const char* name = NULL) const;

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

    /**
     *
     * subtract  a subvector with length @subvec_size from the class member @_vector
     * from the location @iloc
     *
     * @param iloc the starting location to subtract the subvector
     * @param subvec_size the size of the subvector
     */
    void subtract_subvector(int iloc, int subvec_size, const double* subvector);

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
    double getOneNorm() const;

    /** calculate the infinity norm of the member _vector*/
    double getInfNorm() const;

    double times(std::shared_ptr<Vector> rhs);


    const double* negative_of_values() const {
        std::shared_ptr<Vector> negative_values = std::make_shared<Vector>(size_);
        for (int i = 0; i < size_; i++)
            negative_values->values()[i] = -values_[i];


        return negative_values->values();
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

    /**
     * @brief copy constructor
     * @param rhs
     */
    void operator=(const Vector &rhs) {
        size_ = rhs.size_;
        values_ = rhs.values_;
    }

    void write_to_file(FILE* file_to_write, const char* const name) {
        fprintf(file_to_write, "real_t %s[] = {", name);
        for (int i = 0; i < Dim(); i++) {
            if (i % 10 == 0 && i > 1)
                fprintf(file_to_write, "\n");
            if (i == Dim() - 1)
                fprintf(file_to_write, "%10e};\n\n", values_[i]);
            else
                fprintf(file_to_write, "%10e, ", values_[i]);
        }
    }

    bool isAllocated() {
        return isAllocated_;
    }

private:

    /** Default Constructor*/
    Vector();

    //        /** Copy Constructor */
    //        Vector(const Vector &);
    //
    //        /** Overloaded Equals Operator */
    //        void operator=(const Vector &);
private:
    bool isAllocated_;
    int size_; /* the size of an vector*/
    double* values_;/*the array stored vector information*/
};

}

#endif



