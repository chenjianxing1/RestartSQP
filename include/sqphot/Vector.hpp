/* Copyright (C) 2019
 * All Rights Reserved.
 *
 * Authors: Xinyi Luo
 * Date:2019-06
 */
#ifndef SQPHOTSTART_VECTOR_HPP
#define SQPHOTSTART_VECTOR_HPP

#include <memory>
#include <iostream>
#include <cassert>

namespace SQPhotstart {
    
    
    /** 
     * This is part of SQPhotstart.
     * @Vector This is a class for handling vectors and some operations that are used in the 
     * SQP algorithms
     */
    
    class Vector {
    public:
        /** Default Constructor*/
        Vector() :
        size_(0),
        values_(NULL) {};
        
        /**
         * Default constructor
         * @name initialize the size of the vector and allocate the memory to the array
         */
        Vector(int vector_size);
        
        /** 
         * constructor that initializes the _size and make a deep copy from vector_value to
         * class member _vector
         */
        Vector(int vector_size, const double* vector_value);
        
        /** Default destructor*/
        virtual ~Vector();
        
        
        /** assign a sub-vector into the class member vector without shifting elements' positions
         * */
        void assign(int Location, int subvector_size, const double* subvector);
        
        inline void setValueAt(int location, double value) {
            values_[location] = value;
        }
        
        /** assign a sub-vector with a specific size to the class member vector at a specific location.
         * The sub-vector will have all elements equal to the scaling factor
         * */
        void assign_n(int Location, int subvector_size, double scaling_factor);
        
        /** print the vector*/
        void print() const;
        
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
         * subtract  a subvector with length @subvec_size from the class member @_vector
         * from the location @iloc
         *
         * @param iloc the starting location to subtract the subvector
         * @param subvec_size the size of the subvector
         */
        void subtract_subvector(int iloc, int subvec_size, const double* subvector);
        
        /**copy all the entries from another vector*/
        void copy_vector(const double* copy);
        
        /** 
         * copy a subvector from member_vector from (Location) to (Location+subvector_size) 
         * to the pointer (results)
         */
        void get_subVector(int Location, int subvector_size, std::shared_ptr<Vector> rhs) const;
        
        /** set all entries to be 0*/
        void set_zeros();
        
        /** calculate one norm of the member _vector*/
        double getOneNorm() const;
        
        /** calculate the infinity norm of the member _vector*/
        double getInfNorm() const;
        
        /**
         * @brief
         * @param rhs
         * @return
         */
        Vector operator+(const Vector &rhs);
        
        /**
         * @brief
         * @param rhs
         * @return
         */
        const Vector operator-(const Vector &rhs);
        
        bool setValue2Max(const Vector rhs, double compared_const);
        
        bool setValue2Min(const Vector &rhs, double compared_const);
        /**
         * get the class member
         */
        
        inline int Dim() { return size_; }
        
        inline int Dim() const { return size_; }
        
        inline double* values() { return values_; }
        
        inline const double* values() const { return values_; }
        
        /**
         * @brief copy constructor
         * @param rhs
         */
        void operator = (const Vector& rhs){
            size_ = rhs.size_;
            values_ = rhs.values_;
        }
        
    private:
        
        //        Vector():
        
        //        /** Copy Constructor */
        //        Vector(const Vector &);
        //
        //        /** Overloaded Equals Operator */
        //        void operator=(const Vector &);
        
        int size_; /* the size of an vector*/
        double* values_;/*the array stored vector information*/
    };
    
}

#endif



