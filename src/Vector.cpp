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
Vector::Vector(int vector_size, bool allocate)
    :
    isAllocated_(allocate),
    size_(vector_size),
    values_(NULL) {
    if (allocate)
        allocate_memory();
}


/**
 * constructor that initializes the size of the vector
 * and initializes @vector_ to be @vector_value
 */

Vector::Vector(int vector_size, const double* vector_value)
    :
    isAllocated_(true),
    size_(vector_size),
    values_(NULL) {
    allocate_memory();
    std::copy(vector_value, vector_value + size_, values_);
}

/** Default destructor*/
Vector::~Vector() {
    if (isAllocated())
        free();
}


/** assign a sub-vector into the class member vector
 * without shifting elements' positions*/
void Vector::assign(int Location, int subvector_size, const double* subvector) {
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
void Vector::print(const char* name, Ipopt::SmartPtr<Ipopt::Journalist> jnlst,
                   Ipopt::EJournalLevel level, Ipopt::EJournalCategory category) const {

    if(!Ipopt::IsValid(jnlst)) {
        if (name != nullptr)
            printf("%s = :\n",name);

        for (int i = 0; i < size_; i++) {
            if (i == 0)
                printf("{ ");
            printf("%23.16e ", values_[i]);
        }
        printf("}\n\n");
    }
    else {
        if (name != nullptr) {
            jnlst->Printf(level, category, name);
            jnlst->Printf(level, category, " =: \n");
        }
        for (int i = 0; i < size_; i++) {
            if (i == 0)
                jnlst->Printf(level, category, "{ ");
            jnlst->Printf(level, category, "%23.16e ", values_[i]);
        }
        jnlst->Printf(level, category, "}\n\n");
    }
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
void Vector::add_vector(const double* rhs) {
    for (int i = 0; i < size_; i++) {
        values_[i] += rhs[i];
    }
}

/** subtract a vector @rhs from the class member @_vector*/
void Vector::subtract_vector(const double* rhs) {
    for (int i = 0; i < size_; i++) {
        values_[i] -= rhs[i];
    }
}

void Vector::subtract_vector_to(const double* rhs) {
    for(int i = 0; i< size_; i++) {
        values_[i] = rhs[i] - values_[i];
    }
}



/**
 * subtract  a subvector with length @subvec_size from the class member @_vector
 * from the location @iloc
 *
 * @param iloc the starting location to subtract the subvector
 * @param subvec_size the size of the subvector
 */
void Vector::subtract_subvector(int iloc, int subvec_size, const double* subvector) {
    for (int i = 0; i < subvec_size; i++) {
        values_[i + iloc - 1] -= subvector[i];
    }
}

/*copy all the entries from another vector*/
void Vector::copy_vector(const double* rhs) {
    for (int i = 0; i < size_; i++) {
        values_[i] = (double) rhs[i];
    }
}

/*copy all the entries from another vector*/
void Vector::copy_vector(std::shared_ptr<const Vector> rhs) {
    assert(size_ <= rhs->Dim());
    for (int i = 0; i < size_; i++) {
        values_[i] = rhs->values()[i];
    }
}


/** copy a subvector from member_vector from (Location) to (Location+subvector_size) to the pointer (results)*/
void Vector::get_subVector(int Location, int subvector_size,
                           std::shared_ptr<Vector> rhs) const {
    //TODO::test it! not sure if it will work...
    double* tmp = &(values_[Location - 1]);
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

void Vector::write_to_file(const char* name,
                           Ipopt::SmartPtr<Ipopt::Journalist> jnlst,
                           Ipopt::EJournalLevel level,
                           Ipopt::EJournalCategory category,
                           Solver qpsolver) {
#if DEBUG
#if PRINT_QP_IN_CPP
    const char * var_type;
    var_type = (qpsolver==QPOASES)?"real_t":"double const";

    jnlst->Printf(level, category, "%s %s[%i] = {", var_type,name,Dim());
    for (int i = 0; i < Dim(); i++) {
        if (i % 10 == 0 && i > 1)
            jnlst->Printf(level, category, "\n");
        if (i == Dim() - 1)
            jnlst->Printf(level, category, "%23.16e};\n\n", values_[i]);
        else
            jnlst->Printf(level, category, "%23.16e, ", values_[i]);
    }
#else
    //print in file
    for (int i = 0; i < Dim(); i++) {
        jnlst->Printf(level, category, "%23.16e\n", values_[i]);
    }
#endif
#endif
}
//
//Vector Vector::operator+(const Vector& rhs) {
//    assert(this->size_==rhs.size_);
//    auto added_result = Vector(this->size_);
//    for (int i = 0; i < this->size_; i++)
//        added_result.setValueAt(i, this->values()[i] + rhs.values()[i]);
//    return added_result;
//}
//
//const Vector Vector::operator-(const Vector &rhs) {
//    assert(this->size_==rhs.size_);
//    auto added_result = Vector(this->size_);
//    for (int i = 0; i < this->size_; i++)
//        added_result.setValueAt(i, this->values()[i] - rhs.values()[i]);
//    return added_result;
//}


double Vector::times(std::shared_ptr<Vector> rhs) {
    double product = 0;
    for (int i = 0; i < size_; i++) {
        product += values_[i] * rhs->values()[i];
    }
    return product;
}

double Vector::getOneNorm() const {
    double oneNorm = 0;
    for(int i=0; i<size_; i++) {
        oneNorm += fabs(values_[i]);
    }
    return oneNorm;
}


}

