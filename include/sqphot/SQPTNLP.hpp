/* Copyright (C) 2019
 * All Rights Reserved.
 *
 * Authors: Xinyi Luo
 * Date:2019-06
 */
#ifndef SQPHOTSTART_SQPTNLP_HPP
#define SQPHOTSTART_SQPTNLP_HPP

#include <memory>
#include "sqphot/Utils.hpp"
#include "IpTNLP.hpp"
#include "sqphot/Vector.hpp"
#include "sqphot/SpTripletMat.hpp"

namespace SQPhotstart {
/**
 * This is part of SQPhotstart
 *
 * This class enables user to read data from NLP class object with more friendly
 * names and the use of Matrix and Vector objects for data.
 *
 */
class SQPTNLP {

public:

    /** @brief constructor that copies nlp to _nlp as a local data reader*/
    SQPTNLP(Ipopt::SmartPtr<Ipopt::TNLP> nlp);

    /** Default destructor*/
    virtual ~SQPTNLP();

    /**
     *@brief get the bounds information from the NLP object
     */
    virtual bool Get_bounds_info(std::shared_ptr<Vector> x_l,
                                 std::shared_ptr<Vector> x_u,
                                 std::shared_ptr<Vector> c_l,
                                 std::shared_ptr<Vector> c_u);

    /*
     * @brief Get the starting point from the NLP object.
     * TODO: add options_ to enable user to choose if to use default input or not
     */
    virtual bool
    Get_starting_point(std::shared_ptr<Vector> x_0, std::shared_ptr<Vector> lambda_0);

    /**
     *@brief Evaluate the objective value
     */
    virtual bool Eval_f(std::shared_ptr<const Vector> x, double& obj_value);

    /**
     * @brief Evaluate the constraints at point x
     *
     */
    virtual bool
    Eval_constraints(std::shared_ptr<const Vector> x, std::shared_ptr<Vector> constraints);

    /**
     *@brief Evaluate gradient at point x
     */
    virtual bool
    Eval_gradient(std::shared_ptr<const Vector> x, std::shared_ptr<Vector> gradient);

    /**
     * @brief Get the matrix structure of the Jacobian
     * Always call this before the first time using @Eval_Jacobian
     */
    virtual bool Get_Strucutre_Jacobian(std::shared_ptr<const Vector> x,
                                        std::shared_ptr<SpTripletMat> Jacobian);

    /**
     *@brief Evaluate Jacobian at point x
     */
    virtual bool
    Eval_Jacobian(std::shared_ptr<const Vector> x, std::shared_ptr<SpTripletMat> Jacobian);


    /**
     * @brief Get the structure of the Hessian
     * Always call this before the first time using @Eval_Hessian
     */
    virtual bool
    Get_Structure_Hessian(std::shared_ptr<const Vector> x, std::shared_ptr<const Vector> lambda,
                          std::shared_ptr<SpTripletMat> Hessian);

    /**
     *@brief Evaluate Hessian of Lagragian function at  (x, lambda)
     */
    virtual bool
    Eval_Hessian(std::shared_ptr<const Vector> x, std::shared_ptr<const Vector> lambda,
                 std::shared_ptr<SpTripletMat> Hessian);

    /**
     * @brief This function shifts the initial starting point to be feasible to the bound constraints
     * @param x initial starting point
     * @param x_l lower bound constraints
     * @param x_u upper bound constraints
     */
    virtual bool
    shift_starting_point(std::shared_ptr<Vector> x, std::shared_ptr<const Vector> x_l,
                         std::shared_ptr<const Vector> x_u);

public:
    NLPInfo nlp_info_; /**< the struct record the number of variables, number of
                               constraints, number of nonzeoro entry of Hessian and that of Jacobian
                               Please check Types.hpp for details*/
    Ipopt::SmartPtr<Ipopt::TNLP> nlp_;/**< a local nlp reader */

private:
    /** Default constructor*/
    SQPTNLP();


    /** Copy Constructor */
    SQPTNLP(const SQPTNLP&);

    /** Overloaded Equals Operator */
    void operator=(const SQPTNLP&);
    //@}
};
}

#endif //CPP_SQPTNLP_HPP
