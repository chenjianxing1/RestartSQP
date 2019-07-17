/* Copyright (C) 2019
 * All Rights Reserved.
 *
 * Authors: Xinyi Luo
 * Date:2019-06
 */
#ifndef SQPHOTSTART_SQPTNLP_HPP
#define SQPHOTSTART_SQPTNLP_HPP

#include <memory>
#include <sqphot/Utils.hpp>
#include <sqphot/MyNLP.hpp>
#include <coin/IpTNLP.hpp>
#include <sqphot/Vector.hpp>
#include <sqphot/Matrix.hpp>

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
        SQPTNLP(SmartPtr<Ipopt::TNLP> nlp);
        
        /** Default destructor*/
        virtual ~SQPTNLP();
        
        /**
         *@brief get the bounds information from the NLP object
         */
        virtual bool Get_bounds_info(shared_ptr<Vector> x_l,
                                     shared_ptr<Vector> x_u,
                                     shared_ptr<Vector> c_l,
                                     shared_ptr<Vector> c_u);
        
        /*
         * @brief Get the starting point from the NLP object.
         * TODO: add options to enable user to choose if to use default input or not
         */
        virtual bool Get_starting_point(shared_ptr<Vector> x_0, shared_ptr<Vector> lambda_0);
        
        /**
         *@brief Evaluate the objective value
         */
        virtual bool Eval_f(shared_ptr<const Vector> x, Number &obj_value);
        
        /**
         * @brief Evaluate the constraints at point x
         *
         */
        virtual bool Eval_constraints(shared_ptr<const Vector> x, shared_ptr<Vector> constraints);
        
        /**
         *@brief Evaluate gradient at point x
         */
        virtual bool Eval_gradient(shared_ptr<const Vector> x, shared_ptr<Vector> gradient);
        
        /**
         * @brief Get the matrix structure of the Jacobian
         * Always call this before the first time using @Eval_Jacobian
         */
        virtual bool Get_Strucutre_Jacobian(shared_ptr<const Vector> x, shared_ptr<SpTripletMat> Jacobian);
        
        /**
         *@brief Evaluate Jacobian at point x
         */
        virtual bool Eval_Jacobian(shared_ptr<const Vector> x, shared_ptr<SpTripletMat> Jacobian);
        
        
        /**
         * @brief
         * @param x the
         * @param Jacobian
         * @return
         */
        virtual bool Get_Structure_Hessian(shared_ptr<const Vector> x, shared_ptr<const Vector> lambda,
                                           shared_ptr<SpTripletMat> Hessian);
        
        /**
         *@brief Evaluate Hessian of Lagragian function at  (x, lambda)
         */
        virtual bool
        Eval_Hessian(shared_ptr<const Vector> x, shared_ptr<const Vector> lambda, shared_ptr<SpTripletMat> Hessian);
        
        /**
         * @brief
         * @param x
         * @param x_l
         * @param x_u
         * @return
         */
        virtual bool
        shift_starting_point(shared_ptr<Vector> x, shared_ptr<const Vector> x_l, shared_ptr<const Vector> x_u);
        
    public:
        Index_info nlp_info_; /**< the struct record the number of variables, number of
                               constraints, number of nonzeoro entry of Hessian and that of Jacobian
                               Please check Types.hpp for details*/
        SmartPtr<Ipopt::TNLP> nlp_;/**< a local nlp reader */
        
    private:
        /** Default constructor*/
        SQPTNLP();
        
        
        /** Copy Constructor */
        SQPTNLP(const SQPTNLP &);
        
        /** Overloaded Equals Operator */
        void operator=(const SQPTNLP &);
        //@}
    };
}

#endif //CPP_SQPTNLP_HPP
