/* Copyright (C) 2019
* All Rights Reserved.
*
* Authors: Xinyi Luo
* Date:2019-07
*/
#include <sqphot/SQPTNLP.hpp>


namespace SQPhotstart {
    /** Default constructor*/
    SQPTNLP::SQPTNLP(SmartPtr<TNLP> nlp) {
        nlp_ = nlp;
        Ipopt::TNLP::IndexStyleEnum index_style;
        nlp_->get_nlp_info(nlp_info_.nVar, nlp_info_.nCon, nlp_info_.nnz_jac_g,
                           nlp_info_.nnz_h_lag, index_style);
        //assert(index_style == Ipopt::TNLP::C_STYLE);
    }


    /** Default constructor*/
    SQPTNLP::~SQPTNLP() {
    }

    /**
     *@name get the bounds information from the NLP object
     */
    bool SQPTNLP::Get_bounds_info(shared_ptr<Vector> x_l, shared_ptr<Vector> x_u,
                                  shared_ptr<Vector> c_l, shared_ptr<Vector> c_u) {

        nlp_->get_bounds_info(nlp_info_.nVar, x_l->values(), x_u->values(),
                              nlp_info_.nCon, c_l->values(), c_u->values());

        return true;
    }

    /*
     * @name Get the starting point from the NLP object.
     * TODO: add options to enable user to choose if to use default input or not
     */
    bool
    SQPTNLP::Get_starting_point(shared_ptr<Vector> x_0, shared_ptr<Vector> lambda_0) {
        nlp_->get_starting_point(nlp_info_.nVar, true, x_0->values(),
                                 false, NULL, NULL, nlp_info_.nCon, true,
                                 lambda_0->values());

        return true;
    }

    /**
     *@name Evaluate the objective value
     */
    bool SQPTNLP::Eval_f(shared_ptr<const Vector> x, Number& obj_value) {
        nlp_->eval_f(nlp_info_.nVar, x->values(), true, obj_value);
        return true;
    }

    /**
     * @name Evaluate the constraints at point x
     *
     */
    bool SQPTNLP::Eval_constraints(shared_ptr<const Vector> x,
                                   shared_ptr<Vector> constraints) {
        nlp_->eval_g(nlp_info_.nVar, x->values(), true, nlp_info_.nCon,
                     constraints->values());
        return true;
    }

    /**
     *@name Evaluate gradient at point x
     */
    bool SQPTNLP::Eval_gradient(shared_ptr<const Vector> x, shared_ptr<Vector> gradient) {
        nlp_->eval_grad_f(nlp_info_.nVar, x->values(), true, gradient->values());
        return true;
    }

    /**
     *@name Evaluate Jacobian at point x
     */

    bool SQPTNLP::Get_Strucutre_Jacobian(shared_ptr<const Vector> x,
                                         shared_ptr<SpTripletMat> Jacobian) {
        nlp_->eval_jac_g(nlp_info_.nVar, x->values(), true, nlp_info_.nCon,
                         nlp_info_.nnz_jac_g,
                         Jacobian->RowIndex(), Jacobian->ColIndex(), NULL);
        return true;
    }


    /**
     *@name
     *@param x
     *@param Jacobian
     */
    bool SQPTNLP::Eval_Jacobian(shared_ptr<const Vector> x,
                                shared_ptr<SpTripletMat> Jacobian) {
        nlp_->eval_jac_g(nlp_info_.nVar, x->values(), true, nlp_info_.nCon,
                         nlp_info_.nnz_jac_g,
                         Jacobian->RowIndex(), Jacobian->ColIndex(), Jacobian->MatVal());
        return true;
    }

    /**
     *
     * @param x
     * @param lambda
       * @param Hessian
     * @return
     */
    bool SQPTNLP::Get_Structure_Hessian(shared_ptr<const Vector> x,
                                        shared_ptr<const Vector> lambda,
                                        shared_ptr<SpTripletMat> Hessian) {
        nlp_->eval_h(nlp_info_.nVar, x->values(), true, 1.0, nlp_info_.nVar,
                     lambda->values(), true,
                     nlp_info_.nnz_h_lag, Hessian->RowIndex(), Hessian->ColIndex(), NULL);

        return true;
    }


    /**
     *@name Evaluate Hessian of Lagragian function at  (x, lambda)
     */
    bool
    SQPTNLP::Eval_Hessian(shared_ptr<const Vector> x, shared_ptr<const Vector> lambda,
                          shared_ptr<SpTripletMat> Hessian) {
        nlp_->eval_h(nlp_info_.nVar, x->values(), true, 1, nlp_info_.nVar,
                lambda->values(), true,
                     nlp_info_.nnz_h_lag, Hessian->RowIndex(), Hessian->ColIndex(), Hessian->MatVal());

        return true;
    }

    /**
     * @name This function shifts the initial starting point to be feasible to the bound constraints
     * @param x initial starting point
     * @param x_l lower bound constraints
     * @param x_u upper bound constraints
     */
    bool SQPTNLP::shift_starting_point(shared_ptr<Vector> x, shared_ptr<const Vector> x_l,
                                       shared_ptr<const Vector> x_u) {
        for (int i = 0; i < x->Dim(); i++) {
            assert(x_l->values()[i] <= x_u->values()[i]);
            if (x_l->values()[i] > x->values()[i]) {
                x->setValueAt(i, x_l->values()[i]);
            } else if (x->values()[i] > x_u->values()[i]) {
                x->setValueAt(i, x_u->values()[i]);
            }
        }
        return true;
    }

}
