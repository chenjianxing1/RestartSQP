/* Copyright (C) 2019
* All Rights Reserved.
*
* Authors: Xinyi Luo
* Date:2019-07-23
*/

#include <sqphot/OptTest.hpp>
#include <coin/IpTNLP.hpp>
#include <sqphot/Options.hpp>
#include <sqphot/Vector.hpp>

namespace SQPhotstart {

    NLP_OptTest::NLP_OptTest(shared_ptr<SQPhotstart::Options> options,
                             Index_info nlp_index_info) {
        opt_tol_ = options->opt_tol;
        opt_compl_tol_ = options->opt_compl_tol;
        opt_dual_fea_tol_ = options->opt_dual_fea_tol;
        opt_prim_fea_tol_ = options->opt_prim_fea_tol;
        opt_second_tol_ = options->opt_second_tol;
        nVar_ = nlp_index_info.nVar;
        nCon_ = nlp_index_info.nCon;
    }

    bool NLP_OptTest::Check_KKTConditions() {
        Check_Feasibility();
        Check_Complementarity();
        Check_Dual_Feasibility();
        Check_Stationarity();
        return true;
    }

    bool NLP_OptTest::Check_Complementarity() {
        return true;
    }

    bool NLP_OptTest::Check_Feasibility() {

        return true;
    }


    bool NLP_OptTest::Check_Dual_Feasibility() {
        return true;
    }


    bool NLP_OptTest::Check_Dual_Feasibility(const double* multiplier,
                                             ConstraintType nlp_cons_type) {

        return true;
    }

    bool NLP_OptTest::Check_SecondOrder() {
        if (complementarity_ && dual_feasibility_ && stationarity_ &&
            primal_feasibility_) {
            //Check the 2nd order condition
        } else {
            //print error message
        }
        return true;
    }


    bool NLP_OptTest::Check_Stationarity() {
        return true;
    }

    bool NLP_OptTest::Check_Feasibility(double infea_measure) {
        if (primal_feasibility_ == false && infea_measure < opt_prim_fea_tol_)
            primal_feasibility_ = true;
    }


    bool QP_OptTest::Check_KKTConditions() {
        Check_Feasibility();
        Check_Complementarity();
        Check_Dual_Feasibility();
        Check_Stationarity();
        return true;
    }

    bool QP_OptTest::Check_Complementarity() {
        return true;
    }

    bool QP_OptTest::Check_Feasibility() {

        return true;
    }


    bool QP_OptTest::Check_Dual_Feasibility() {
        return true;
    }


    bool QP_OptTest::Check_SecondOrder() {
        if (complementarity_ && dual_feasibility_ && stationarity_ &&
            primal_feasibility_) {
            //Check the 2nd order condition
        } else {
            //print error message
        }
        return true;

    }


    bool QP_OptTest::Check_Stationarity() {
        return true;
    }


}
