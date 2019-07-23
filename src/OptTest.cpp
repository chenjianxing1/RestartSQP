/* Copyright (C) 2019
* All Rights Reserved.
*
* Authors: Xinyi Luo
* Date:2019-07-23
*/

#include <sqphot/OptTest.hpp>
#include <coin/IpTNLP.hpp>
#include <sqphot/Vector.hpp>

namespace SQPhotstart {


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


    bool NLP_OptTest::Check_SecondOrder() {
        if (complementarity_ && dual_feasibility_ && stationarity_&&
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
