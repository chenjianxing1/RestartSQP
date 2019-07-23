/* Copyright (C) 2019
* All Rights Reserved.
*
* Authors: Xinyi Luo
* Date:2019-07-23
*/
#include <sqphot/Utils.hpp>
#include <sqphot/Types.hpp>
#include <coin/IpTNLP.hpp>
#include "Vector.hpp"

#ifndef _SQPHOTSTART_OPTTEST_HPP
#define _SQPHOTSTART_OPTTEST_HPP

namespace SQPhotstart {

    enum TestOption {
        NO_TEST = -99,
        TEST_1ST_ORDER = 1,
        TEST_2ND_ORDER = 2,
        TEST_STATIONARITY = -1,
        TEST_COMPLEMENTARITY = -2,
        TEST_FEASIBILITY = -3,
        TEST_DUAL_FEASIBILITY = -4,
        TEST_ALL = 0
    };

    /**
     * @brief This is a pure virtual base class for Optimality Test
     * */
    class OptTest {

    public:
        /** Default constructor*/
        OptTest() = default;

        /**@brief Default Destructor*/
        virtual ~OptTest() = default;

        /**
         * @brief Test the KKT conditions
         */
        virtual bool Check_KKTConditions() { return false; }

        /**
         * @brief Test the Second-order optimality conditions
         */
        virtual bool Check_SecondOrder() { return false; }

        /**
         * @brief Check the Feasibility conditions;
         */
        virtual bool Check_Feasibility() { return false; }

        /**
         * @brief Check the sign of the multipliers
         */
        virtual bool Check_Dual_Feasibility() { return false; }

        /**
         * @brief Check the complementarity condition
         *
         */
        virtual bool Check_Complementarity() { return false; }

        /**
         * @brief Check the Stationarity condition
         */
        virtual bool Check_Stationarity() { return false; }

        /** Public class members*/
    public:
        bool primal_feasibility_ = false;
        bool dual_feasibility_ = false;
        bool complementarity_ = false;
        bool Stationarity_ = false;
        bool Second_order_opt_ = false;

    private:

        //@{
        /** Copy Constructor */
        OptTest(const OptTest&);

        /** Overloaded Equals Operator */
        void operator=(const OptTest&);
        //@}

    };

    /**
     * @brief Derived class for OptTest. Test the optimality condition for a given NLP
     * problem
     */
    class NLP_OptTest : public OptTest {
    public:

        NLP_OptTest(Ipopt::SmartPtr<Ipopt::TNLP> nlp, shared_ptr<Vector> x_k);

        NLP_OptTest() = default;

        virtual ~NLP_OptTest() = default;

        /**
         * @brief Test the KKT conditions
         */
        bool Check_KKTConditions() override;

        /**
         * @brief Test the Second-order optimality conditions
         */
        bool Check_SecondOrder() override;

        /**
         * @brief Check the Feasibility conditions;
         */
        bool Check_Feasibility() override;

        /**
         * @brief Check the sign of the multipliers
         */
        bool Check_Dual_Feasibility() override;

        /**
         * @brief Check the complementarity condition
         *
         */
        bool Check_Complementarity() override;

        /**
         * @brief Check the Stationarity condition
         */
        bool Check_Stationarity() override;


        /** Public class members*/
    public:
        bool primal_feasibility_ = false;
        bool dual_feasibility_ = false;
        bool complementarity_ = false;
        bool stationarity_ = false;
        bool second_order_opt_ = false;

    private:
        //@{
        /** Copy Constructor */
        NLP_OptTest(const NLP_OptTest&);

        /** Overloaded Equals Operator */
        void operator=(const NLP_OptTest&);
        //@}
    };


    /**
     * @brief Derived class for OptTest. Test the optimality condition for a given QP
     * problem
     */
    class QP_OptTest : public OptTest {
    public:
        QP_OptTest() = default;

        virtual ~QP_OptTest() = default;

        /**
         * @brief Test the KKT conditions
         */
        bool Check_KKTConditions();

        /**
         * @brief Test the Second-order optimality conditions
         */
        bool Check_SecondOrder();

        /**
         * @brief Check the Feasibility conditions;
         */
        bool Check_Feasibility();

        /**
         * @brief Check the sign of the multipliers
         */
        bool Check_Dual_Feasibility();

        /**
         * @brief Check the complementarity condition
         *
         */
        bool Check_Complementarity();

        /**
         * @brief Check the Stationarity condition
         */
        bool Check_Stationarity();


        /** Public class members*/
    public:
        bool primal_feasibility_ = false;
        bool dual_feasibility_ = false;
        bool complementarity_ = false;
        bool stationarity_ = false;
        bool second_order_opt = false;

    private:
        //@{
        /** Copy Constructor */
        QP_OptTest(const QP_OptTest&);

        /** Overloaded Equals Operator */
        void operator=(const QP_OptTest&);
        //@}
    };
}

#endif //SQPHOTSTART_OPTTEST_HPP

