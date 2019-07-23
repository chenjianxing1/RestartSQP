//
// Created by Xinyi Luo on 7/23/19.
//
#include <sqphot/Utils.hpp>

#ifndef _SQPHOTSTART_OPTTEST_HPP
#define _SQPHOTSTART_OPTTEST_HPP

namespace SQPhotstart {
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


    public:
        bool primal_feasibility_;
        bool dual_feasibility_;
        bool complementarity_;
        bool Stationarity_;
        bool Second_order_opt_;

    private:

        //@{
        /** Copy Constructor */
        OptTest(const OptTest&);

        /** Overloaded Equals Operator */
        void operator=(const OptTest&);
        //@}

    };

    class NlP_OptTest : public OptTest {
    public:
        NlP_OptTest() = default;

        virtual ~NlP_OptTest() = default;

        /**
         * @brief Test the KKT conditions
         */
        bool Check_KKTConditions() { return false; }

        /**
         * @brief Test the Second-order optimality conditions
         */
        bool Check_SecondOrder() { return false; }

        /**
         * @brief Check the Feasibility conditions;
         */
         bool Check_Feasibility() { return false; }

        /**
         * @brief Check the sign of the multipliers
         */
         bool Check_Dual_Feasibility() { return false; }

        /**
         * @brief Check the complementarity condition
         *
         */
         bool Check_Complementarity() { return false; }

        /**
         * @brief Check the Stationarity condition
         */
         bool Check_Stationarity() { return false; }


    public:
        bool primal_feasibility_;
        bool dual_feasibility_;
        bool complementarity_;
        bool Stationarity_;
        bool Second_order_opt_;

    private:
        //@{
        /** Copy Constructor */
        NlP_OptTest(const NLP_OptTest&);

        /** Overloaded Equals Operator */
        void operator=(const NlP_OptTest&);
        //@}
    };


    class QP_OptTest : public OptTest {
    public:
        QP_OptTest() = default;

        virtual ~QP_OptTest() = default;

        /**
         * @brief Test the KKT conditions
         */
        bool Check_KKTConditions() { return false; }

        /**
         * @brief Test the Second-order optimality conditions
         */
        bool Check_SecondOrder() { return false; }

        /**
         * @brief Check the Feasibility conditions;
         */
        bool Check_Feasibility() { return false; }

        /**
         * @brief Check the sign of the multipliers
         */
        bool Check_Dual_Feasibility() { return false; }

        /**
         * @brief Check the complementarity condition
         *
         */
        bool Check_Complementarity() { return false; }

        /**
         * @brief Check the Stationarity condition
         */
        bool Check_Stationarity() { return false; }


    public:
        bool primal_feasibility_;
        bool dual_feasibility_;
        bool complementarity_;
        bool Stationarity_;
        bool Second_order_opt_;

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

