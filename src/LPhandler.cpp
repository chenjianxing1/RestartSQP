/* Copyright (C) 2019
 * All Rights Reserved.
 *
 * Authors: Xinyi Luo
 * Date:2019-07
 */

#include <sqphot/LPhandler.hpp>

namespace SQPhotstart {
/**
 * Default constructor
 */

LPhandler::LPhandler() {};

/**
 *Default destructor
 */
LPhandler::~LPhandler() {
}


/**
 * Get the optimal solution from the LPhandler_interface
 *
 *This is only an interface for user to avoid call interface directly.
 * @param p_k       the pointer to an empty array with the length equal to the size of the QP subproblem
 */
void LPhandler::GetOptimalSolution(double* p_k) {
    solverInterface_->get_optimal_solution(p_k);

}


/**
 * @brief solve the LP problem with bounds and objectives defined by the class
 * members
 */

void LPhandler::solveLP(shared_ptr<SQPhotstart::Stats> stats,
                        shared_ptr<Options> options) {
    solverInterface_->optimizeLP(stats, options);
    lp_obj_ = solverInterface_->get_obj_value();

}


/**
 * @brief This function initializes all class members which
 * will be used in _qp_interface
 *
 * @param nlp_info
 *          it contains the information about number
 *          of variables, number of constraints, number of
 *          elements in the Hessian and that of
 *          Jacobian. The definition of Index_info is
 *          in Types.hpp
 * @param Variable_type
 *          it specifies how the variables are bounded.
 *          Please check @ClassifyConstraintType in Algorithm.hpp
 *          for more details
 * @param qptype
 *          it specifies which type of QP is going to
 *          be solved. It can be either LP, or QP, or SOC
 */
void LPhandler::allocate(SQPhotstart::Index_info nlp_info) {
    solverInterface_ = make_shared<qpOASESInterface>(nlp_info, LP);

}


} // namespace SQPhotstart



