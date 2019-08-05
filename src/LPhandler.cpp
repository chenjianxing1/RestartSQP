/* Copyright (C) 2019
 * All Rights Reserved.
 *
 * Authors: Xinyi Luo
 * Date:2019-07
 */

#include <sqphot/LPhandler.hpp>

namespace SQPhotstart {

LPhandler::LPhandler(Index_info nlp_info):
    QPhandler(nlp_info),
    nlp_info_(nlp_info)
{

    solverInterface_ = make_shared<qpOASESInterface>(nlp_info,LP);
}
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
}




} // namespace SQPhotstart



