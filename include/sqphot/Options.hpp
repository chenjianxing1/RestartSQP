/**All Rights Reserved.
*
* Authors: Xinyi Luo
* Date:2019-05
*/
#ifndef SQPHOTSTART_OPTIONS_HPP
#define SQPHOTSTART_OPTIONS_HPP

#include <string>
#include <sqphot/Types.hpp>

namespace SQPhotstart {
class Options {
public:
    /**Default Constructor */
    Options();

    /**Destructor**/
    ~Options();

    /** Sets all options to the default values */
    int setToDefault();

public:
    /** Public Member Variables*/

    bool EnablePertubation;
    bool second_order_correction;
    double tol;
    int iter_max;
    int printLevel;


    /**solver choice*/
    //@{
    QPSolver QPsolverChoice;
    LPSolver LPsolverChoice;

    //@}

    /** trust-region update parameters*/
    //@{
    double delta;
    double delta_max;
    double delta_min;
    double eta_c;
    double eta_e;
    double eta_s;
    double gamma_c;
    double gamma_e;
    //@}

    /** optimality test parameters */
    //@{
    bool auto_gen_tol = false;
    double active_set_tol;
    double opt_compl_tol;
    double opt_dual_fea_tol;
    double opt_prim_fea_tol;
    double opt_second_tol;
    double opt_tol;
    //@}

    /**QPsolver options */
    //@{

    int lp_maxiter;
    int qpPrintLevel;
    int qp_maxiter;
    //@}

    /** penalty update parameters*/
    //@{

    //double decrease_parm;
    //double decrease_tol;
    //int update_method

    bool penalty_update;
    double eps1;//FIXME: it may need to be changed inside the main loop
    double eps1_change_parm;
    double eps2;
    double increase_parm;
    double penalty_update_tol;
    double rho;
    double rho_max;
    int penalty_iter_max;
    //@}
};//ENDCLASS
}
#endif /* SQPHOTSTART_OPTIONS_HPP */

