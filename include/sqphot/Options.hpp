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

    double tol;
    int iter_max;
    int printLevel;
    bool second_order_correction;
    bool EnablePertubation;


    /**solver choice*/
    //@{
    std::string QPsolverChoice;
    std::string LPsolverChoice;

    //@}

    /** trust-region update parameters*/
    //@{
    double eta_c;
    double eta_s;
    double eta_e;
    double gamma_c;
    double gamma_e;
    double delta;
    double delta_max;
    //@}

    /** optimality test parameters */
    //@{
    TestOption testOption_NLP;
    bool auto_gen_tol = false;
    double active_set_tol;
    double opt_tol;
    double opt_compl_tol;
    double opt_dual_fea_tol;
    double opt_prim_fea_tol;
    double opt_second_tol;
    //@}

    /**QPsolver options */
    //@{
    TestOption testOption_QP;
    int qp_maxiter;
    int lp_maxiter;
    int qpPrintLevel;
    //@}

    /** penalty update parameters*/
    //@{
    double penalty_update_tol;
    bool penalty_update;
    double eps1;//FIXME: it may need to be changed inside the main loop
    double eps2;
    double rho;
    double increase_parm;
    int rho_max;
    int penalty_iter_max;

    //int update_method
    //double decrease_parm;
    //double decrease_tol;
    //@}
};//ENDCLASS
}
#endif /* SQPHOTSTART_OPTIONS_HPP */
