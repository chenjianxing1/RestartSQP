#ifndef SQPHOTSTART_STATS_HPP
#define SQPHOTSTART_STATS_HPP
namespace SQPhotstart {
    class Stats {
    public:
        /* Default constructor*/
        Stats() {
            qp_iter = 0;
            iter = 0;
            qp_break_down = 0;
            penalty_change_trial = 0;
            penalty_change_Fail = 0;
            penalty_change_Succ = 0;
            soc_iter = 0;
        };

        /* Destructor*/

        ~Stats() {};

        /* add 1 to the value of class member iter*/
        inline void iter_addone() {
            iter++;
        }

        inline void iter_addValue(int n) {
            iter += n;
        };


        /* add 1 to the value of class member qp_iter*/
        inline void qp_iter_addone() {
            qp_iter++;
        };

        inline void qp_iter_addValue(int n) {
            qp_iter += n;
        };


        /* add 1 to the value of class member qp_break_down*/
        inline void qp_break_down_addone() {
            qp_break_down++;
        };

        /* add n to the value of class member qp_break_down*/
        inline void qp_break_down_addValue(int n) {
            qp_break_down += n;
        };

        /* add 1 to the value of class member penalty_change_Fail*/
        inline void penalty_change_Fail_addone() {
            penalty_change_Fail++;
        };

        /* add n to the value of class member qp_break_down*/
        inline void penalty_change_Fail_addValue(int n) {
            penalty_change_Fail += n;
        };


        /* add 1 to the value of class member penalty_change_trial*/
        void penalty_change_trial_addone() {
            penalty_change_trial++;
        };

        /* add n to the value of class member penalty_change_trial*/
        inline void penalty_change_trial_addValue(int n) {
            penalty_change_trial += n;
        };


        /* add 1 to the value of class member penalty_change_Succ*/
        inline void penalty_change_Succ_addone() {
            penalty_change_Succ++;
        };

        /* add n to the value of class member penalty_change_Succ*/
        inline void penalty_change_Succ_addValue(int n) {
            penalty_change_Succ += n;
        };


        /* add 1 to the value of class member soc_iter*/
        inline void soc_iter_addone() {
            soc_iter++;
        };

        /* add n to the value of class member soc_iter*/
        inline void soc_iter_addValue(int n) {
            soc_iter += n;
        };


        /* Member Variables */
    public:
        int qp_iter;
        int iter;
        int qp_break_down;
        int penalty_change_trial;
        int penalty_change_Fail;
        int penalty_change_Succ;
        int soc_iter;
    };

}//END_NAMESPACE_SQPHOTSTART

#endif /* SQPHOTSTART_STATS_HPP*/
