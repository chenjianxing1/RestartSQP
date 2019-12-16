#include <AmplTNLP.hpp>
#include <IpIpoptApplication.hpp>
#include <IpTNLP.hpp>
#include <cstddef>
#include <iostream>
#include <memory>
#include <sqphot/QPhandler.hpp>
#include <sqphot/SqpAlgorithm.hpp>
#include <sqphot/SqpTNlp.hpp>
#include <sqphot/Utils.hpp>
#include <stdio.h>
#include <string.h>
#include <string>
#include <sys/stat.h>

using namespace Ipopt;
using namespace SQPhotstart;
using namespace std;

inline bool exist(const std::string& name)
{
  struct stat buffer;
  return (stat(name.c_str(), &buffer) == 0);
}

class Table_Writer
{

public:
  Table_Writer(const std::string& filename)
  {

    if (!exist(filename)) {
      file = fopen(filename.c_str(), "a"); // append at the end
      fprintf(file, "%10s   %10s    %10s    %10s    %10s    %10s    ", "name",
              "nVar", "nConstr", "iter", "QP_iter", "exitflag");

      fprintf(file, "%23s    %23s    %23s    %23s    %23s    %23s\n",
              "objective", "||p||", "primal_violation", "dual_violation",
              "stationarity_violation", "compl_violation");

    } else {
      file = fopen(filename.c_str(), "a"); // append at the end
    }
  }

  ~Table_Writer()
  {
    fclose(file);
  }

  /**
   * @brief Write a brief summary for each problem being solved,
   */
  void write_in_brief(const std::string& pname, /**<name of the problem*/
                      SqpAlgorithm& alg)
  {

    std::size_t found = pname.find_last_of("/\\");
    OptimalityStatus opt_status = alg.get_opt_status();
    shared_ptr<Stats> stats;
    stats = alg.get_stats();

    fprintf(file, "%10s   %10d    %10d    %10d    %10d    %10d    ",
            pname.substr(found + 1).c_str(), alg.get_num_var(),
            alg.get_num_constr(), stats->iter, stats->qp_iter,
            alg.get_exit_flag());
    fprintf(file,
            "%23.16e    %23.16e    %23.16e    %23.16e    %23.16e    %23.16e\n",
            alg.get_final_objective(), alg.get_norm_p(),
            opt_status.primal_violation, opt_status.dual_violation,
            opt_status.stationarity_violation, opt_status.compl_violation);
  }

private:
  FILE* file;
};

int main(int argc, char** args)
{
  // Create the SQP algorithm object
  SqpAlgorithm alg;

  // Create an Ipopt::AmplTNLP object to handle the AMPL model
  SmartPtr<OptionsList> dummy_options = new OptionsList();
  SmartPtr<TNLP> ampl_tnlp =
      new AmplTNLP(ConstPtr(alg.get_jnlst()), dummy_options, args);
  shared_ptr<SqpTNlp> sqp_nlp = make_shared<SqpTNlp>(ampl_tnlp);

  // Solve the AMPL model
  alg.initialize(sqp_nlp, args[1]);
  alg.Optimize();

  shared_ptr<Table_Writer> writer = make_shared<Table_Writer>("result_table");
  writer->write_in_brief(args[1], alg);

  return 0;
}
