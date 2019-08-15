#include <sqphot/QOREInterface.hpp>

namespace SQPhotstart {

QOREInterface::QOREInterface(Index_info nlp_info)
{
	qp_=0;
      int rv =QPNew(&qp_, nlp_info.nVar, nlp_info.nCon, nlp_info.nnz_jac_g, nlp_info.nnz_h_lag);
//      std::cout <<" QOREInterface is created"<<std::endl;

}
}
