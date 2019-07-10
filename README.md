DEPENDENCIES
-------
cmake: http://www.cmake.org (Version 3.2 or better)

Ipopt: https://projects.coin-or.org/Ipopt

Included ThirdParty Libraries
-------

qpOASES: https://projects.coin-or.org/qpOASES/wiki/WikiStart

OPTIONAL
-------

Eigen: https://github.com/eigenteam/eigen-git-mirror

Gurobi: http://www.gurobi.com

Cplex: https://www-01.ibm.com/software/commerce/optimization/cplex-optimizer/

-------

Follow these simple instructions:
1) `cd SQPhotstart`

2) If your solver is not installed in a default location, you need to add the following to your `~/.bash_profile` (edit or create the file named `bash_profile` under your home and add the following lines):
* `export SOLVERNAME_ROOT_DIR="your_location"`, where `your_location` contains the corresponding include/headers and lib/libraries.

For instance: 

* `export IPOPT_ROOT_DIR="/Users/yourname/Dev/CoinIpopt/build"`


Note: If Cplex and Gurobi were installed using the standard installer, Cmake should automatically locate the latest version on your system, no need to add anything to your bash_profile.


3) You're now ready to compile everything, just enter:

* `mkdir build`

* `cd build`

* `cmake ..`

All dependencies are switched off by default, to enable a solver that is installed on your system,  
append the flag `-D$Solvername$=ON`, e.g., `cmake -DIpopt=ON -DGurobi=ON -DCplex=ON ..`.

Note: To build an Xcode project append `-G Xcode` to the command above

* `make -j4`

This will build the solver library and the tests under `SQPhotstart/test`.

The corresponding binaries will then appear under `SQPhotstart/bin/Release`.

