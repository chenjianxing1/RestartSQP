<p align="center">
<img src="media/fire.svg" width="150">
<H2 align="center"> SQPhotstart: A Sequential Quadratic Programming Solver </H2>
</p>
<p align="center"> for Constrained Nonlinear Optimization </p>

DEPENDENCIES
-------
cmake: http://www.cmake.org (Version 3.2 or better)

Ipopt: https://projects.coin-or.org/Ipopt

Included ThirdParty Libraries
-------

qpOASES: https://projects.coin-or.org/qpOASES/wiki/WikiStart

The qpOASES library included is not compiled with HSL http://www.hsl.rl.ac.uk/catalogue/ linear solvers.
If you compile qpOASES with HSL solvers, you just need to replace the libraries libqpOASES.dylib (for MacOS) and libqpOASES.so (for Linux) found under `ThirdParty/qpOASES-3.2.1/bin`.

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

All dependencies are switched off by default, except for Ipopt, to enable a solver that is installed on your system,  
append the flag `-D$Solvername$=ON`, e.g., `cmake -DGurobi=ON -DCplex=ON ..`.

Note: To build an Xcode project append `-G Xcode` to the command above

* `make -j4`

This will build the solver library and the tests under `SQPhotstart/test`.

The corresponding binaries will then appear under `SQPhotstart/bin/Release`.

