[![Build Status](https://travis-ci.org/lanl-ansi/SQPhotstart.svg?branch=master)](https://travis-ci.org/lanl-ansi/SQPhotstart/)
[![Code Coverage](https://codecov.io/gh/lanl-ansi/SQPhotstart/branch/master/graphs/badge.svg)](https://codecov.io/gh/lanl-ansi/SQPhotstart)
<p align="center">
<img src="media/fire.png" width="150">
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


###### Logo credits: <div>Icons made by <a href="https://www.freepik.com/?__hstc=57440181.4a6c3f6caafe83440921e359e8f146ff.1562721776119.1562764946690.1562881160669.3&__hssc=57440181.4.1562881160669&__hsfp=3697938389" title="Freepik">Freepik</a> from <a href="https://www.flaticon.com/"                 title="Flaticon">www.flaticon.com</a> is licensed by <a href="http://creativecommons.org/licenses/by/3.0/"                 title="Creative Commons BY 3.0" target="_blank">CC 3.0 BY</a></div>
