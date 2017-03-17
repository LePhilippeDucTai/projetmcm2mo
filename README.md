# projetmcm2mo

## Summary ##
This repository contains our Monte-Carlo project on calculating a Value-at-Risk and Conditional Value-at-Risk (expected shortfall) using an Importance Sampling method in order to reduce the variance on the VaR Monte-Carlo estimate. We rely our work on this article : Recursive Computation of Value-at-Risk and Conditional Value-at-Risk using MC and QMC, written by 
Olivier Bardou, Noufel Frikha, and Gilles Pagès.

## Requirements ##
Make sure to have installed a `g++` compiler (C++14) installed on your OS.
Make sure to install the armadillo library, found at : http://arma.sourceforge.net

## Compiling ##
Once this has been done, run on a terminal the following commands :
```
make clean
make
```
## Run ##
To run the program, type :
```
./mc
```
