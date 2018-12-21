# acopf
acopf.m is a Matlab function that finds the AC optimal power flow of problems that have the Matpower case format using the IPOPT solver. The requierements for runing this function is to have installed IPOPT package in Matlab. Matpower package is not required.
For more information about IPOPT, visit: https://projects.coin-or.org/Ipopt/wiki/MatlabInterface

acopf.m needs as an input the name of the case file. This case file has to have the Matpower case format, visit http://www.pserc.cornell.edu/matpower/ to download Matpower and their case examples.

acopf.m creates a file (.out) containing the detailed solution founded by the IPOPT solver, as well as a Matlab file (.mat) containing the variables of the optimal solution
