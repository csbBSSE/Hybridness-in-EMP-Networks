# Hybridness-in-EMT-networks
This repository contains a MATLAB Code, a python Code, and 6 R codes.


MATLAB code description
Matlab codes takes input as continous-time simulations of any EMP network. The struucture of input data file is like having as many columns as the number of nodes in the network and as many rows as the number of parameter sets used to simulate the model. Each row represents a stable steady state. The matlab code first normalizes each data column using Z- Score standardization, discretizes the data columns, and finally calculates the frequency of each steady state. It also comapres the frequency distribution of each network-perturbed data file with the corresponding wild-type file. The output of the code is a data file that contains steady states of the perturbed file and that of wild type and and their respective frequencies.     

R Codes have their descriptions in the comments.


Python Code:
Python codes takes network topology files as input and calculates the number of psoitive/negative feedback loops in the network.

