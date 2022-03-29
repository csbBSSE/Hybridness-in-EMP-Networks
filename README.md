# Hybridness-in-EMT-networks
This repository contains a MATLAB Code, a python Code, and few R codes.


MATLAB code description
Matlab codes takes input as continous-time simulations of any EMP network. The struucture of input data file is like having as many columns as the number of nodes in the network and as many rows as the number of parameter sets used to simulate the model. Each row represents a stable steady state. The matlab code first normalizes each data column using Z- Score standardization, discretizes the data columns, and finally calculates the frequency of each steady state. It also comapres the frequency distribution of each network-perturbed data file with the corresponding wild-type file. The output of the code is a data file that contains steady states of the perturbed file and that of wild type and and their respective frequencies.     

R Code Description

R Code #1: JSD
Thsi code takes input as steady state frequency file of a perturbed network and compares it with the steady state frequency file of the corresponding wild type file to calculate JSD using the JSD function.

R Code #2: J metric
This code takes continuous-time simulation file of a network as input and it gives J metric as output. It finds J metric in the following steps: Firstly, it finds the correlation (Pearson) matrix of the solution matrix. Secondly, it calculates the upper traingular matrix of the correlation matrix. Finally, it sums up the elements of the upper triangular matrix (without diagonal elements) by taking the absolute values of the elements.
