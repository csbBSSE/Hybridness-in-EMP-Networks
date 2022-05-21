
## Network Topology Metrics Explaining Enrichment of Hybrid Epithelial Mesenchymal Phenotypes in Metastasis

This resource provides the codes to reproduce all the key results described in [**Network Topology Metrics Explaining Enrichment of Hybrid Epithelial Mesenchymal Phenotypes in Metastasis**]([https://www.biorxiv.org/content/10.1101/2022.05.16.492000v1](https://www.biorxiv.org/content/10.1101/2022.05.16.492000v1))

The analysis done can be briefly described as:

* **Discrete Modelling**: Asynchonous Update on Ising model of 13 Wild-Type regulatory networks underlying Epithelial Mesenchymal Plasticity (EMP) and their single-edge perturbation counterparts, generated as described in the paper. This is done by the set of codes given [here]([https://github.com/csbBSSE/CSB-SCLC/tree/master/Additional_Codes/Fast-Bool](https://github.com/askhari139/Boolean.jl)).
* **Continuous Modelling**: **RACIPE** is used to generate and simulate an ensemble of ODE-Based models of the same set of networks. For this we have used **RACIPE-1.0** package which can be found [here](https://github.com/simonhb1990/RACIPE-1.0).
* **Network Metrics**: For the netowrks used in the study, custom codes are provided here. The code for hiLoop metrics is also provided here. 

### Topo files (raw data)

The topo files for **Edge Perturbation** can be found in the **Boolean** folder. The topo files for **Random Networks** can be found in the **RandomNets** folder. 

### Simulation Data
**Boolean Simulation data** for **Edge Perturbation** is provided in the **Boolean** folder and that for **Edge-weight analysis** is provided in **EdgeWeightPert** folder. Since **RACIPE** simulation data files are quite huge, they are uploaded to this [drive link](https://drive.google.com/drive/folders/1PKs5vHkXCoJm9Wcg7P4nBPdPrFJCxJ5B?usp=sharing).

### Additonal Codes
This folder contains all the codes required for data analysis the networks. The codes can be classified into two categories:

#### Simulation data analysis:
This involves discretizing the RACIPE outputs, frequency calculation for steady states, multistability calculation, and calculating frustration, hybridness and J_metric.
Discretization and frequency calculation is done with the matlab code *Z-score_SS_freq.m*. The matlab code needs input as continous-time simulations of any EMP network. The struucture of input data file is like having as many columns as the number of nodes in the network and as many rows as the number of parameter sets used to simulate the model. Each row represents a stable steady state. The matlab code first normalizes each data column using Z- Score standardization, discretizes the data columns, and finally calculates the frequency of each steady state. It also comapres the frequency distribution of each network-perturbed data file with the corresponding wild-type file. The output of the code is a data file that contains steady states of the perturbed file and that of wild type and and their respective frequencies. 

Frustration, hybridness and multistability can be calculated using *frustration.R* , *hybridness.R*, *multistability.R* and *J_metric.R* respectively.

#### Network topology metrics:
Loop data needs to be generated using *dataGen.R* script. All the network metrics, including hiloops can be calculated using *NetworkMetric.R* file. 

The R scripts have detailed comments describing their usage.

### Requirements
Python(Tested on Version 3.8.5)
R (Tested on Version 4.1.2)
Matlab (Tested on Version R2018b)
