import numpy as np
import pandas as pd
# import matplotlib.pyplot as plt
from numba import jit
import os
# from mpi4py import MPI
import sys
# import networkx as nx
import csv
# from scipy.stats import pearsonr
import glob

# comm = MPI.COMM_WORLD # this is a communicators
# rank = comm.Get_rank() # this tells you the id of the current instance

# num_id = 1 + rank
# file = str(sys.argv[num_id]) 
# print("my rank is " + str(rank) + " and my file is " + file)
# net = file.partition('_del')[0]

def compute_hamming(rand_net_file, wt_net_file, column):
    rand_net = pd.read_table(rand_net_file, sep = " ")
    wt_net = pd.read_table(wt_net_file, sep = " ")
    hamming = np.sum([0 if (rand_net.iloc[i,column] == wt_net.iloc[i,column]) else 1 for i in range(rand_net.shape[0])])
    return hamming

def topo2interaction(file, rule = 0, weight=100):
    table = pd.read_table(file, sep=" ")
    tmp = table[["Source", "Target"]].values.reshape(-1)
    node_labels = sorted(list(set(tmp)))
    N = len(node_labels)
    idx2label = dict(enumerate(node_labels))
    label2idx = {v: k for k, v in idx2label.items()}
    T2J_dict = {1: 1, 2: -1}
    if rule != 0:
        T2J_dict[rule] = T2J_dict[rule]*weight
    J = np.zeros((N, N), np.float64)
    for u, v, t in table.values:
        j = label2idx[u]
        i = label2idx[v]
        J[i, j] = T2J_dict[t]
    return J, np.array(node_labels)

@jit(nopython=True)
def _run(J, J_pseudo,initial_pert_ss,maxT,mode, can_be_updated):
    s = initial_pert_ss.copy()
    s_check = np.sign(J_pseudo@s)
    if (np.all(s_check==s)):
        convergence = True
        return convergence, s
    convergence = False
    CT = 0
    for ct in range(maxT):
        n = ct//100
        if mode == "async":
            k = np.random.choice(can_be_updated)
            sk_new = np.sign((J_pseudo@s)[k])
            if s[k] != sk_new:
                s[k] = sk_new
        else:
            s = np.sign(J_pseudo@s)
        s_check = np.sign(J_pseudo@s)
        if (np.all(s_check==s)):
            CT = n
            convergence = True
            return convergence, s
    return convergence, s

def state_run(net_file,states_file,rule=0,mode="async",maxT=1000):
    state_df = pd.read_csv(states_file)
    # try:
    flag_state_df = state_df[state_df['flag'] == 1]
    flag_state_df = flag_state_df[flag_state_df['Avg0'].notna()]
    states = flag_state_df['states'].tolist()
    freq = flag_state_df['Avg0'].tolist()
    # except:
    #     flag_state_df = state_df
    #     flag_state_df = flag_state_df[flag_state_df['async_0_Mean'].notna()]
    #     states = flag_state_df['State'].tolist()
    #     freq = flag_state_df['async_0_Mean'].tolist()
    J,node_labels = topo2interaction(net_file,rule=rule)
    num_nodes, M = J.shape
    can_be_updated = np.array([a for a, b in enumerate(node_labels)])
    J_pseudo = np.identity(num_nodes) + 2 * J
    w=0
    coherence_state = []
    for state in states:
        #print(w)
        state_s = state[1:len(state)-1]
        state_n = []
        for s in state_s:
            state_n.append(int(s))
        state_n = np.array(state_n, dtype=np.float64)
        indices_one = state_n == 1
        indices_zero = state_n == 0
        state_n[indices_one] = 1 
        state_n[indices_zero] = -1
        steady_state = state_n
        coherence_node = 0
        for i in range(len(state_n)):
            state_pert = state_n.copy()
            state_pert[i] = state_pert[i]*(-1)
            for j in range(100):
                convergence,s = _run(J, J_pseudo, state_pert,maxT, mode, can_be_updated)
                if convergence:
                    if np.all(s == steady_state):
                        coherence_node = coherence_node + 1
        coherence_state.append(coherence_node/(len(state_n)*100))
    df = pd.DataFrame({'states':states, 'coherence':coherence_state})
    return df



topoFiles = glob.glob('*.topo')
nets = [s.replace(".topo", "") for s in topoFiles]
stateFiles = [s.replace(".topo", "_finFlagFreq.csv") for s in topoFiles]

for i in range(len(topoFiles)):
    net_file = topoFiles[i]
    net = nets[i]
    states_file = stateFiles[i]
    df = state_run(net_file, states_file)
    df.to_csv(net + "_coherence.csv")
    print(net)

