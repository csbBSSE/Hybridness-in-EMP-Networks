"""
Created on Tue 12.December.2019

@author: Burhan, Kishore

code_description: This code contains a function that requires a topo file as an input and outputs the number of postive
and negative cycles, the edges that they traverse through, and the length of each feedback loop.

"""
import numpy as np 
import os
import networkx as nx
import sys
import pandas as pd
import glob


def networkNcycles(name):
    G = nx.DiGraph()
    if '.topo' in name:
        f = open(str(name),'r')
        content = f.read()
        f.close()
        content = content.split('\n')
        content.remove(content[0])
        content = [content[i].split() for i in range(len(content))]
        cdict = {}
        content = [x for x in content if not x == []]
        for i in range(len(content)):
            G.add_edge(content[i][0], content[i][1])
            if content[i][2] == '2':
                cdict[str(content[i][0] + ','+content[i][1])] = -1
            elif content[i][2] == '1': 
                cdict[str(content[i][0] + ','+content[i][1])] = 1
        #print(cdict)
        cycles= list(nx.cycle_basis(G.to_undirected()))
        return [G, cycles, cdict]
        # counting the number of positive and negative loop
    else:
        raise Exception("Input should be a topo file")


def loopNEdgeCounter(cycles, cdict):
    p_fbc = 0
    n_fbc = 0
    for i in range(len(cycles)):
        prod = 1
        for j in range(len(cycles[i])):
            prod = prod * cdict[cycles[i][j]+','+cycles[i][(j+1)%len(cycles[i])]]
        if prod == -1:
            print(cycles[i])
            p_fbc = p_fbc+1
        else:
            n_fbc = n_fbc + 1
    return [p_fbc, n_fbc]

def edgeCounter(cycles, cdict, G):
    edict = {'Cycles': [], 'Edge_count':[]};
    for i in range(len(cycles)):
        cyc = ",".join(cycles[i])
        edict['Cycles'].append(cyc)
        edict["Edge_count"].append(len(cycles[i]))
        #print(edict)
    return(pd.DataFrame(edict))




topo_file_list = glob.glob('*.topo')
my_folder = sys.argv[1]
if not os.path.exists(my_folder):
        os.makedirs(my_folder)
for i in topo_file_list:
    list1 = networkNcycles(i)
    name = i.replace(".topo", "")
    cycles = list1[1]
    G = list1[0]
    cdict = list1[2]
    edict = edgeCounter(cycles, cdict, G)
    edict.to_csv(("undirectedLoops/" + name + "_undirectedLoops.csv"), index=False)
    print(i)
