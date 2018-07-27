import numpy as np
import matplotlib.pyplot as plt
import networkx as nx
import pydot as pd

def load_graph(fname,numflag=0):

    G = nx.Graph()
    with open(fname,'r') as f:
        n, m = (int(num) for num in f.readline().split())

        if numflag == 0:
            vtxs = list(range(0,n))
            print(vtxs)
        else:
            vtxs = list(range(1,n+1))
            print(vtxs)

        G.add_nodes_from(vtxs)

        for i in vtxs:
            edges = [(i,int(num)) for num in f.readline().split()]
            print(edges)
            G.add_edges_from(edges)

    return G

def export_graph(fname,G):

    with open(fname,'w') as f:
        nx.drawing.nx_pydot.write_dot(G,f)

if __name__ == '__main__':
    
    # G = load_graph("../graphs/4elt.graph")
    # G = load_graph("fdm.graph")

    G = load_graph("test1.graph",1)

    # export_graph("fdm.dot",G)
    # print("loaded graph")

    plt.subplot(111)

    nx.draw(G,with_labels=True)

    plt.show()
