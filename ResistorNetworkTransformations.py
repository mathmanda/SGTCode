from scipy import linalg
from sympy import *


def GraphtoWeighted(G):
    A = G.adjacency_matrix()
    G = Graph(A,format = 'weighted_adjacency_matrix', loops = True, multiedges = True)
    return G
        
def series(G,v):
    if G.degree(v) ==2:
        edges = G.edges_incident(v)
        edge1 = edges[0]
        edge2 = edges[1]
        vert1 = list(edge1[0:2])
        vert1.remove(v)
        vert2 = list(edge2[0:2])
        vert2.remove(v)
        newlabel = edge1[2] + edge2[2]
        newedge = (vert1[0], vert2[0])
        print newedge
        G.add_edge(newedge, label = str(newlabel))
        G.delete_vertex(v)
        return G
    else: 
        return G
    
def YDelta(G,v):
    if G.degree(v) ==3:
        edges = G.edges_incident(v)
        #Y-edges
        edge1 = edges[0]
        edge2 = edges[1]
        edge3 = edges[2]
        #Y-edge weights
        w1 = edge1[2]
        w2 = edge2[2]
        w3 = edge3[2]
        #Delta-edge weights
        s = w1*w2 + w2*w3 + w3*w1
        nw1 = s/w1
        nw2 = s/w2
        nw3 = s/w3
        
        #Which vertices do we care about? 
        vert1 = list(edge1[0:2])
        vert1.remove(v)
        vert2 = list(edge2[0:2])
        vert2.remove(v)
        vert3 = list(edge3[0:2])
        vert3.remove(v)
        print vert1, vert2, vert3
        
        #Delta edges
        ne3 = (vert1[0], vert2[0])
        ne2 = (vert1[0], vert3[0])
        ne1 = (vert3[0], vert2[0])
        
        #Add new edges
        G.add_edge(ne1, label = str(nw1))
        G.add_edge(ne2, label = str(nw2))
        G.add_edge(ne3, label = str(nw3))
        
        #Delete old edges
        G.delete_vertex(v)
        return G
    else: 
        return G
