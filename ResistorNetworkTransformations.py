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
        #print newedge
        G.add_edge(newedge, label = newlabel)
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
        #print vert1, vert2, vert3
        
        #Delta edges
        ne3 = (vert1[0], vert2[0])
        ne2 = (vert1[0], vert3[0])
        ne1 = (vert3[0], vert2[0])
        
        #Add new edges
        G.add_edge(ne1, label = nw1)
        G.add_edge(ne2, label = nw2)
        G.add_edge(ne3, label = nw3)
        
        #Delete old edges
        G.delete_vertex(v)
        return G
    else: 
        return G

    
def DeltaY(G,v1, v2, v3):
    #Delta-edges
    el3 = [edge for edge in G.edges() if set([edge[0], edge[1]]) == set([v1, v2])]
    el2 = [edge for edge in G.edges() if set([edge[0], edge[1]]) == set([v1, v3])]
    el1 = [edge for edge in G.edges() if set([edge[0], edge[1]]) == set([v3, v2])]
    
    m = min(len(el1), len(el2), len(el3))
    
    if m!=0:
        
        #Delta-edges
        e1 = el1[0]
        e2 = el2[0]
        e3 = el3[0]

        #Delta-edge weights
        w1 = e1[2]
        w2 = e2[2]
        w3 = e3[2]
        #print w2
        
        #print w1, w2, w3
        #Y-edge weights
        s = w1 + w2+ w3
        nw1 = w2*w3/s
        nw2 = w1*w3/s
        nw3 = w1*w2/s
        #print nw1
        
        #Which vertices do we care about? 
        vnew = max(G.vertices())+1

        
        #print v1, v2, v3, vnew
        
        #Y edges
        ne1 = (v1,vnew)
        ne2 = (v2,vnew)
        ne3 = (v3,vnew)
        #print ne1, ne2, ne3
        
        #Add new edges
        G.add_edge(ne1, label = nw1)
        G.add_edge(ne2, label = nw2)
        G.add_edge(ne3, label = nw3)
        
        #Delete old edges
        G.delete_edges([e1, e2, e3])
        return G
    #else: 
        #return G

        
def DeltaY2(G,e1,e2,e3):
    #Delta-edges
    Gsub = Graph()
    Gsub.add_edges([e1,e2,e3])
    
    if e1 in G.edges() and e2 in G.edges() and e3 in G.edges() and Gsub.is_clique():

        #Delta-edge weights
        w1 = e1[2]
        w2 = e2[2]
        w3 = e3[2]
        #print w2
        
        #print w1, w2, w3
        #Y-edge weights
        s = w1 + w2+ w3
        nw1 = w2*w3/s
        nw2 = w1*w3/s
        nw3 = w1*w2/s
        #print nw1
        
        #Which vertices do we care about? 
        oldverts = Gsub.vertices()
        v1 = copy(oldverts)
        v1.remove(e1[0])
        v1.remove(e1[1])
        v1 = v1[0]
        v2 = copy(oldverts)
        v2.remove(e2[0])
        v2.remove(e2[1])
        v2 = v2[0]

        v3 = copy(oldverts)
        v3.remove(e3[0])
        v3.remove(e3[1])
        v3 = v3[0]
        vnew = max(G.vertices())+1

        
        #print v1, v2, v3, vnew
        
        #Y edges
        ne1 = (v1,vnew)
        ne2 = (v2,vnew)
        ne3 = (v3,vnew)
        #print ne1, ne2, ne3
        
        #Add new edges
        G.add_edge(ne1, label = nw1)
        G.add_edge(ne2, label = nw2)
        G.add_edge(ne3, label = nw3)
        
        #Delete old edges
        G.delete_edges([e1, e2, e3])
        return G
    #else: 
        #return G
