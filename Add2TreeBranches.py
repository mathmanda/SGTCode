def AddLinear2TreeBranch(G,v1,v2,l):
    vert1 = v1
    vert2 = v2
    for i in range(l):
        vnew = max(G.vertices()) + 1
        G.add_edges([[vert1, vnew],[vert2,vnew]])
        vert1 = vert2
        vert2 = vnew
    return G
        
        
def Pinwheel(n):
    
    G = Graph()
    G.add_edges([[0,1],[1,2],[2,0]])
    G = AddLinear2TreeBranch(G,0,1,n)
    G = AddLinear2TreeBranch(G,1,2,n)
    G = AddLinear2TreeBranch(G,2,0,n)
    return G
        

def APinwheel(n):
    
    G = Graph()
    G.add_edges([[0,1],[1,2],[2,0]])
    G = AddLinear2TreeBranch(G,0,1,n)
    G = AddLinear2TreeBranch(G,1,2,n)
    G = AddLinear2TreeBranch(G,0,2,n)
    return G
            