from scipy import linalg
from sympy import *


def Katzmatrix(G, alpha):
    A = G.adjacency_matrix()
    In = matrix.identity(G.num_verts())
    M = In - alpha * A
    K =  M.inverse()-In
    Missing = allonesmatrix(G.num_verts()) - A- In
    #print "Missing"
    #print Missing 
    return elemwiseprod(K,Missing)

def Katzmethod(G,alpha):
    #Make lists of edges
    all_edges = list([ (e[0],e[1]) for e in Combinations(G.vertices(),k=2)])
    Gedges = list([ (e[0],e[1]) for e in G.edges()])
    notedges=all_edges
    for e in Gedges:
        notedges.remove(e)
        
    returnlist = {}
    K = Katzmatrix(G,alpha)
    #print K
    for e in notedges:
        se = K[e[0],e[1]]
        returnlist[e] = se
        #print e, se
        
    returnlist =[(a,b) for b,a in returnlist.items()]
    returnlist = sorted(returnlist)
    return returnlist
    
  
def elemwiseprod( M, N):
    nc, nr = M.ncols(), M.nrows()
    A = copy(M.parent().zero())
    for r in xrange(nr):
        for c in xrange(nc):
            #try:
            #    A[r,c] = nsimplify(M[r,c]*N[r,c])
            #except:
            A[r,c] = M[r,c]*N[r,c]
    return A

def allonesmatrix(n):
    M = matrix(n)
    A = copy(M.parent().zero())
    for r in xrange(n):
        for c in xrange(n):
            #try:
            #    A[r,c] = nsimplify(M[r,c]*N[r,c])
            #except:
            A[r,c] = 1
    return A
    
    
        


            