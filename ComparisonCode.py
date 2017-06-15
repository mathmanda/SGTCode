from scipy import linalg
from sympy import *


def Gcreate(i):
    atlas = Graph_atlas()
    G = Graph(atlas[i])
    return G
        
def Graph_atlas():
    from sage.graphs.graph import Graph
    import networkx.generators.atlas
    return networkx.generators.atlas.graph_atlas_g()
#def Fiedlers(G):

def n_Atlas(i):
    Atlas=[]
    for G in graphs.nauty_geng(str(i)):
        if G.is_connected():
            Atlas.append(G)
    return Atlas
    
def maxcn(G):
    #Make lists of edges
    all_edges = list([ (e[0],e[1]) for e in Combinations(G.vertices(),k=2)])
    Gedges = list([ (e[0],e[1]) for e in G.edges()])
    notedges=all_edges
    for e in Gedges:
        notedges.remove(e)
        
    #How many common neighbors does each 'edge' have?       
    dist2 = []
    n = len(G.vertices())
    neighbor = {e:0 for e in notedges}
    
    
    for e in notedges:
        if G.distance(e[0],e[1]) == 2:
            dist2.append(e)
            for m in range(n):
                if G.distance(e[0],m) == 1 and G.distance(e[1],m) == 1:
                    neighbor[e] +=1
    
    #What's the max?
    if len(notedges) == 0:
        returnlist = []
    else:
        mx = max(neighbor.values())
    
        #Now we return only the biggest. 
        returnlist = []
    
        for e in dist2:
            if neighbor[e] == mx:
                returnlist.append(e)   
    
    return returnlist

def lowresist(G, numeric = True, rd = 20):
    #Make lists of edges
    all_edges = list([ (e[0],e[1]) for e in Combinations(G.vertices(),k=2)])
    Gedges = list([ (e[0],e[1]) for e in G.edges()])
    notedges=all_edges
    for e in Gedges:
        notedges.remove(e)
        
    #Laplacian, PseudoInverse
    L = G.laplacian_matrix()
    #L = matrix(QQbar, L)
    
    M = matrix(linalg.pinv(L))

    #M = L.pseudoinverse
    
    #Construct v
    id = matrix.identity(len(G.vertices()))
    returnlist = {}
    for e in notedges:
        ve = id[e[0]] -id[e[1]]
        if numeric == False:
            se = ve*M*ve
            se = nsimplify(se)
        else:
	    	se = round(ve*M*ve,30)
        returnlist[e] = se
        
    returnlist =[(a,b) for b,a in returnlist.items()]
    returnlist = sorted(returnlist)
    return returnlist
    
def comparegraphlist(list, numeric = True ,rd = 20):
    badcutoffscores = []
    exceptions = []
    
    for G in list:
        mxcn = maxcn(G)
        m = len(mxcn)
        lowscorelist = lowresist(G, numeric = numeric, rd = rd)
        scores = [a[0] for a in lowscorelist]
        scoreitems = [a[1] for a in lowscorelist]
        n = len(G.vertices())
        if m<len(scores) and scores[m-1] == scores[m]:
            badcutoffscores.append(G)
        else:
            scoreitems = [scoreitems[i] for i in range(m)]
            if set(mxcn) != set(scoreitems):
                exceptions.append(G)
    print len(badcutoffscores), " bad cutoff scores"
    print len(exceptions), " exceptions"
    
    
    return [badcutoffscores, exceptions]
    
def comparegraphsofsize(n, numeric = True, rd = 20):
    biglist = n_Atlas(n)
    returnlist = comparegraphlist(biglist, numeric = numeric, rd = rd)
    return returnlist
    

def cubemethod(G):
    #Make lists of edges
    n = len(G.vertices())
    hcn = []
    for i in range(n):
    	for j in range(i):
    		if G.distance(i,j) == 3:
    			hcn.append((i,j))
    
    A = G.adjacency_matrix()
    A3 = A**3
    pathlengths = [ A3[i][j] for (i,j) in hcn]
    if len(hcn) == 0:
        returnlist = []
    #What's the max?

    else:
        mx = max(pathlengths)
    
        #Now we return only the biggest. 
        returnlist = []
    
        for (i,j) in hcn:
            if A3[i][j] == mx:
                returnlist.append((i,j))   
    
    return returnlist
	
def lowresistb(G, numeric = True, rd = 20):
    #Make lists of edges
    all_edges = list([ (e[0],e[1]) for e in Combinations(G.vertices(),k=2)])
    Gedges = list([ (e[0],e[1]) for e in G.edges()])
    notedges=all_edges
    for e in Gedges:
        notedges.remove(e)

    candidates = [e for e in notedges if G.distance(e[0],e[1]) % 2 != 0]
    #Laplacian, PseudoInverse
    L = G.laplacian_matrix()
    #L = matrix(QQbar, L)
    
    M = matrix(linalg.pinv(L))

    #M = L.pseudoinverse
    
    #Construct v
    id = matrix.identity(len(G.vertices()))
    returnlist = {}
    for e in candidates:
        ve = id[e[0]] -id[e[1]]
        if numeric == False:
	        se = ve*M*ve
        else:
	    	se = round(ve*M*ve,30)
        returnlist[e] = se
        
    returnlist =[(a,b) for b,a in returnlist.items()]
    returnlist = sorted(returnlist)
    return returnlist

                
def comparebipartitegraphs(n, numeric = True, rd = 20):
    biglist = n_Atlas(str(n) + " -b")
    #small=[graph for graph in biglist if graph.is_bipartite()]
    returnlist = comparebipartitelist(biglist, numeric = numeric, rd = rd)
    return returnlist
    
def comparebipartitelist(list, numeric = True ,rd = 20):
    badcutoffscores = []
    exceptions = []
    nopredict =[]
    for G in list:
        mxcn = cubemethod(G)
        m = len(mxcn)
        lowscorelist = lowresistb(G, numeric = numeric, rd = rd)
        scores = [a[0] for a in lowscorelist]
        scoreitems = [a[1] for a in lowscorelist]
        n = len(G.vertices())
        if len(scores) == 0:
        	nopredict.append(G)
        elif m ==0:
        	exceptions.append(G)
        elif m<len(scores) and scores[m-1] == scores[m]:
            badcutoffscores.append(G)
        else:
            scoreitems = [scoreitems[i] for i in range(m)]
            scoreitems = [ [e[0],e[1]] for e in scoreitems]
            for e in scoreitems:
				e.sort()
            set3paths = [[e[0],e[1]] for e in mxcn]
            for e in set3paths:
				e.sort()
            set3paths.sort()
            scoreitems.sort()
            print 'paths', set3paths
            print 'resistance', scoreitems
            print set3paths == scoreitems

            if set3paths != scoreitems:
                exceptions.append(G)
        		
    print len(badcutoffscores), " bad cutoff scores"
    print len(exceptions), " exceptions"
    print len(nopredict), " with no possible links"
    
    return [badcutoffscores, exceptions,nopredict]
        
    
def minmaxedges(list):
	min = 100000
	max = 0 
	for G in list:
		m = G.num_edges()
		if m<min:
			min = m         
		elif m>max:
			max = m
	return[min,max]
    

def elementwiseprod( M, N):
    #assert(M.parent() == N.parent())
    nc, nr = M.ncols(), M.nrows()
    A = copy(M.parent().zero_element())
    for r in xrange(nr):
        for c in xrange(nc):
            A[r,c] = M[r,c]*N[r,c]
    return A


def effresistance_matrix(G, numeric = True, rd = 20, zeroadjacents = True):
    #Make lists of edges
    #all_edges = list([ (e[0],e[1]) for e in Combinations(G.vertices(),k=2)])
    #Gedges = list([ (e[0],e[1]) for e in G.edges()])
    #notedges=all_edges
    #for e in Gedges:
    #    notedges.remove(e)
        
    #Laplacian, PseudoInverse
    L = G.laplacian_matrix()
    #L = matrix(QQbar, L)
    
    M = matrix(linalg.pinv(L))
    
    d = matrix(M.diagonal())
    d = d.transpose()
    
    e = allonesvector(G.num_verts())
    
    R = d*e.transpose() + e*d.transpose() - 2* M
    
    m = len(G.spanning_trees())
    S = R * m
    
    n = G.num_verts()
    for i in range(n):
        for j in range(n):
            S[i,j] = nsimplify(S[i,j])
            
    if zeroadjacents == True:
        B = allonesmatrix(n) - G.adjacency_matrix()
        
        S = elementwiseprod(S,B)    

    #M = L.pseudoinverse
    
    #Construct v
    #id = matrix.identity(len(G.vertices()))
    #returnlist = {}
    #for e in notedges:
    #    ve = id[e[0]] -id[e[1]]
    #    if numeric == False:
    #        se = ve*M*ve
    #        se = nsimplify(se)
    #    else:
    #        se = round(ve*M*ve,30)
    #    returnlist[e] = se
        
    #returnlist =[(a,b) for b,a in returnlist.items()]
    #returnlist = sorted(returnlist)
    return S, m 


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

    


def allonesvector(n):
    v = matrix(n,1)
    #A = copy(M.parent().zero())
    for r in xrange(n):
            #try:
            #    A[r,c] = nsimplify(M[r,c]*N[r,c])
            #except:
            v[r,0] = 1
    return v

def resistance_distance(G, (a,b), numeric = True, rd = 20):
    #Make lists of edges
        
    #Laplacian, PseudoInverse
    L = G.laplacian_matrix()
    #L = matrix(QQbar, L)
    
    M = matrix(linalg.pinv(L))

    #M = L.pseudoinverse
    
    #Construct v
    id = matrix.identity(len(G.vertices()))
    ve = id[a] -id[b]
    if numeric == False:
        se = ve*M*ve
        se = nsimplify(se)
    else:
        se = round(ve*M*ve,30)
        
    return se


'''    
def WBalg(G,listofnodelists):
    #This takes a graph and an equitable partition (given as a list of nodelists)
    #and should return the 'S' matrix, and the resulting product S**(-1)*L*S
    
    n = G.num_verts()
    L = G.laplacian_matrix()
    S = []
    for nodeset in listofnodelists:
        nodevec = []
        for i in range(n):
            if i in nodeset:
                nodevec.append(1)
            else:
                nodevec.append(0)
        S.append(nodevec)
    
    V = VectorSpace(QQ,n)
    W = V.subspace(S)
    S1 = list(W.basis_matrix())
    Wperp = W.complement()
    S2 = list(Wperp.basis_matrix())
    S3 = S2+S1
    S = matrix(S3).transpose()
    result = S**(-1)*L*S
    return [S, result]

def Adjacency(G,listofnodelists):
    #This takes a graph and an equitable partition (given as a list of nodelists)
    #and should return the 'S' matrix, and the resulting product S**(-1)*L*S
    
    n = G.num_verts()
    A = G.adjacency_matrix()
    S = []
    for nodeset in listofnodelists:
        nodevec = []
        for i in range(n):
            if i in nodeset:
                nodevec.append(1)
            else:
                nodevec.append(0)
        S.append(nodevec)
    
    V = VectorSpace(QQbar,n)
    W = V.subspace(S)
    S1 = list(W.basis_matrix())
    Wperp = W.complement()
    B = Wperp.basis_matrix()
    M,G = B.gram_schmidt()
    S2 = list(M)
    S3 = S2+S1
    S = matrix(S3).transpose()
    result = S**(-1)*A*S
    return [S, result]

def Adjacency2(G,listofnodelists):
    #This takes a graph and an equitable partition (given as a list of nodelists)
    #and should return the upper block matrix, and the lower block matrix in S**-1 A S
    
    n = G.num_verts()
    A = G.adjacency_matrix()
    K = [i for i in range(n)]
    S = []
    block2set=[]
    for nodeset in listofnodelists:
        nodevec = []
        block2set.append(nodeset[0])
        K.remove(nodeset[0])
        for i in range(n):
            if i in nodeset:
                nodevec.append(1)
            else:
                nodevec.append(0)
        S.append(nodevec)
        
    
    V = VectorSpace(QQbar,n)
    W = V.subspace(S)
    S1 = list(W.basis_matrix())
    Wperp = W.complement()
    B = Wperp.basis_matrix()
    M,G = B.gram_schmidt()
    S2 = list(M)
    S3 = S2+S1
    S = matrix(S3).transpose()
    result = S**(-1)*A*S
    princsmA = A.matrix_from_rows_and_columns(K,K)
    p = len(listofnodelists)
    J1 = [i for i in range(n-p)]
    J2 = [i  for i in range(n-p,n)]
    
    upperblock = result.matrix_from_rows_and_columns(J1,J1)
    lowerblock = result.matrix_from_rows_and_columns(J2,J2)
    
    return [result,princsmA, upperblock,lowerblock]



G = Graph()
G.add_edges([(0,1),(0,2),(1,2)])
G.add_edges([(0,3),(1,6),(2,9),(0,4),(1,7),(2,10),(0,5),(1,8),(2,11)])
G.add_edges([(3,4),(4,5),(3,5),(6,7),(7,8),(8,6),(9,10),(9,11),(10,11)])
G.add_edges([(3,9),(9,7),(7,5),(5,11),(11,6),(6,4),(4,10),(8,10),(8,3)])
G.show()

A=G.adjacency_matrix()

k = CyclotomicField(9)
z = k.gen()
w = z**3
S3 = matrix([[1, 1, 1],[1,w,w**2],[1,w**2,w]])
Sone = matrix([[1, 1, 1],[1, 1, 1],[1, 1, 1]])
B = diagonal_matrix([1,w,w**2])
SB1 = block_matrix([[Sone*B**0, B**(1)*Sone*B**0,B**2*Sone*B**0]])
SB2 = block_matrix([[Sone*B**1, z*B**(1)*Sone*B**1,z**2*B**2*Sone*B**1]])
SB3 = block_matrix([[Sone*B**2, z**2*B**(1)*Sone*B**2,z*B**2*Sone*B**2]])
SB = block_matrix([[SB1],[SB2],[SB3]])

#SS = matrix([[1, 1, 1,1,1,1,1,1,1],[1,z,  z**2,z**3,z**4,z**5,z**6,z**7,z**8],[1,z**2,z**4,z**6,z**8,z**1,z**3,z**5,z**7],[1,z**3,z**6,z**9,z**3,z**6,z**9,z**3,z**6],[1,z**4,z**8,z**3,z**7,z**2,z**6,z**1,z**5],[1,z**5,z**1,z**6,z**2,z**7,z**3,z**8,z**4],[1,z**6,z**3,z**9,z**6,z**3,z**9,z**6,z**3],[1,z**7,z**5,z**3,z**1,z**8,z**6,z**4,z**2],[1,z**8,z**7,z**6,z**5,z**4,z**3,z**2,z**1]])

partition = [[0,5],[1,4],[2,3]]
result = Adjacency2(G,partition)
print latex(G.adjacency_matrix())
print 'new'
print result[0]
print '$'
print latex(result[3])
print '$&$ '
print latex(result[2])
print '$&$ '
print latex(result[1])
'''


            