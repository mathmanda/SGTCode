from scipy import linalg
from sympy import *


def SuperTriangle(h):
    G = Graph()
    ej=0
    prevrowstart = 0
    for i in range(h):
        #This row starts on next vertex after end of last row
        rowstart = ej+1
        #Set up starting vertex of this row, and vertex to attach to in prevrow
        ej = rowstart
        dj = prevrowstart
        for j in range(i+1):
            #Add one triangle
            G.add_edges([[ej, ej+1],[ej,dj],[ej+1,dj]])
            #Advance to next triangle
            dj+=1
            ej+=1
        #Now this row becomes previous row
        prevrowstart = rowstart
    return G
        

            