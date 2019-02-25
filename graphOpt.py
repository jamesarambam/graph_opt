"""
++++++++++++++++++++++++++++++++++
Author : James Arambam
Date   : 28 Jan 2016
++++++++++++++++++++++++++++++++++
"""
# =============================================================================== #
import sys
import os
import readBinaryWcspSP  as rwcspBin_minSP
import ctypes
import numpy as np
import auxLib as ax
import parameters as param

# =============================================================================== #

# ------------------------------ Parameters ------------------------------ #
ppath = os.getcwd()+"/"
conVioThres = param.Constraint_Violation_Threshold
outerConvWindow = param.Outer_Convergence_Window
fname = sys.argv[1][5:len(sys.argv[1]) - 5]
opFile = "Result/"+fname+"/"+fname+".log"
ax.opFile = opFile
zeroThreshold = param.zeroThreshold
rcdcopC = ctypes.cdll.LoadLibrary("lib/./rcdcop.so")
rcdcopC.gateWay.restype = ctypes.c_double
rcdcopC.test.restype = ctypes.c_double

# --------------------------------- Methods --------------------------------- #

def addDummyQP():

    global file, dummyQPedges
    dummyQPedges = 0
    new_node = len(G.nodes())
    for n in G.nodes():
        if len(G.node[n]['Nbr_q']) == 0:
            G.add_node(new_node)
            G.add_edge(n, new_node)
            dummyQPedges += 1
            # print getDomain(file)
            # exit()
            G.node[new_node]['Domain'] = 1 #getDomain(file)
            G.node[new_node]['Nbr_l'] = []
            G.node[new_node]['Nbr_q'] = [n]
            G.node[n]['Nbr_q'] = [new_node]
            G[n][new_node]['Potential'] = {}
            for xi in range(G.node[n]['Domain']):
                for xj in range(G.node[new_node]['Domain']):
                    G[n][new_node]['Potential'][(xi, xj)] = [0]
            G[n][new_node]['Type'] = 'QP'
            new_node += 1
    return G

def thetaMax():

    temp = []
    for e in G.edges():
        for xi in range(G.node[e[0]]['Domain']):
            for xj in range(G.node[e[1]]['Domain']):
                temp.append(theta(e[0], e[1], xi, xj))
    return max(temp)

def thetaMin():

    minP = []
    for e in G.edges():
        for xi in range(G.node[e[0]]['Domain']):
            for xj in range(G.node[e[1]]['Domain']):
                minP.append(theta(e[0], e[1], xi, xj))
    return min(minP)

def theta(i, j, xi, xj):
    if i < j:
        return sum(G[i][j]['Potential'][(xi, xj)])
    else:
        return sum(G[j][i]['Potential'][(xj, xi)])

def thetaHat(i, j, xi, xj):

    global maxTheta, minTheta
    val = (theta(i, j, xi, xj) - minTheta + 1) / float(maxTheta - minTheta)
#     val = ( G[i][j]['Potential'][(xi,xj)] ) / float( maxTheta - minTheta )
    return val

def get_LP_QP():

    lp = []
    qp = []
    for e in G.edges():
        if G[e[0]][e[1]]['Type'] == 'LP':
            lp.append(e)
        if G[e[0]][e[1]]['Type'] == 'QP':
            qp.append(e)
    return lp, qp

def get(fname, edge, defCost, st, en):

    pot2 = {}
    with open(fname) as f:
        for i, line in enumerate(f):
            if st <= i <= en:
                temp2 = line.split(" ")
                tmp = [int(temp2[i]) for i in range(0, len(temp2) - 1)]
                pot2[tuple(tmp)] = int(temp2[len(temp2)-1].strip())
    pot2['defCost'] = defCost
    return {tuple(edge) : pot2}

def getNaryPotentials(file):

    pot = {}
    counter = 2
    with open(file) as f:
        for i, x in enumerate(f):
            if i >= 2 and i == counter:
                i_tmp = i
                temp = x.split(" ")
                ary = int(temp[0])
                if ary > 0:
                    edge = [int(temp[n]) for n in range(1, len(temp) - 2)]
                    defCost = int(temp[len(temp) - 2])
                    pot.update(get(file, edge, defCost, i+1, i+int(temp[len(temp) - 1].strip())))
                counter = int(temp[len(temp) - 1].strip()) + i_tmp + 1
    return pot

def minimizeObj():

    global minObj, maxThetaObj
    minObj = True
    maxThetaObj = thetaMax() + 5

    for e in G.edges():
        for xi in range(G.node[e[0]]['Domain']):
            for xj in range(G.node[e[1]]['Domain']):
                if len(G[e[0]][e[1]]['Potential'][(xi, xj)]) > 1:
                    temp = sum(G[e[0]][e[1]]['Potential'][(xi, xj)])
                    G[e[0]][e[1]]['Potential'][(xi, xj)] = [temp]

    for e in G.edges():
        for xi in range(G.node[e[0]]['Domain']):
            for xj in range(G.node[e[1]]['Domain']):
                G[e[0]][e[1]]['Potential'][(xi, xj)][0] = maxThetaObj - G[e[0]][e[1]]['Potential'][(xi, xj)][0]

def init():

    # ------------ Init ------------ #
    if os.path.exists(ax.opFile):
        os.remove(ax.opFile)
    ax.createDir(ppath+"Result/", fname)
    ax.createDir(ppath+"Result/"+fname+"/", "pklDumps")

    # ------------------------------ #

def graphDetails(G):

    e = len(G.edges())
    n = len(G.nodes())
    d = float((2 * e)) / (n * (n - 1))
    return n, e, round(d, 2)

def getThetas(tmp, tmp2):

    for (i, j) in G.edges():
        for xi in range(G.node[i]['Domain']):
            for xj in range(G.node[j]['Domain']):
                tmp[i][j][xi][xj] = theta(i, j, xi, xj)
                tmp[j][i][xj][xi] = theta(i, j, xi, xj)
                tmp2[i][j][xi][xj] = thetaHat(i, j, xi, xj)
                tmp2[j][i][xj][xi] = thetaHat(i, j, xi, xj)
    return tmp, tmp2

def getLogFiles():

    file = "file.log"
    with open(file) as f:
        for line in f:
            ax.write(line)

# -------------------------------- Main() -------------------------------- #

def main():

    global file, G, minObj, LP, QP, fname, MaxItr, totRuntime, naryPotential
    global pklLoc, tbarAdjust, dummyQPedges, ObjAdjust, tSP, lpO, qpO
    global maxTheta, minTheta
    st = ax.now()
    init()
    MaxItr = 100
    minObj = False
    totRuntime = 0.0
    file = str(sys.argv[1])
    G = rwcspBin_minSP.Main(file, param.SpanningTree)
    n, e, d = graphDetails(G)
    naryPotential = getNaryPotentials(file)

    # --------- Dump for future analysis if required --------- #
    ax.dumpDataStr("Result/"+fname+"/pklDumps/G", G)
    tbarAdjust = len([keys for keys in naryPotential])
    lpO, qpO = get_LP_QP()

    # ------------------------------------------------------ #
    ax.writeln("# -------------------- START --------------------- #")
    ax.writeln("Graph > ")
    ax.writeln("Nodes : " + str(n))
    ax.writeln("Edges : " + str(e))
    ax.writeln("Density : " + str(d))
    ax.writeln("QP Edges : "+str(round(float(len(qpO)*100)/(len(lpO)+len(qpO)), 2))+" %")
    ax.writeln("LP Edges : "+str(round(float(len(lpO)*100)/(len(lpO)+len(qpO)), 2))+" %")
    ax.writeln("")

    # ---------------------------------- #
    G = addDummyQP()
    LP, QP = get_LP_QP()
    param.Optimize = sys.argv[2][1:len(sys.argv[2])]
    if param.Optimize == "min":
        minimizeObj()

    # ---------------------- Go To C ------------------------- #
    maxTheta = thetaMax()
    minTheta = thetaMin()
    Nodes = len(G.nodes())
    dList = map(lambda i :  G.node[i]['Domain'], G.nodes())
    Domain = max(map(lambda i :  G.node[i]['Domain'], G.nodes()))
    thetaM = np.full((Nodes, Nodes, Domain, Domain), 0, dtype = np.double)
    thetaHatM = np.full((Nodes, Nodes, Domain, Domain), 0, dtype=np.double)
    thetaM, thetaHatM = getThetas(thetaM, thetaHatM)
    D = np.array(dList, dtype = np.double)
    nl = np.full((Nodes, Nodes), -1, dtype = np.double)
    nq = np.full((Nodes, Nodes), -1, dtype = np.double)
    for i in G.nodes():
        nbr = G.node[i]['Nbr_q']
        for j in nbr:
            nq[i][j] = j
    for i in G.nodes():
        nbr = G.node[i]['Nbr_l']
        for j in nbr:
            nl[i][j] = j

    # ------------- #
    # Note : Make sure default cost is added in eP
    # ------------- #
    eP = np.full((Nodes, Nodes, Domain, Domain), 0, dtype = np.double)
    for e in naryPotential:
        i = e[0]
        j = e[1]
        for d in naryPotential[e]:
            if d <> 'defCost':
                xi = d[0]
                xj = d[1]
                eP[i][j][xi][xj] = naryPotential[e][d]

    # --------------- Extract Pointers ---------------- #
    Pij_test = np.full((Nodes, Nodes, Domain, Domain), 0, dtype=np.double)
    Pstar_ij_test = np.full((Nodes, Nodes, Domain, Domain), 0, dtype=np.double)
    D_ptr = D.ctypes.data_as(ctypes.POINTER(ctypes.c_double))
    nl_ptr = nl.ctypes.data_as(ctypes.POINTER(ctypes.c_double))
    nq_ptr = nq.ctypes.data_as(ctypes.POINTER(ctypes.c_double))
    eP_ptr = eP.ctypes.data_as(ctypes.POINTER(ctypes.c_double))
    thetaM_ptr = thetaM.ctypes.data_as(ctypes.POINTER(ctypes.c_double))
    thetaHM_ptr = thetaHatM.ctypes.data_as(ctypes.POINTER(ctypes.c_double))
    Pij_ptr = Pij_test.ctypes.data_as(ctypes.POINTER(ctypes.c_double))
    Pstar_ij_test_ptr = Pstar_ij_test.ctypes.data_as(ctypes.POINTER(ctypes.c_double))

    # -------------- Pass it To C --------------------- #
    rcdcopC.gateWay(Nodes, Domain, D_ptr, nl_ptr, nq_ptr, eP_ptr, thetaM_ptr, thetaHM_ptr, Pij_ptr, Pstar_ij_test_ptr, param.Total_Itr, param.maxRuntime)
    getLogFiles()

    # -------------------------------- #
    en = ax.now()
    ax.writeln("Instance : "+str(fname))
    ax.writeln("QP Edges : "+str((len(qpO)*100)/(len(lpO)+len(qpO)))+" %")
    ax.writeln("LP Edges : "+str((len(lpO)*100)/(len(lpO)+len(qpO)))+" %")
    ax.writeln("Total Execution Runtime : "+str(ax.getRuntime(st, en)) + " Seconds")
    ax.writeln("# --------------------- END ---------------------- #")

# =============================================================================== #

if __name__ == "__main__":main()
