"""
++++++++++++++++++++++++++++++++++
Author : James Arambam
Date   : 28 Jan 2016

++++++++++++++++++++++++++++++++++
"""

# ===============================================================================

import networkx as nx
import os
from pprint import pprint
import random as rnd
import auxLib as ax

# --------------------- Variables ------------------------------

ppath = os.getcwd()+"/"         # Project Path Location

# --------------------- Methods ------------------------------ #

def addPotential(G, e0, e1, st, en, defCost):

    G.add_edge(e0, e1)
    xixj = []
    f = open(fname)
    if not G[e0][e1].has_key('Potential'):
        G[e0][e1]['Potential'] = {}
    for xi in range(G.node[e0]['Domain']):
        for xj in range(G.node[e1]['Domain']):
            if not G[e0][e1]['Potential'].has_key((xi, xj)):
                G[e0][e1]['Potential'][(xi, xj)] = []
    for i, x in enumerate(f):
        if st <= i <= en:
            temp = x.split(" ")
            xixj.append((int(temp[0]), int(temp[1])))
            G[e0][e1]['Potential'][(int(temp[0]), int(temp[1]))].append(int(temp[2]))
    f.close()
    # --------- Add Default Costs ---------- #
    for xi in range(G.node[e0]['Domain']):
        for xj in range(G.node[e1]['Domain']):
            if not (xi, xj) in xixj:
                G[e0][e1]['Potential'][(xi, xj)].append(int(defCost))
    return G


def addPotentialRev(G, e0, e1, st, en, defCost):

    e0, e1 = e1, e0
    G.add_edge(e0, e1)
    xixj = []
    f = open(fname)
    if not G[e0][e1].has_key('Potential'):
        G[e0][e1]['Potential'] = {}
    for xi in range(G.node[e0]['Domain']):
        for xj in range(G.node[e1]['Domain']):
            if not G[e0][e1]['Potential'].has_key((xi, xj)):
                G[e0][e1]['Potential'][(xi, xj)] = []
    for i, x in enumerate(f):
        if st <= i <= en:
            temp = x.split(" ")
            xixj.append((int(temp[1]), int(temp[0])))
            G[e0][e1]['Potential'][(int(temp[1]), int(temp[0]))].append(int(temp[2]))
    f.close()
    # --------- Add Default Costs ---------- #
    for xi in range(G.node[e0]['Domain']):
        for xj in range(G.node[e1]['Domain']):
            if not (xi, xj) in xixj:
                G[e0][e1]['Potential'][(xi, xj)].append(int(defCost))
    return G


def addDomain(G, oVar, n, st, en):

    G.node[n]['Domain'] = {}
    G.node[n]['Potential'] = {}
    dCount = 0
    with open(fname) as f:
        for i, line in enumerate(f):
            if st <= i <= en:
                temp = line.split(" ")
                xi = []
                for j in range(len(temp) - 1):
                    xi.append(int(temp[j]))
                G.node[n]['Domain'][dCount] = {}
                for k in range(len(xi)):
                    G.node[n]['Domain'][dCount][oVar[k]] = xi[k]
                G.node[n]['Potential'][dCount] = int(temp[len(temp) - 1])
                dCount += 1


def getSharedVar(di, dj):

    sv = list(set(di[0].keys()).intersection(set(dj[0].keys())))
    return sv


def sharedValueViolate(i, j, xi, xj, sv):

    for d in sv:
        if i[xi][d] <> j[xj][d]:
            return False
    return True


def getLocalMax(fname, st, en):

    tempMax = 0
    with open(fname) as f:
        for i, line in enumerate(f):
            if st <= i <= en:
                temp2 = line.split(" ")
                tmp = float(temp2[len(temp2) - 1].strip())
                if tempMax < tmp:
                    tempMax = tmp
    return tempMax


def getGlobalMax(file):

    counter = 2
    maxP = 0
    with open(file) as f:
        for i, x in enumerate(f):
            if i >= 2 and i == counter:
                i_tmp = i
                temp = x.split(" ")
                maxTmp = int(temp[len(temp) - 2])
                if maxP < maxTmp:
                    maxP = maxTmp
                max2 = getLocalMax(file, i+1, i+int(temp[len(temp) - 1].strip()))
                if maxP < max2:
                    maxP = max2
                counter = int(temp[len(temp) - 1].strip()) + i_tmp + 1
    return maxP


def read(fname, G):

    f = open(fname)
    for i, x in enumerate(f):
        if i == 0:
            temp = x.split(" ")
            for n in range(int(temp[1])):
                G.add_node(n)
        if i == 1:
            temp = x.split(" ")
            for n in G.nodes():
                G.node[n]['Domain'] = int(temp[n])
    f.close()
    f = open(fname)
    counter = 2
    for i, x in enumerate(f):
        if i >= 2 and i == counter:
            i_tmp = i
            temp = x.split(" ")
            aryity = int(temp[0])
            convar = [int(temp[j]) for j in range(1, len(temp) - 2)]
            if aryity == 2:
                if sorted(convar) == convar:
                    G = addPotential(G, int(temp[1]), int(temp[2]), i + 1, i + int(temp[4].strip()), int(temp[3]))
                else:
                    G = addPotentialRev(G, int(temp[1]), int(temp[2]), i + 1, i + int(temp[4].strip()), int(temp[3]))
            else:
                ax.writeln("*************Nary Constraints in file !")
                ax.writeln("Exiting !")
                exit()
            counter = int(temp[len(temp) - 1].strip()) + i_tmp + 1
    f.close()
    return G


def getLP(G, n):

    # e = []
    # intrv = len(G.nodes()) / n
    # for i in range(0, len(G.nodes()), intrv):
    #
    #     source = i  #rnd.randint(0, len(G.nodes())-1)
    #     T = ax.prim_mst_Source(G, source)
    #     e.extend(T.edges())
    # e = list(set(e))
    # return e

    rnd.seed(100)


    e = []
    for i in range(n):
        source = rnd.randint(0, len(G.nodes())-1)
        T = ax.prim_mst_Source(G, source)
        e.extend(T.edges())
    e = list(set(e))
    return e


def add_LP_QP(G):

    LP_tmp = getLP(G, tSP)
    QP_tmp = list(set(G.edges()) - set(LP_tmp))

    # if param.QP is False:
    #     LP_tmp = G.edges()
    #     QP_tmp = []
    # if param.QP is True:
    #     LP_tmp = []
    #     QP_tmp = G.edges()

    # ----------- Compute LP ------------ #
    for e in LP_tmp:
        G[e[0]][e[1]]['Type'] = 'LP'

    # ----------- Compute QP ------------ #
    for e in QP_tmp:
        G[e[0]][e[1]]['Type'] = 'QP'
    return G


def find_QP_Nbr(G, n):

    temp = []
    for nbr in G.neighbors(n):
        if G[n][nbr]['Type'] == "QP":
            temp.append(nbr)
    return temp


def find_LP_Nbr(G, n):

    temp = []
    for nbr in G.neighbors(n):
        if G[n][nbr]['Type'] == "LP":
            temp.append(nbr)
    return temp


def xtraAttr(G):

    G = add_LP_QP(G)
    for n in G.nodes():
        G.node[n]['Nbr_l'] = find_LP_Nbr(G, n)
        G.node[n]['Nbr_q'] = find_QP_Nbr(G, n)
    return G


def print_graph(G):

    for n in G.nodes():
        print "Node : ",
        print n
        pprint(G.node[n])
    for e in G.edges():
        print "Edge : ",
        print e
        pprint(G[e[0]][e[1]])


def Main(file, tmpSP):

    global fname, instID, constraints, globalMax, tSP

    tSP = tmpSP
    constraints = {}
    fname = ppath + file
    instID = file[5:len(file) - 5]
    ax.opFile = "Result/" + instID + "/" + instID + "_" + str(tSP) + ".log"
    globalMax = getGlobalMax(fname)
    G = nx.Graph()
    G = read(fname, G)
    G = xtraAttr(G)
    return G


# --------------------------------------------------------------

#if __name__ == '__main__':Main("data/test22.wcsp")
