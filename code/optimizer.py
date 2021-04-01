import math
import time

from config import *


try:
    import numpy as np
    from gurobipy import *

except Exception as e:
    logging.exception("External required libraries importing exception occurred")
    raise e


def run_optimization(protein_seq: str, cpb: list, fitness_values: list, threshold: float) -> str:
    """
    Take amino acid sequence of protein as input and mathematically optimize to DNA sequence. Fully based on MAXCPBstCAI
    function from software implementation of approach described in "Codon Optimization: A Mathematical Programming
    Approach" article. The function maximizes CPB (a parameter depending on the occurrence of codon pairs) when the CAI
    (Codon Adaptation Index) does not fall below the specified value (threshold).

    :param protein_seq: sequence of protein for optimization
    :param cpb: codon pairs bias table
    :param fitness_values: fitness values list
    :param threshold: the min value for CAI
    :return: string with DNA sequence, optimized for input protein
    """
    start_time = time.time()
    array = list(protein_seq)
    N = len(array)

    aminoacids = [0] * N
    for i in range(N):
        aminoacids[i] = array[i]

    # Y matrix whose entry yij  is equal to 1 if ith amino acid in the protein is the j th amino acid in our list.
    Y = np.zeros((N, len(AMINOACIDS)), dtype=np.int)
    for i in range(len(aminoacids)):
        for j in range(len(AMINOACIDS)):
            if aminoacids[i] == AMINOACIDS[j]:
                Y[i, j] = 1

    test2 = np.sum(Y, axis=1)

    # M matrix whose entry mjk  is equal to 1 if jth amino acid can be represented by codon k.
    if len(test2) != N:
        print('warning! there exists undefined AA letter abbreviation')

    # M matrix whose entry mjk  is equal to 1 if jth amino acid can be represented by codon k.
    M = np.zeros((len(AMINOACIDS), len(CODONS)), dtype=np.int)
    M[0, 0:3] = [1, 1, 1]
    M[1, 3:9] = [1, 1, 1, 1, 1, 1]
    M[2, 9:13] = [1, 1, 1, 1]
    M[3, 13:15] = [1, 1]
    M[4, 15:16] = [1]
    M[5, 16:18] = [1, 1]
    M[6, 18:22] = [1, 1, 1, 1]
    M[7, 22:26] = [1, 1, 1, 1]
    M[8, 26:30] = [1, 1, 1, 1]
    M[9, 30:34] = [1, 1, 1, 1]
    M[10, 34:40] = [1, 1, 1, 1, 1, 1]
    M[11, 40:42] = [1, 1]
    M[12, 42:43] = [1]
    M[13, 43:45] = [1, 1]
    M[14, 45:47] = [1, 1]
    M[15, 47:49] = [1, 1]
    M[16, 49:51] = [1, 1]
    M[17, 51:53] = [1, 1]
    M[18, 53:55] = [1, 1]
    M[19, 55:61] = [1, 1, 1, 1, 1, 1]
    M[20, 61:64] = [1, 1, 1]

    # R matrix for possible codonset of amino acids in the protein
    R = np.zeros((len(aminoacids), len(CODONS)), dtype=np.int)
    for i in range(len(aminoacids)):
        for j in range(len(AMINOACIDS)):
            if aminoacids[i] == AMINOACIDS[j]:
                R[i] = M[j]
                break

    logFitnessValues = np.zeros(64, dtype=np.double)
    for i in range(len(fitness_values)):
        logFitnessValues[i] = np.log(float(fitness_values[i]))

    NumAApairs = np.zeros((len(AMINOACIDS), len(AMINOACIDS)), dtype=np.int)
    for j in range(len(AMINOACIDS)):
        for k in range(len(AMINOACIDS)):
            for i in range(len(aminoacids) - 1):
                if Y[i, j] == 1 and Y[i + 1, k] == 1:
                    NumAApairs[j, k] = NumAApairs[j, k] + 1

    etapair = np.zeros((len(CODONS), len(CODONS)), dtype=np.int)
    for j in range(len(CODONS)):
        for k in range(len(CODONS)):
            for i in range(len(AMINOACIDS)):
                for m in range(len(AMINOACIDS)):
                    if M[i, j] == 1 and M[m, k] == 1:
                        etapair[j, k] = NumAApairs[i, m]

    NumAminoAcidPairCodoPairPossibility = np.zeros((len(AMINOACIDS), len(AMINOACIDS)), dtype=np.int)
    for i in range(len(AMINOACIDS)):
        for m in range(len(AMINOACIDS)):
            NumAminoAcidPairCodoPairPossibility[i, m] = np.sum(M[i]) * np.sum(M[m])

    # %%###################################Build the Model##################################################################
    m = Model("Codon Optimization")

    # Variable Xik=1 if ith amino acid of the protein is assigned to k th codon
    X = {}
    for i in range(len(aminoacids)):
        for j in range(len(CODONS)):
            if R[i, j] == 1:
                X[i, j] = m.addVar(vtype=GRB.BINARY, name="X%s" % str([i, j]))
    m.update()

    # Constraint that every aminoacid is assigned to exactly one codon
    for i in range(len(aminoacids)):
        xvar = []
        for j in range(len(CODONS)):
            if R[i, j] == 1:
                xvar.append(j)
        m.addConstr(sum(X[i, k] * R[i, k] for k in xvar), GRB.EQUAL, 1)
    m.update()

    # Variable Zijk=1 if codon pair jk is used for amino acids i and i+1
    coef = {}
    Z = {}
    for i in range(len(aminoacids) - 1):
        for j in range(len(CODONS)):
            for k in range(len(CODONS)):
                if (R[i, j] == 1 & R[i + 1, k] == 1):
                    Z[i, j, k] = m.addVar(vtype=GRB.BINARY, name="Z%s" % str([i, j, k]))
                    coef[i, j, k] = cpb[j][k]

    m.update()

    for i in range(len(aminoacids) - 1):
        jvar = []
        kvar = []
        for j in range(len(CODONS)):
            if R[i, j] == 1:
                jvar.append(j)
        for k in range(len(CODONS)):
            if R[i + 1, k] == 1:
                kvar.append(k)
        m.addConstr(sum(Z[i, l, n] for l in jvar for n in kvar), GRB.EQUAL, 1)
        for l in jvar:
            for n in kvar:
                m.addConstr(X[i, l] + X[i + 1, n] >= 2 * Z[i, l, n])
                m.addConstr(X[i, l] + X[i + 1, n] <= Z[i, l, n] + 1)

    m.update()

    minCAI = N * math.log(threshold)

    m.addConstr(sum(X[i, j] * logFitnessValues[j] for (i, j) in X), GRB.GREATER_EQUAL, minCAI)
    m.update()

    # Set objective function
    obj = sum(Z[i, j, k] * coef[i, j, k] for (i, j, k) in Z) / (N - 1)

    m.setObjective(obj, GRB.MAXIMIZE)

    m.optimize()
    if m.status == GRB.Status.OPTIMAL:
        print('\nobjective function value: %g' % m.objVal)

    objValue = m.objVal

    for v in m.getVars():
        if v.x >= 1:
            print('%s %g' % (v.varName, v.x))

    codons = [j for (i, j) in X if X[i, j].X == 1]
    ans = ""
    for j in range(len(codons)):
        ans += CODONS[codons[j]]
    elapsed_time = time.time() - start_time
    print(elapsed_time)
    print("END")
    return f"Objective function value: {objValue}\n" + ans + "\n"
