import math
import time

from config import *


try:
    import numpy as np
    from gurobipy import *

except Exception as e:
    logging.exception("External required libraries importing exception occurred")
    raise e


def max_cpb_st_cai_optimization(protein_seq: str, fitness_values: list, cps: list, threshold: float) -> str:
    """
    Take amino acid sequence of protein as input and mathematically optimize DNA sequence. Fully based on MaxCPBstCAI
    function from software implementation of approach described in "Codon Optimization: A Mathematical Programming Approach"
    article. The function maximizes Codon Pair Bias (CPB) index (that depends of the occurrence of codon pairs)
    when the CAI (Codon Adaptation Index) does not fall below the specified value (threshold).

    :param protein_seq: sequence of protein for optimization
    :param fitness_values: fitness values list
    :param cps: Codon Pair Score (CPS) table
    :param threshold: the min value for CAI
    :return: string with DNA sequence, optimized for input protein
    """
    start_time = time.time()
    aminoacids = list(protein_seq)
    N = len(aminoacids)

    Y = _create_Y(aminoacids)

    # M matrix whose entry mjk  is equal to 1 if jth amino acid can be represented by codon k.
    test2 = np.sum(Y, axis=1)
    if len(test2) != N:
        print('warning! there exists undefined AA letter abbreviation')

    M = _create_M()
    R = _create_R(aminoacids, M)

    model, X = _create_model(aminoacids, R)

    # Variable Zijk=1 if codon pair jk is used for amino acids i and i+1
    coef = {}
    Z = {}
    for i in range(len(aminoacids) - 1):
        for j in range(len(CODONS)):
            for k in range(len(CODONS)):
                if R[i, j] == 1 & R[i + 1, k] == 1:
                    Z[i, j, k] = model.addVar(vtype=GRB.BINARY, name="Z%s" % str([i, j, k]))
                    coef[i, j, k] = cps[j][k]
    model.update()

    for i in range(len(aminoacids) - 1):
        jvar = []
        kvar = []
        for j in range(len(CODONS)):
            if R[i, j] == 1:
                jvar.append(j)
        for k in range(len(CODONS)):
            if R[i + 1, k] == 1:
                kvar.append(k)
        model.addConstr(sum(Z[i, l, n] for l in jvar for n in kvar), GRB.EQUAL, 1)
        for l in jvar:
            for n in kvar:
                model.addConstr(X[i, l] + X[i + 1, n] >= 2 * Z[i, l, n])
                model.addConstr(X[i, l] + X[i + 1, n] <= Z[i, l, n] + 1)
    model.update()

    minCAI = N * math.log(threshold)
    logFitnessValues = np.zeros(64, dtype=np.double)
    for i in range(len(fitness_values)):
        logFitnessValues[i] = np.log(fitness_values[i])
    model.addConstr(sum(X[i, j] * logFitnessValues[j] for (i, j) in X), GRB.GREATER_EQUAL, minCAI)
    model.update()

    # Set objective function
    obj = sum(Z[i, j, k] * coef[i, j, k] for (i, j, k) in Z) / (N - 1)
    model.setObjective(obj, GRB.MAXIMIZE)
    model.optimize()

    objValue = model.objVal

    codons = [j for (i, j) in X if X[i, j].X == 1]
    ans = ""
    for j in range(len(codons)):
        ans += CODONS[codons[j]]
    elapsed_time = time.time() - start_time
    print(elapsed_time)
    print("END")
    return f"Objective function value: {objValue}\n" + ans + "\n"


def min_rcpb_st_rcb_optimization(protein_seq: str, freq_codons: list, freq_codon_pair: list, threshold: float):
    """
    Take amino acid sequence of protein as input and mathematically optimize DNA sequence. Fully based on MinRCPBstRCB
    function from software implementation of approach described in "Codon Optimization: A Mathematical Programming Approach"
    article. The function minimize Relative Codon Pair Bias (RCPB) index when the Relative Codon Bias (RCB) index
    does not rise above the specified value (threshold).

    :param protein_seq: sequence of protein for optimization
    :param freq_codons: observed frequency of each codon ([0.421418,	0.538557,	0.578582, ... ])
    :param freq_codon_pair: observed frequency of codon pair
    :param threshold: the max value for RCB
    :return: string with DNA sequence, optimized for input protein
    """

    start_time = time.time()

    aminoacids = list(protein_seq)
    N = len(aminoacids)

    Y = _create_Y(aminoacids)

    # M matrix whose entry mjk  is equal to 1 if jth amino acid can be represented by codon k.
    test2 = np.sum(Y, axis=1)
    if len(test2) != N:
        print('warning! there exists undefined AA letter abbreviation')

    M = _create_M()
    R = _create_R(aminoacids, M)

    model, X = _create_model(aminoacids, R)

    Z = {}
    for i in range(len(aminoacids) - 1):
        for j in range(len(CODONS)):
            for k in range(len(CODONS)):
                if R[i, j] == 1 and R[i + 1, k] == 1:
                    Z[i, j, k] = model.addVar(vtype=GRB.BINARY, name="Z%s" % str([i, j, k]))
    model.update()

    for i in range(len(aminoacids) - 1):
        jvar = []
        kvar = []
        for j in range(len(CODONS)):
            if R[i, j] == 1:
                jvar.append(j)
        for k in range(len(CODONS)):
            if R[i + 1, k] == 1:
                kvar.append(k)
        model.addConstr(sum(Z[i, l, n] for l in jvar for n in kvar), GRB.EQUAL, 1)
        for l in jvar:
            for m in kvar:
                model.addConstr(X[i, l] + X[i + 1, m] >= 2 * Z[i, l, m])
                model.addConstr(X[i, l] + X[i + 1, m] <= Z[i, l, m] + 1)
    model.update()

    codondev = {}
    for i in range(len(CODONS)):
        codondev[i] = model.addVar(vtype=GRB.CONTINUOUS, name="codondev%s" % str([i]))
    model.update()

    eta = np.sum(Y, axis=0)
    for j in range(len(CODONS)):
        xvar = []
        etatemp = []
        for k in range(len(AMINOACIDS)):
            if M[k, j] == 1:
                etatemp = eta[k]
        for i in range(len(aminoacids)):
            if R[i, j] == 1:
                xvar.append(i)
        if etatemp != 0:
            model.addConstr(100 * sum(X[i, j] for i in xvar) / etatemp, GRB.LESS_EQUAL,
                        100 * (freq_codons[j] + codondev[j]))
            model.addConstr(100 * sum(X[i, j] for i in xvar) / etatemp, GRB.GREATER_EQUAL,
                        100 * (freq_codons[j] - codondev[j]))
    model.update()

    codonpairdev = {}
    for j in range(len(CODONS)):
        for k in range(len(CODONS)):
            codonpairdev[j, k] = model.addVar(vtype=GRB.CONTINUOUS, name="codonpairdev%s" % str([j, k]))
    model.update()

    for j in range(len(CODONS)):
        for k in range(len(CODONS)):
            ivar = []
            for i in range(len(aminoacids) - 1):
                if R[i, j] == 1 and R[i + 1, k] == 1:
                    ivar.append(i)
            if len(ivar) != 0:
                etapairtemp = len(ivar)
                model.addConstr(100 * sum(Z[i, j, k] for i in ivar) / etapairtemp, GRB.LESS_EQUAL,
                            100 * (freq_codon_pair[j][k] + codonpairdev[j, k]))
                model.addConstr(100 * sum(Z[i, j, k] for i in ivar) / etapairtemp, GRB.GREATER_EQUAL,
                            100 * (freq_codon_pair[j][k] - codonpairdev[j, k]))
    model.update()

    AAdev = {}
    for i in range(len(AMINOACIDS)):
        AAdev[i] = model.addVar(vtype=GRB.CONTINUOUS, name="AAdev%s" % str([i]))
    model.update()

    NumAminoAcidCodonPossibility = (np.sum(M, axis=1))
    for i in range(len(AMINOACIDS)):
        jvar = []
        for j in range(len(CODONS)):
            if M[i, j] == 1:
                jvar.append(j)
        model.addConstr(100 * sum(codondev[j] for j in jvar) / NumAminoAcidCodonPossibility[i], GRB.EQUAL, 100 * AAdev[i])
    model.update()

    AApairdev = {}
    for i in range(len(AMINOACIDS)):
        for j in range(len(AMINOACIDS)):
            AApairdev[i, j] = model.addVar(vtype=GRB.CONTINUOUS, name="AApairdev%s" % str([i, j]))
    model.update()

    NumAminoAcidPairCodoPairPossibility = np.zeros((len(AMINOACIDS), len(AMINOACIDS)), dtype=np.int)
    for i in range(len(AMINOACIDS)):
        for m in range(len(AMINOACIDS)):
            NumAminoAcidPairCodoPairPossibility[i, m] = np.sum(M[i]) * np.sum(M[m])

    for i in range(len(AMINOACIDS)):
        for j in range(len(AMINOACIDS)):
            kvar = []
            lvar = []
            for k in range(len(CODONS)):
                if M[i, k] == 1:
                    kvar.append(k)
            for l in range(len(CODONS)):
                if M[j, l] == 1:
                    lvar.append(l)
            model.addConstr(
                100 * sum(codonpairdev[k, l] for k in kvar for l in lvar) / NumAminoAcidPairCodoPairPossibility[i, j],
                GRB.EQUAL, 100 * AApairdev[i, j])
    model.update()

    NumAA = np.sum(Y, axis=0)
    maxRCB = threshold * N
    model.addConstr(100 * sum(AAdev[j] * NumAA[j] for (j) in AAdev), GRB.LESS_EQUAL, 100 * maxRCB)
    model.update()

    NumAApairs = np.zeros((len(AMINOACIDS), len(AMINOACIDS)), dtype=np.int)
    for j in range(len(AMINOACIDS)):
        for k in range(len(AMINOACIDS)):
            for i in range(len(aminoacids) - 1):
                if Y[i, j] == 1 and Y[i + 1, k] == 1:
                    NumAApairs[j, k] = NumAApairs[j, k] + 1

    obj = 100 * sum(AApairdev[i, j] * NumAApairs[i, j] for (i, j) in AApairdev)
    model.setObjective(obj, GRB.MINIMIZE)
    model.optimize()

    objValue = model.objVal / (100 * (N - 1))

    # Write codons into the file
    codons = [(j) for (i, j) in X if X[i, j].X == 1]

    ans = ""
    for j in range(len(codons)):
        ans += CODONS[codons[j]]
    elapsed_time = time.time() - start_time
    print(elapsed_time)
    print("END")
    return f"Objective function value: {objValue}\n" + ans + "\n"


def _create_Y(aminoacids):
    # Y matrix whose entry yij  is equal to 1 if ith amino acid in the protein is the j th amino acid in our list.
    N = len(aminoacids)
    Y = np.zeros((N, len(AMINOACIDS)), dtype=np.int)
    for i in range(len(aminoacids)):
        for j in range(len(AMINOACIDS)):
            if aminoacids[i] == AMINOACIDS[j]:
                Y[i, j] = 1
    return Y


def _create_M():
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
    return M


def _create_R(aminoacids, M):
    # R matrix for possible codonset of amino acids in the protein
    R = np.zeros((len(aminoacids), len(CODONS)), dtype=np.int)
    for i in range(len(aminoacids)):
        for j in range(len(AMINOACIDS)):
            if aminoacids[i] == AMINOACIDS[j]:
                R[i] = M[j]
                break
    return R


def _create_model(aminoacids, R):
    # Build the Model
    model = Model("Codon Optimization")

    # Variable Xik=1 if ith amino acid of the protein is assigned to k th codon
    X = {}
    for i in range(len(aminoacids)):
        for j in range(len(CODONS)):
            if R[i, j] == 1:
                X[i, j] = model.addVar(vtype=GRB.BINARY, name="X%s" % str([i, j]))
    model.update()

    # Constraint that every aminoacid is assigned to exactly one codon
    for i in range(len(aminoacids)):
        xvar = []
        for j in range(len(CODONS)):
            if R[i, j] == 1:
                xvar.append(j)
        model.addConstr(sum(X[i, k] * R[i, k] for k in xvar), GRB.EQUAL, 1)
    model.update()
    return model, X
