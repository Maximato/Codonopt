import math
from config import CODONS, CODON2AA, AA2CODON, AMINOACIDS


def calculate_fitness_values(codon_frequencies: dict) -> dict:
    """
    Calculate fitness values

    :param codon_frequencies: codon frequencies dictionary
    :return: fitness values dictionary
    """
    max_frequencies = {}
    for a in AMINOACIDS:
        max_frequencies[a] = max(list(map(lambda c: codon_frequencies[c], AA2CODON[a])))

    fitness_values = {}
    for codon, frequency in codon_frequencies.items():
        aa = CODON2AA[codon]
        fitness_values[codon] = round(frequency / max_frequencies[aa], 2)
    return fitness_values


def calculate_cps(codon_frequencies: dict, codon_pair_frequencies: dict) -> dict:
    """
    Calculate Codon Pair Score (CPS) table that needed for calculating Codon Pair Bias (CPB) index.

    :param codon_frequencies: codon frequencies dictionary
    :param codon_pair_frequencies: codon pair frequencies dictionary
    :return: Codon Pair Score table
    """
    aa_frequences = {}
    for a in AMINOACIDS:
        aa_frequences[a] = sum(list(map(lambda x: codon_frequencies[x], AA2CODON[a])))

    aa_pairs_frequences = _calc_aa_pairs_frequencies(codon_pair_frequencies)

    # zero value lead to error in logarithm calculations
    _remove_zeros(codon_pair_frequencies)

    cps = {c: {} for c in CODONS}
    for c1 in CODONS:
        for c2 in CODONS:
            aa_pair = CODON2AA[c1] + CODON2AA[c2]
            cps[c1][c2] = math.log(codon_pair_frequencies[c1+c2] * aa_frequences[CODON2AA[c1]] * aa_frequences[CODON2AA[c2]] /
                                   codon_frequencies[c1] / codon_frequencies[c2] / aa_pairs_frequences[aa_pair])
    return cps


def _calc_aa_pairs_frequencies(codon_pair_frequencies: dict) -> dict:
    """
    Calculate amino acid pairs frequencies based on codon pairs frequencies.

    :param codon_pair_frequencies: codon pair frequencies dictionary
    :return: amino acid pair frequencies dictionary
    """
    aa_pair_frequencies = {a: {} for a in AMINOACIDS}
    for a1 in AMINOACIDS:
        for a2 in AMINOACIDS:
            aa_pair_frequency = 0
            for c1 in AA2CODON[a1]:
                for c2 in AA2CODON[a2]:
                    aa_pair_frequency += codon_pair_frequencies[c1 + c2]
            aa_pair_frequencies[a1+a2] = aa_pair_frequency
    return aa_pair_frequencies


def _remove_zeros(codon_pair_frequencies: dict) -> dict:
    """
    Remove zeros from codon pair frequencies
    :param codon_pair_frequencies: codon pair frequencies
    :return: new codon pair frequencies
    """
    for key, value in codon_pair_frequencies.items():
        if value == 0:
            print(f"Warning! Zero codon pair frequency for {key}. Changed to 0.001")
            codon_pair_frequencies[key] = 0.001
    return codon_pair_frequencies
