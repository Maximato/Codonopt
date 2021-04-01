import math
from config import CODONS, CODON2AA, AA2CODON


def calculate_fitness_values(codon_frequencies: dict) -> dict:
    """
    Calculate fitness values

    :param codon_frequencies: codon frequencies dictionary
    :return: fitness values dictionary
    """
    max_frequencies = {}
    for a in AA2CODON:
        max_frequencies[a] = max(list(map(lambda c: codon_frequencies[c], AA2CODON[a])))

    fitness_values = {}
    for codon, frequency in codon_frequencies.items():
        fitness_values[codon] = round(frequency / max_frequencies[CODON2AA[codon]], 2)
    return fitness_values


def calculate_cps(codon_frequencies: dict, codon_pair_frequencies: dict) -> dict:
    """
    Calculate codon pair bias

    :param codon_frequencies: codon frequencies dictionary
    :param codon_pair_frequencies: codon pair frequencies dictionary
    :return: codon pair bias table
    """
    aa_frequences = {}
    for a in AA2CODON:
        aa_frequences[a] = sum(list(map(lambda x: codon_frequencies[x], AA2CODON[a])))

    aa_pairs_frequences = _calc_aa_pair_frequencies(codon_pair_frequencies)

    cps = {c: {} for c in CODONS}
    for c1 in CODONS:
        for c2 in CODONS:
            aa_pair = CODON2AA[c1] + CODON2AA[c2]
            cps[c1][c2] = math.log(codon_pair_frequencies[c1+c2] * aa_frequences[CODON2AA[c1]] * aa_frequences[CODON2AA[c2]] /
                                   codon_frequencies[c1] / codon_frequencies[c2] / aa_pairs_frequences[aa_pair])
    return cps


def normalize_frequencies(frequencies: dict) -> dict:
    sum_frequencies = sum(frequencies.values())
    for key in frequencies:
        frequencies[key] = frequencies[key]/sum_frequencies*124
    return frequencies


def _calc_aa_pair_frequencies(codon_pairs_frequencies):
    aa_pair_frequencies = {a: {} for a in AA2CODON}
    for a1 in AA2CODON:
        for a2 in AA2CODON:
            aa_pair_frequence = 0
            for c1 in AA2CODON[a1]:
                for c2 in AA2CODON[a2]:
                    aa_pair_frequence += codon_pairs_frequencies[c1+c2]
            aa_pair_frequencies[a1+a2] = aa_pair_frequence
    return aa_pair_frequencies
