from os.path import join

from config import DB_COLUMN_DELIMITER, DB_DIR


def extract_db_frequencies(taxid: str or int, db_path: str, frequency_col_from: int) -> dict:
    """
    Extracting frequencies (codon or codon pair) from database
    :param taxid: taxid ID of organism
    :param db_path: path to database with frequencies information
    :param frequency_col_from: column id for first frequencies
    :return: frequencies dictionary
    """
    taxid = str(taxid)
    with open(db_path, "r") as f:
        header = f.readline().rstrip().split(DB_COLUMN_DELIMITER)
        key_order = list(map(lambda c: c.upper(), header[frequency_col_from:]))
        frequencies = {codon: 0 for codon in key_order}
        for line in f:
            row = line.split(DB_COLUMN_DELIMITER)
            db_taxid = row[2]
            if db_taxid == taxid:
                for codon, frequency in zip(key_order, row[frequency_col_from:]):
                    frequencies[codon] += int(frequency)
    return frequencies


def extract_codonopt_data(filename) -> (str, float, list):
    """
    Extracting information for running optimization

    :param filename: path to filename
    :return: organism, threshold, proteins sequences
    """
    with open(filename, "r") as r:
        organism = r.readline().split()[1]
        threshold = float(r.readline().split()[1])

        seqs = []
        for line in r:
            seqs.append(line.rstrip())

    fitness_values, cps = _extract_builded_data(DB_DIR, organism)
    return fitness_values, cps, threshold, seqs


def extract_builder_data(filename) -> (str, str, str):
    """
    Extracting information for building db for new organism

    :param filename: path to filename
    :return: taxid id, new organism name
    """
    with open(filename, "r") as r:
        taxid = r.readline().rstrip().split()[1]
        organism = r.readline().rstrip().split()[1]
    return taxid, organism


def _extract_builded_data(db_path: str, organism: str) -> (float, float):
    """
    Extracting information from db for optimization

    :param db_path: path do db
    :param organism: organism
    :return: fitness values, codon pair bias
    """
    file_cpb = join(db_path, organism, "cps.txt")
    cps_str = open(file_cpb).readlines()
    cps = []
    for i, line in enumerate(cps_str):
        cps[i] = line.split(",")

    file_fv = join(db_path, organism, "fv.txt")
    with open(file_fv) as f:
        fitness_values = f.read().split(',')

    return fitness_values, cps
