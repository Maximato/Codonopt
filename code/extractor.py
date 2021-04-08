from os.path import join

from config import DB_COLUMN_DELIMITER, DB_TAXID_INDEX, DB_DIR


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
            db_taxid = row[DB_TAXID_INDEX]
            if db_taxid == taxid:
                for codon, frequency in zip(key_order, row[frequency_col_from:]):
                    frequencies[codon] += int(frequency)
    return frequencies


def extract_codonopt_data(filename) -> (str, float, list):
    """
    Extracting information for running optimization

    :param filename: path to filename
    :return: fitness values, Codon Pair Score (CPS) table, threshold CAI, sequences
    """
    with open(filename, "r") as r:
        organism = r.readline().split()[1]
        method = r.readline().split()[1]
        threshold = float(r.readline().split()[1])

        seqs = []
        for line in r:
            seqs.append(line.rstrip())

    line_data, matrix_data = _extract_built_data(DB_DIR, organism, method)
    return line_data, matrix_data, method, threshold, seqs


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


def _extract_built_data(db_path: str, organism: str, method: str) -> (list, list):
    """
    Extracting built information from db for optimization

    :param db_path: path do database directory
    :param organism: organism (name of directory with built data)
    :param method: method for optimization (MaxCPBstCAI or MinRCPBstRCB)
    :return: line data (fitness values, observed frequencies), matrix data (CPS table, observed codon pair frequencies)
    """

    if method == "MaxCPBstCAI":
        codon_data = "fv.txt"
        codon_pair_data = "cps.txt"
    elif method == "MinRCPBstRCB":
        codon_data = "ocf.txt"
        codon_pair_data = "opf.txt"
    else:
        raise Exception(f"Unknown method {method}")

    file_codon_data = join(db_path, organism, codon_data)
    with open(file_codon_data) as f:
        line_data = f.read().rstrip().split(',')

    file_codon_pair_data = join(db_path, organism, codon_pair_data)
    matrix_data_str = open(file_codon_pair_data).readlines()
    matrix_data = []
    for i, line in enumerate(matrix_data_str):
        matrix_data.append(list(map(float, line.rstrip().split(","))))

    return list(map(float, line_data)), matrix_data
