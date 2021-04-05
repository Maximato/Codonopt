from os import path, mkdir

from config import CODONS


def create_directory(path_dir: str):
    if path.isdir(path_dir):
        return
    mkdir(path_dir)


def write_ordering_line(filename: str, data_1d: dict):
    """
    Write data (codon frequencies, fitness values) to file in specific order.
    :param filename: filename to write
    :param data_1d: 1 dimensional data
    """
    codons_list = []
    for codon in CODONS:
        codons_list.append(str(data_1d[codon]))

    with open(filename, "w") as f:
        f.write(", ".join(codons_list))


def write_ordering_matrix(filename: str, data_2d: dict):
    """
    Write data (codon pairs frequencies, codon pairs bias) to file in specific order.
    :param filename: filename to write
    :param data_2d: 1 dimensional data
    """
    with open(filename, "w") as f:
        for key1 in CODONS:
            bufer = ""
            for key2 in CODONS:
                bufer += f"{round(data_2d[key1][key2], 3)}, "
            f.write(bufer[:-2] + "\n")
