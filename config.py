import logging
from os.path import join

logging.basicConfig(filename='.log', filemode='w', level=logging.ERROR)

DB_DIR = "db"
CODON_DB_PATH = join(DB_DIR, "codon_refseq.tsv")
BICODON_DB_PATH = join(DB_DIR, "bicodon_refseq.tsv")
DB_COLUMN_DELIMITER = "\t"
DB_CODON_FREQUENCY_START = 12
DB_CODON_PAIR_FREQUENCY_START = 8


"""
20 aminoacids and one stop symbol
aminoacids=['Isoleucine', 'Leucine','Valine','Phenylalanine','Methionine','Cysteine','Alanine','Glycine','Proline',
            'Threonine','Serine','Tyrosine','Tryptophan','Glutamine','Asparagine','Histidine','Glutamic acid',
            'Aspartic acid','Lysine','Arginine','Stop codons'];
"""

# 20 aminoacids and one stop symbol
AMINOACIDS = ['I', 'L', 'V', 'F', 'M', 'C', 'A', 'G', 'P', 'T', 'S', 'Y', 'W', 'Q', 'N', 'H', 'E', 'D', 'K', 'R', 'X']

# 64 codons
CODONS = ['ATT', 'ATA', 'ATC', 'CTA', 'CTC', 'CTG', 'CTT', 'TTA', 'TTG', 'GTT', 'GTA', 'GTC', 'GTG', 'TTT', 'TTC',
          'ATG', 'TGT', 'TGC', 'GCA', 'GCC', 'GCG', 'GCT', 'GGT', 'GGC', 'GGA', 'GGG', 'CCT', 'CCC', 'CCA', 'CCG',
          'ACT', 'ACC', 'ACA', 'ACG', 'TCT', 'TCC', 'TCA', 'TCG', 'AGT', 'AGC', 'TAT', 'TAC', 'TGG', 'CAA', 'CAG',
          'AAT', 'AAC', 'CAT', 'CAC', 'GAA', 'GAG', 'GAT', 'GAC', 'AAA', 'AAG', 'CGT', 'CGC', 'CGA', 'CGG', 'AGA',
          'AGG', 'TAA', 'TAG', 'TGA']

CODON2AA = {
    'ATA': 'I', 'ATC': 'I', 'ATT': 'I', 'ATG': 'M',
    'ACA': 'T', 'ACC': 'T', 'ACG': 'T', 'ACT': 'T',
    'AAC': 'N', 'AAT': 'N', 'AAA': 'K', 'AAG': 'K',
    'AGC': 'S', 'AGT': 'S', 'AGA': 'R', 'AGG': 'R',
    'CTA': 'L', 'CTC': 'L', 'CTG': 'L', 'CTT': 'L',
    'CCA': 'P', 'CCC': 'P', 'CCG': 'P', 'CCT': 'P',
    'CAC': 'H', 'CAT': 'H', 'CAA': 'Q', 'CAG': 'Q',
    'CGA': 'R', 'CGC': 'R', 'CGG': 'R', 'CGT': 'R',
    'GTA': 'V', 'GTC': 'V', 'GTG': 'V', 'GTT': 'V',
    'GCA': 'A', 'GCC': 'A', 'GCG': 'A', 'GCT': 'A',
    'GAC': 'D', 'GAT': 'D', 'GAA': 'E', 'GAG': 'E',
    'GGA': 'G', 'GGC': 'G', 'GGG': 'G', 'GGT': 'G',
    'TCA': 'S', 'TCC': 'S', 'TCG': 'S', 'TCT': 'S',
    'TTC': 'F', 'TTT': 'F', 'TTA': 'L', 'TTG': 'L',
    'TAC': 'Y', 'TAT': 'Y', 'TAA': 'X', 'TAG': 'X',
    'TGC': 'C', 'TGT': 'C', 'TGA': 'X', 'TGG': 'W'
}

AA2CODON = {'I': ['ATA', 'ATC', 'ATT'], 'M': ['ATG'], 'T': ['ACA', 'ACC', 'ACG', 'ACT'], 'N': ['AAC', 'AAT'],
            'K': ['AAA', 'AAG'], 'S': ['AGC', 'AGT', 'TCA', 'TCC', 'TCG', 'TCT'],
            'R': ['AGA', 'AGG', 'CGA', 'CGC', 'CGG', 'CGT'], 'L': ['CTA', 'CTC', 'CTG', 'CTT', 'TTA', 'TTG'],
            'P': ['CCA', 'CCC', 'CCG', 'CCT'], 'H': ['CAC', 'CAT'], 'Q': ['CAA', 'CAG'],
            'V': ['GTA', 'GTC', 'GTG', 'GTT'], 'A': ['GCA', 'GCC', 'GCG', 'GCT'], 'D': ['GAC', 'GAT'],
            'E': ['GAA', 'GAG'], 'G': ['GGA', 'GGC', 'GGG', 'GGT'], 'F': ['TTC', 'TTT'], 'Y': ['TAC', 'TAT'],
            'X': ['TAA', 'TAG', 'TGA'], 'C': ['TGC', 'TGT'], 'W': ['TGG']}
