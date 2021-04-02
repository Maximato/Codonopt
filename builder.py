from os.path import join

from config import DB_CODON_FREQUENCY_START, DB_CODON_PAIR_FREQUENCY_START, CODON_DB_PATH, BICODON_DB_PATH

from code.extractor import extract_db_frequencies, extract_builder_data
from code.calculator import calculate_fitness_values, calculate_cps
from code.writer import create_directory, write_ordering_line, write_ordering_matrix


taxid, organism = extract_builder_data("builder_input.txt")

# build fitness values
codon_frequencies = extract_db_frequencies(taxid, CODON_DB_PATH, DB_CODON_FREQUENCY_START)
fitness_values = calculate_fitness_values(codon_frequencies)

create_directory(join("db", organism))
write_ordering_line(join("db", organism, "fv.txt"), fitness_values)


# build Codon Pair Score (CPS) table
codon_pairs_frequencies = extract_db_frequencies(taxid, BICODON_DB_PATH, DB_CODON_PAIR_FREQUENCY_START)
cps = calculate_cps(codon_frequencies, codon_pairs_frequencies)

write_ordering_matrix(join("db", organism, "cps.txt"), cps)
