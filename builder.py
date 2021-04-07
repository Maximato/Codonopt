from os.path import join

from code.calculator import calculate_fitness_values, calculate_cps, calculate_observed_codon_frequencies, calculate_observed_pair_frequencies
from code.extractor import extract_db_frequencies, extract_builder_data
from code.writer import create_directory, write_ordering_line, write_ordering_matrix
from config import DB_CODON_FREQUENCY_START, DB_CODON_PAIR_FREQUENCY_START, CODON_DB_PATH, BICODON_DB_PATH, DB_DIR

taxid, organism = extract_builder_data("builder_input.txt")

# build fitness values and observed codon frequencies
codon_frequencies = extract_db_frequencies(taxid, CODON_DB_PATH, DB_CODON_FREQUENCY_START)
fitness_values = calculate_fitness_values(codon_frequencies)
observed_codon_frequencies = calculate_observed_codon_frequencies(codon_frequencies)

create_directory(join(DB_DIR, organism))
write_ordering_line(join(DB_DIR, organism, "fv.txt"), fitness_values)
write_ordering_line(join(DB_DIR, organism, "ocf.txt"), observed_codon_frequencies)


# build Codon Pair Score (CPS) table and observed codon pair frequencies table
codon_pairs_frequencies = extract_db_frequencies(taxid, BICODON_DB_PATH, DB_CODON_PAIR_FREQUENCY_START)
cps = calculate_cps(codon_frequencies, codon_pairs_frequencies)
observed_pair_frequencies = calculate_observed_pair_frequencies(codon_pairs_frequencies)

write_ordering_matrix(join(DB_DIR, organism, "cps.txt"), cps)
write_ordering_matrix(join(DB_DIR, organism, "opf.txt"), observed_pair_frequencies)
