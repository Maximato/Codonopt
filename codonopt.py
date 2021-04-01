from config import *

from datetime import datetime

from code.extractor import extract_codonopt_data
from code.optimizer import run_optimization

try:
    fitness_values, cps, threshold, seqs = extract_codonopt_data("codonopt_input.txt")

    with open("output.txt", "w") as w:
        w.write("last update: " + datetime.now().ctime() + "\n\n")
        for seq in seqs:
            w.write(run_optimization(seq, cps, fitness_values, threshold).replace("U", "T") + "\n")

except Exception as e:
    logging.exception("Runtime exception occurred")
    raise e
