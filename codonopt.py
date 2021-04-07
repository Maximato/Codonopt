from config import *

from datetime import datetime

from code.extractor import extract_codonopt_data
from code.optimizer import MAXCPBstCAI_optimization, MinRCPBstRCB_optimization

try:
    line_data, matrix_data, method, threshold, seqs = extract_codonopt_data("codonopt_input.txt")

    if method == "MaxCPBstCAI":
        run_optimization = MAXCPBstCAI_optimization
    elif method == "MinRCPBstRCB":
        run_optimization = MinRCPBstRCB_optimization
    else:
        raise Exception("Unknown method")

    with open("output.txt", "w") as w:
        w.write("last update: " + datetime.now().ctime() + "\n\n")
        for seq in seqs:
            optimization_result = run_optimization(seq, line_data, matrix_data, threshold)
            w.write(optimization_result + "\n")

except Exception as e:
    logging.exception("Runtime exception occurred")
    raise e
