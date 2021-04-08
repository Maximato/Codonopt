from config import *

from datetime import datetime

from code.extractor import extract_codonopt_data
from code.optimizer import max_cpb_st_cai_optimization, min_rcpb_st_rcb_optimization

try:
    line_data, matrix_data, method, threshold, seqs = extract_codonopt_data("codonopt_input.txt")

    if method == "MaxCPBstCAI":
        run_optimization = max_cpb_st_cai_optimization
    elif method == "MinRCPBstRCB":
        run_optimization = min_rcpb_st_rcb_optimization
    else:
        raise Exception(f"Unknown method {method}")

    with open("output.txt", "w") as w:
        w.write("last update: " + datetime.now().ctime() + "\n\n")
        for seq in seqs:
            optimization_result = run_optimization(seq, line_data, matrix_data, threshold)
            w.write(optimization_result + "\n")

except Exception as e:
    logging.exception("Runtime exception occurred")
    raise e
