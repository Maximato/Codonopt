## About

Codonopt program is designed for codon optimization of protein sequences. The core of it based on the mathematical codon optimization proposed by Alper Şen, Kamyar Kargar, Esma Akgün and Mustafa Çelebi Pınar in [Codon optimization: A mathematical programming approach](http://alpersen.bilkent.edu.tr/codonoptimization/CodonPaper.pdf) article. Codonopt is software that extend `MaxCPBstCAI` and `MinRCPBstRCB` functions by providing optimization of specified set of sequences. Moreover, now it is available to build fitness values, Codon Pair Score (CPS) table and observed codon/codon-pair frequencies for your own organism or for organism from database. To do this you need to create codon usage database for your organism or download prepared codon usage database from [HIVE platform](https://hive.biochemistry.gwu.edu/cuts/about) that contains a lot of organism.

### Installation

To run the program Python version 3.8 (or higher) and Gurobi software are required. Two python packages is also needed: numpy (mathematical library) and gurobipy (Gurobi Optimization library). They can be installed with pip manager by commands:

        pip install numpy    
        python -m pip install -i https://pypi.gurobi.com gurobipy

### codonopt.py script

Read protein sequence (or set of protein sequences) from **codonopt_input.txt** and mathematically optimize for specific organism. There are two options:
- based on `MaxCPBstCAI` function that maximizes Codon Pair Bias (CPB) index (that depends of the occurrence of codon pairs) when the CAI (Codon Adaptation Index) does not fall below the specified value (threshold)
- based on `MinRCPBstRCB` function that minimize Relative Codon Pair Bias (RCPB) index when the Relative Codon Bias (RCB) index does not rise above the specified value (threshold)
Both functions are described in "Codon Optimization: A Mathematical Programming Approach" article.

To run this script you need to specify some values in **codonopt_input.txt** file:
- organism --- organism for optimization. This value should match with directory name in database. That directory should contain fitness values and CPS table (**fv.txt** and **cps.txt** files, details in **builder.py** script information)
- method --- method for optimization (`MaxCPBstCAI` or `MinRCPBstRCB`)
- threshold --- threshold, depends on the method (min value of CAI for `MaxCPBstCAI` or max value of RCB for `MinRCPBstRCB`)
- sequences --- protein sequences (on a new line each)

There are several available organisms that was previously prepared (stored in **db** directory) and can be used for optimization (bacillus_anthracis, corynebacterium_diphtheriae, escherichia_coli, lactococcus_lactis, pseudomonas_syringae, staphylococcus_aureus, streptococcus_pneumoniae). If you want to run optimization for other organism you should build fitness values, Codon Pair Score table and observed codon/codon-pair frequencies before (**builder.py** script).

Example:

        organism: escherichia_coli
        method: MaxCPBstCAI   
        threshold: 0.8
        PLKATSTPVSIKSTLLGGGSATVKFKYKGEELEVDISK
        LNIEDEHRLHETSKEPDVSLGSTWLSDFPQAWAETGGMGLAVRQAPLIIPLKATS

Run **codonopt.py** script from command line. All results will be saved in **output.txt** file.

### builder.py script

Builder script is designed to build fitness values, Codon Pair Score (CPS) table and observed codon/codon-pair frequencies for a specific taxid. These values are calculated based on the Codon and Codon-Pair Usage Tables stored in the database which must be created and configured to run the script. The calculated fitness values, Codon Pair Score (CPS) table and observed codon/codon-pair frequencies are saved in a separate directory and can then be used for optimization protein sequences. To do this, you must specify the name of the directory with calculated values (fitness values, Codon Pair Score (CPS) table and observed codon/codon-pair frequencies) as the value **organism** in the **codonopt_input.txt** file.

#### Running

First of all you need to created and configured Codon and Codon-Pair Usage Tables files. This files can be downloaded from [HIVE platform](https://hive.biochemistry.gwu.edu/cuts/about) or created manually. It must contain taxid column and columns named as codon/codon-pair (see example).

Codon usage example:

        Taxid	TTT	TTC	TTA	TTG ...
        2161	28430	10752	31143	7667 ...
        1204725	15851	13702	18373	6926 ...

Codon pair usage example:

        Taxid	tttttt	tttttc	ttttta	tttttg ...
        2161	4450	50672	37181	9284 ...
        1204725	15631	40625	31317	1566 ...

Then you need to configure the files by specifying constants in **config.py** file:
- `DB_DIR` --- root directory for database ("db" default)
- `CODON_DB` --- name codon usage file ("codon_db.tsv" default)
- `BICODON_DB` --- name codon pair usage file ("bicodon_db.tsv" default)
- `DB_COLUMN_DELIMITER` --- columns delimiter in codon/codon-pair usage files("\t" default)
- `DB_TAXID_INDEX` --- index of taxid column (0 default)
- `DB_CODON_FREQUENCY_START` --- index of first column with frequencies in codon usage file (2 default)
- `DB_CODON_PAIR_FREQUENCY_START` --- index of first column with frequencies in codon pair usage file (2 default)

If there are some troubles with database creating you can use small example Codon and Codon-Pair Usage Tables with some organisms. This files come with Codonopt program (**codon_db.tsv** and **bicodon_db.tsv**) just use default configuration parameters.

Besides, to run **builder.py** script you need to specify some values in **builder_input.txt** file:

- taxid --- taxid ID of organism that stored in the database.
- organism --- name of directory (without spaces) to save fitness values, Codon Pair Score (CPS) table and observed codon/codon-pair frequencies. This should be used as **organism** parameter to run optimization with built data.

Default configuration example:

        taxid: 562
        organism: ecoli_test

Now you can run **builder.py** script from command line and check new built data in `DB_DIR` directory.

#### Built data

The built data are simple text files (**fv.txt**, **cps.txt**, **ocf.txt** and **opf.txt**) that store fitness values, Codon Pair Score (CPS) table and observed codon/codon-pair frequencies. Each value are addressed corresponding to `CODONS` constant from configuration file.

Files **fv.txt** and **ocf.txt** contain Fitness Values and observed codon frequencies of each codon (1 by 64) that separated by comma. Small 1 by 3 example of fitness values:

        1, 0.795767933108405, 0.229059703681972

Files **cps.txt** and **opf.txt** contain Codon Pair Score (CPS) and observed codon-pair frequencies of each codon pair in table (64 by 64). Columns of that separated by comma. Small 3 by 3 example of CPS:

        -0.345,	-0.221,  0.125
        0.341,	 0.204,	-0.011
        0.272,	 0.225,	 0.364
