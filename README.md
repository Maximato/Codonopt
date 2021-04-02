## About

Codonopt program is designed for codon optimization of protein sequences.  
The core of it based on the mathematical codon optimization proposed  
by  Alper Şen, Kamyar Kargar, Esma Akgün and Mustafa Çelebi Pınar in ["Codon optimization: A
mathematical programming approach"](http://alpersen.bilkent.edu.tr/codonoptimization/CodonPaper.pdf) article.
Codonopt is software that extend MAXCPBstCAI function by providing optimization of specified set of sequences.  
Moreover, now it is available to build fitness values and Codon Pair Score (CPS) tables for your own organism  
or from organism in database. To do this you need to create codon usage database for your organism  
or download prepared codon usage database from [HIVE platform](https://hive.biochemistry.gwu.edu/cuts/about)  
that contains a lot of organism.


### Installation

To run the program Python version 3.8 (or higher) and Gurobi software are required.
Two python packages is also needed: numpy (mathematical library) and gurobipy (Gurobi Optimization library).

They can be installed with pip manager by commands:

        pip install numpy    
        python -m pip install -i https://pypi.gurobi.com gurobipy



### codonopt.py script

Read protein sequence (or set of protein sequences) from **codonopt_input.txt** and mathematically optimize for specific organism.
It is fully based on `MAXCPBstCAI` function described in "Codon Optimization: A Mathematical Programming Approach" article.
The function maximizes Codon Pair Bias (CPB) index (that depends of the occurrence of codon pairs)  
when the CAI (Codon Adaptation Index) does not fall below the specified value (threshold).

To run this script you need to specify some values in **codonopt_input.txt** file:
- organism --- organism for optimization. This value should match with directory name in database. That directory should contain with fitness values and CPS table (details in `builder.py` script information)
- minCAI --- threshold for CAI
- sequences --- protein sequences

Example:

        organism: ecoli
        minCAI: 0.8
        PLKATSTPVSIKSTLLGGGSATVKFKYKGEELEVDISK
        LNIEDEHRLHETSKEPDVSLGSTWLSDFPQAWAETGGMGLAVRQAPLIIPLKATS


Run **codonopt.py** script from command line. All results will be saved in **output.txt** file.


### builder.py script

Builder script is designed to build fitness values and Codon Pair Score (CPS) tables for a specific taxid.  
These values are calculated based on the Codon and Codon-Pair Usage Tables stored in the database  
which must be created and configured to run the script.  
The calculated fitness values and codon tables are saved in a separate  
directory and can then be used for optimization protein sequences.  
To do this, you must specify the name of the directory with data as the value **organism** in the **codonopt_input.txt** file.

#### Running

First of all you need to created and configured Codon and Codon-Pair Usage Tables files.
This files can be downloaded from [HIVE platform](https://hive.biochemistry.gwu.edu/cuts/about) or created manually.
It must contain taxid column and columns named as codon/codon-pair (see example).

Codon usage example:

        Taxid	TTT	TTC	TTA	TTG ...
        2161	28430	10752	31143	7667 ...
        1204725	15851	13702	18373	6926 ...

Codon pair usage example:

        Taxid	tttttt	tttttc	ttttta	tttttg ...
        2161	4450	50672	37181	9284 ...
        1204725	15631	40625	31317	1566 ...

Then you need configure the files by specifying constants in **config.py** file:
- `DB_DIR` --- root directory for database ("db" default)
- `CODON_DB` --- name codon usage file ("codon_db.tsv" default)
- `BICODON_DB` --- name codon pair usage file ("bicodon_db.tsv" default)
- `DB_COLUMN_DELIMITER` --- columns delimiter in codon/codon-pair usage files("\t" default)
- `DB_TAXID_INDEX` --- index of taxid column (0 default)
- `DB_CODON_FREQUENCY_START` --- index of first column with frequencies in codon usage file (2 default)
- `DB_CODON_PAIR_FREQUENCY_START` --- index of first column with frequencies in codon pair usage file (2 default)

If there are some troubles with database creating you can use short  
example Codon and Codon-Pair Usage Tables with some popular organisms.  
This files come with Codonopt program (**codon_db.tsv** and **bicodon_db.tsv**) just use default configuration parameters.

Besides, to run this script you need to specify some values in **builder_input.txt** file:

- taxid --- taxid ID of organism that stored in the database.
- organism --- name of directory to save fitness values and Codon Pair Score (CPS) table. This should be used as **organism** parameter to run optimization with built data.

Now you can run `builder.py` script from command line and check new built data in `DB_DIR` directory.

#### Built data

The built data are simple text files (**fv.txt** and **cps.txt**) that store fitness values and Codon Pair Score (CPS) table. Each value are addressed corresponding to `CODONS` constant from configuration file.

File **fv.txt** contains Fitness Values of each codon (1 by 64) that separated by comma. Small 1 by 3 example:

        1, 0.795767933108405, 0.229059703681972

File **cps.txt** contains Codon Pair Score (CPS) of each codon pair in table (64 by 64). Columns of that separated by comma. Small 3 by 3 example:

        -0.345,	-0.221,  0.125
        0.341,	 0.204,	-0.011
        0.272,	 0.225,	 0.364

