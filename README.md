## About

Codonopt program is designed for codon optimization for proteins and consist of two python scripts.

`codonopt.py`
Take protein sequence (or set of protein sequences) as input and mathematically optimize for specific organism.
Fully based on MAXCPBstCAI function from software implementation of approach described
in "Codon Optimization: A Mathematical Programming Approach" article.
The function maximizes CPB (a parameter depending on the occurrence of codon pairs) when the CAI
Codon Adaptation Index) does not fall below the specified value (threshold).

`builder.py`




======== DESCRIPTION =============

This software is a reference implementation of the integer programming formulation 
for the codon optimization proposed by  Alper Şen, Kamyar Kargar, Esma Akgün and Mustafa Çelebi Pınar. The 
details of the formulation is provided in the manuscript titled "Codon optimization: A 
mathematical programming approach". The manuscript is available at:
http://alpersen.bilkent.edu.tr/codonoptimization/CodonPaper.pdf


============ INSTALLATION REQUIREMENTS ===========

To build the software, Python and Gurobi are required. Gurobi Optimization library which is gurobipy, 
numpy, math and time packages have to be installed by a Python package manager. 

установка библиотки guribipy:
python -m pip install -i https://pypi.gurobi.com gurobipy


========== INPUT DATA INFO ===========
В файле input_data.txt необходимо задать:
организм - пока только ecoli
minCAI - минимальное значение CAI
белковая последовательность 1
белковая последовательность 2
белковая последовательность 3
и т.д.


"Please enter the Fitness Value file name:" A text file which includes Fitness Value of each codon (1 by 64) should be addressed. Each value should be separated by "," in the file.
Small 1 by 3 Example:  1, 0.795767933108405, 0.229059703681972	

"Please enter the CPB file name:" A text file which includes Codon pair bias values (64 by 64) should be addressed. Each value should be separated by "," in the file. 
Small 3 by 3 Example:  -0.345,	-0.221,  0.125,	
                        0.341,	 0.204,	-0.011,
	                0.272,	 0.225,	 0.364,	

"Please enter the min value for CAI:" A number stating the desired value for minimum codon adaptation index should be entered. For example "0.55"

We assume the following ordering for amino acids and codons sets:

#20 aminoacids
aminoacidset=['I','L','V','F','M','C','A','G','P','T','S','Y','W','Q','N','H','E','D','K','R','X']


#64 codons
codonset=['AUU','AUA','AUC','CUA','CUC','CUG','CUU','UUA','UUG','GUU','GUA','GUC','GUG','UUU','UUC','AUG','UGU','UGC','GCA','GCC','GCG','GCU','GGU','GGC','GGA','GGG','CCU','CCC','CCA','CCG','ACU','ACC','ACA','ACG','UCU','UCC','UCA','UCG','AGU','AGC','UAU','UAC','UGG','CAA','CAG','AAU','AAC','CAU','CAC','GAA','GAG','GAU','GAC','AAA','AAG','CGU','CGC','CGA','CGG','AGA','AGG','UAA','UAG','UGA'];
