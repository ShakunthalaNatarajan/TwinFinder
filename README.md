# TwinFinder
Automatic detection of syntenic tandem gene duplications

Script details:

The features of TwinFinder are detailed below.

Dependencies: NCBI BLAST, Python modules (pandas, numpy, operator, collections, copy, seaborn)

Input files: Whole genome FASTA files and GFF3 (annotation) files of the query (target) and the reference organisms

Output (Key result) files: List of true tandem duplications in the query; Syntelog tandem duplications; Non-syntelog tandem duplications; Histogram depicting the number of genes in the detected tandem groups and their mutual average similarity scores 

Adjustable parameters: 

•	Similarity threshold value (the minimum percentage identity score that must be satisfied by genes that are detected to belong to the potential tandem duplication groups identified

•	Number of flanking genes to be looked at in query organism (flankq) – number of flanking genes to be considered in the query organism

•	Number of flanking genes to be looked at in reference organism (flankr) – number of flanking genes to be considered in the reference organism

•	Synteny ratio – Number of best forward blast hits/ (2*flankq)

•	Synteny cutoff – minimum cutoff ratio to be fulfilled for an identified tandem group of genes and their corresponding anchor gene in the reference organism to be syntenic 

Recommendations: Defining a perfect fitting value for flankq and flankr and finding the best fit synteny cutoff is a tricky step. The optimum values of each of these parameters could differ from organism to organism depending on its inherent genetic makeup and other complexities. Hence these values could be tweaked around and optimum values can be obtained through iterative trial and error steps. It is recommended that flankr be the same as flankq or be somewhat closer to it in the upper range. This would help reduce false positives in the respective results files.
