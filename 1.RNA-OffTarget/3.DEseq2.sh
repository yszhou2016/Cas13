perl    ~/yszhou/software/trinityrnaseq-Trinity-v2.8.5/Analysis/DifferentialExpression/run_DE_analysis.pl  --matrix gene_count_matrix.csv.xls  --method DESeq2  --output  DESeq2_out  --samples_file  sample.txt  --contrasts  contrast.txt   --dispersion 0.04