wget   http://ccb.jhu.edu/software/stringtie/dl/prepDE.py
chmod 755  prepDE.py
./prepDE.py  -i StringTie/ -g StringTie-new/gene_count_matrix.csv -t StringTie-new/transcript_count_matrix.csv
