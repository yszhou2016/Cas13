#### All to All blast of S protein
blastall -p blastp -i S.pep.fa   -d  S.pep.fa  -e 1e-10  -m 8 -F F -b 50000 -v 50000 -a 20  -o  S-S.pep.tab

#### Filter alignment
cat S-S.pep.tab | awk '{if ($3》50) print}'  >  S-S.pep.new.tab

#### Clustering
mcl   S-S.pep.new.tab --abc - o   S-S.pep.new.tab.output
