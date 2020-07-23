### 1.align
hisat2 -p 18 --dta-cufflinks  -q -x   ~/data/hg19/Homo_sapiens.GRCh38.dna_sm.primary_assembly.fa  -1   $sample\_1.clean.fq.gz  -2   $sample\_2.clean.fq.gz  -S  $sample.sam
samtools view -Su  $sample.sam | samtools sort  -@ 18 - > $sample.sorted.bam

##### 2.REDItool
python /gpfsdata/home/yszhou/software/anaconda3/envs/py2/bin/REDItoolDenovo.py  -o  $sample.REDtools -i $sample.sorted.bam  -f  ~/data/hg19/Homo_sapiens.GRCh38.dna_sm.primary_assembly.fa  -t 24 -e -d -l -U [AG,TC,CT,GA] -p -u -m60 -T5-5 -W -v 1 -n  0.0
