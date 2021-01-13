###   1.align 
hisat2 -p 18 --dta-cufflinks  -q -x     ~/data/Homo_sapiens.GRCh38.dna_sm.primary_assembly.fa  -1   $sample\_1.clean.fq.gz  -2   $sample\_2.clean.fq.gz  -S  $sample.sam 
samtools view -Su -q 30  $sample.sam | samtools sort  -@ 18 - > $sample.sorted.bam  
samtools   index  $sample.sorted.bam

##### 2.Count Read
htseq-count -f bam  $sample.sorted.bam  ~/data/Homo_sapiens.GRCh38.96.chr.gtf  >  $sample.HTSeq.out  
