# This project is about mining Cas13 proteins from metagenomic data and RNA off-target detection of Cas13/Cas13-ADAR2dd.
a.Identify the Cas13 proteins from metagenomic samples.

1. Install the “prodigal.linux” , "bedtools" and “pilercr” in the default environment.
   bedtools:  https://sourceforge.net/projects/bedtools/
   pilecr:  http://www.drive5.com/pilercr/

2. run “perl  0.Cas-Finder  $sample.fasta” to obtain the Cas proteins and generate "$sample.pep.cas.fasta" file.

3. run " perl  1.Cas13-Finder.pl   $sample.pep.cas.fasta" to otain the Cas13 proteins and generate "$sample.pep.cas.RxxxxH.fa" file.

4. Multiple alignment of Cas13 proteins with mafft.
mafft  --maxiterate 1000  --thread 12   --globalpair  Cas13.fa > Cas13.mafft.fasta


b. RNAseq off-Target analysis of Cas13.

1. align the RNAseq to reference genome
hisat2 -p 18 --dta-cufflinks  -q -x     GRCh38.genome.fa  -1   {$sample}_1.clean.fq.gz  -2   {$sample}_2.clean.fq.gz  -S  $sample.sam 
samtools view -Su -q 30  $sample.sam | samtools sort  -@ 18 - > $sample.sorted.bam  
samtools   index  $sample.sorted.bam

2. Calculate the read counts.
htseq-count -f bam  $sample.sorted.bam  GRCh38.gtf  >  $sample.HTSeq.out

3.Convert the read counts into FPKM values.
perl  HTSeq2FPKM.pl   $sample.HTSeq.out    GRCh38.gtf

4. Calculate differently expressed gene using DEseq2
perl    ~/path/trinityrnaseq-Trinity-v2.8.5/Analysis/DifferentialExpression/run_DE_analysis.pl  --matrix gene_count_matrix.csv.xls  --method DESeq2  --output  DESeq2_out  --samples_file  sample.txt  --contrasts  contrast.txt   --dispersion 0.04
perl  4.get_DEG-DESeq2.pl  gene_count_matrix.csv.xls.Control_vs_Treat.DESeq2.DE_results

#5. Predict the off-target site of Hg38 genome and transcriptome of the spacers with no more than eight mismatches.
perl  Mismatch-search.pl  Spacer.fa    Hg38.genome.fa    8
perl  OffTarget_gene.pl   Spacer.Mismatch.fa   Hg38.gtf


c. RNAseq off-target analysis of RNA base editor. 
1.align the RNAseq to reference genome
hisat2 -p 18 --dta-cufflinks  -q -x   GRCh38.genome.fa  -1   $sample\_1.clean.fq.gz  -2   $sample\_2.clean.fq.gz  -S  $sample.sam
samtools view -Su  $sample.sam | samtools sort  -@ 18 - > $sample.sorted.bam

2.Install the REDItool fom "https://sourceforge.net/projects/reditools/" .

python ~/software/anaconda3/envs/py2/bin/REDItoolDenovo.py  -o  $sample.REDtools -i $sample.sorted.bam  -f  GRCh38.genome.fa  -t 24 -e -d -l -U [AG,TC,CT,GA] -p -u -m60 -T5-5 -W -v 1 -n  0.0

