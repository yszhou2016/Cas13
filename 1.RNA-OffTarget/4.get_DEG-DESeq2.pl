$file=$ARGV[0];     # output file of DESeq2
########################################
# example: perl 4.get_DEG-DESeq2.pl  gene_count_matrix.csv.xls.Control_vs_Treat.DESeq2.DE_results
########################################
open(FF,">DEG.list.xls");
print FF "gene_ID\tlog2FC\t-log10(FDR)\tnone-up-down\tname\n";

open(AA,"$file");
open(BB,">DE.$file");
open(CC,">DEplot.$file");
print BB "gene_ID\tlog2FC\t-log10(FDR)\tnone-up-down\tname\n";
print CC "gene_ID\tgene1\tgene2\tnone-up-down\n";

open(DD,">UP-regulated.list");
open(EE,">Down-regulated.list");
$line=<AA>;
$s=0;
while($line=<AA>)
{  chomp $line;
   @bb=split/\s+/,$line;
   if($bb[6]>1 && $bb[10]<0.05){$temp="down";  print EE $bb[0],"\n";}
    elsif($bb[6]<-1 && $bb[10]<0.05){$temp="up";  print DD $bb[0],"\n";}
     else{$temp="none";}
   $fdr=-log($bb[10]+1e-200)/log(10);
     $tt=$bb[6]*-1;
    @b=split/\|/,$bb[0];
   print CC "$b[0]\t$bb[3]\t$bb[4]\t$temp\t$b[1]\n";   
   print BB "$b[0]\t$tt\t$fdr\t$temp\t$b[1]\n"; 

  if($bb[6]>1 && $bb[10]<0.05){$temp="down";  print FF "$b[0]\t$tt\t$fdr\t$temp\t$b[1]\n";}
  elsif($bb[6]<-1 && $bb[10]<0.05){$temp="up";  print FF "$b[0]\t$tt\t$fdr\t$temp\t$b[1]\n";}

   
}
#print $s,"\n";
close AA; close BB; close CC; close DD; close EE; 
