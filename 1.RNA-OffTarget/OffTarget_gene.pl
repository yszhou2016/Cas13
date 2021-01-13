$file=$ARGV[0];      #  output file of  Mismatch-search.pl
$file_gtf=$ARGV[0];  #  Hg38 gtf file
#########################
# Example: perl  OffTarget_gene.pl   Spacer.Mismatch.fa   GRCh38.gtf
#########################
open(AA,"$file");
open(BB,">$file.loc");
while($line=<AA>)
{  $line2=<AA>;
   @bb=split/\s+/,$line;
   if($bb[4]<=8 && $bb[4]>0){print BB "$bb[1]\t$bb[2]\t$bb[2]\t$bb[3]\t$bb[4]\n";}
}
close AA; close BB;

open(AA,"$file_gtf");
open(BB,">gene.loc");
while($line=<AA>)
{ chomp $line;
  $line=~s/"//g;  $line=~s/;//g;
  @bb=split/\s+/,$line;
   if($bb[2] eq "gene" )
    {      print BB  "$bb[0]\t$bb[3]\t$bb[4]\t$bb[9]\t$bb[6]","\n";
    }
}
close AA; close BB; 

############################################################################
`bedtools intersect -wo -a $file.loc  -b  gene.loc  > $file-gene.bed.tmp`;

open(FF,"$file-gene.bed.tmp");
open(FF2,">$file-gene.bed.filterted.tmp");
while($line=<FF>)
{ chomp $line;
  @bb=split/\s+/,$line;
  if($bb[3] ne $bb[9]){
        if(not exists($h{$bb[8]}))
        {print FF2 $line,"\n";}
  }
}
close FF; close FF2;

##############################
open(AA,"$file");
while($line=<AA>)
{  chomp $line;  $line2=<AA>;
   @bb=split/\s+/,$line;  
    chomp $line2;   
         $line2=~tr/ATCGatcg/TAGCTAGC/;      $line2=reverse($line2);
      $h{"$bb[1],$bb[2]"}=$line2;    
}
close AA;

open(FF,"$file-gene.bed.filterted2.tmp");
open(FF2,">$file.gene.list");
while($line=<FF>)
{  chomp $line;
    @bb=split/\s+/,$line;
        if($bb[3] eq "+"){$direction="-";}else{$direction="+";}
    print FF2  $h{"$bb[0],$bb[1]"},"\t","chr$bb[0]:$bb[1]\_$bb[2]\t",$direction,"\t",$bb[4],"\t",$bb[8],"\t",$bb[11],"\n";
}
close FF;   close FF2;
###############################
`rm *tmp`;
