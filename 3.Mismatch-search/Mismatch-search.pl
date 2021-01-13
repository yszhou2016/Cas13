$spacer=$ARGV[0];                    #  fasta file of query sequences
$reference=$ARGV[1];                 #  fasta file of reference genome for search
$MaxMismatch=$ARGV[2];               #  max mismatch between spacer and target 
###################################################
#  Example:  perl   Mismatch-search.pl     Spacer.fa    Hg38.genome.fa    8
###################################################
open(FASTA,"$reference");
open(FA,">genome.tmp");
$s=0;
while($line=<FASTA>)
{  chomp $line;
   if($line=~/^>/)
   { if($s>0){print FA "\n";}
       print FA $line,"\n";
     $s=1;
   } else{  $line=~tr/atcg/ATCG/;  print FA $line;}
 } 
 print FA "\n";
close FASTA; close FA;

##################################
open(AA,"$spacer");
while($line=<AA>)
{   $line2=<AA>;  
    chomp $line; $line=~s/>//;  @tt=split/\s+/,$line;  $target=$tt[0];
    chomp $line2;
    @bb=split//,$line2;
    $len=@bb;
    $s=0;
  ##############
  open(FF,"genome.tmp");
  open(BB,">$target.fa");
  while($lin=<FF>)
  { $lin2=<FF>;  
    chomp $lin;  $lin=~s/>//;  @ar=split/\s+/,$lin;
    chomp $lin2;
    $len2=length($lin2);
    $i=0;
    while($i<$len2)
      {
        $seq=substr($lin2,$i,$len+3);
        $seq2=$seq;   $seq2=~tr/ATCGatcg/TAGCTAGC/;  $seq2=reverse($seq2);
        @bb2=split//,$seq;
        @bb3=split//,$seq2;
        $j=0;  $t1=0;  $t2=0;
        while($j<$len)
            {  if($bb[$j] eq $bb2[$j]){$t1++;}
               if($bb[$j] eq $bb3[$j]){$t2++;}
               $j++;
            }
            $mismatch1=$len-$t1;
            $mismatch2=$len-$t2;  
            if($mismatch1<=$MaxMismatch){$s++;  print BB ">$target-$s\t$ar[0]\t$i\t+\t$mismatch1\n$seq\n";}
            if($mismatch2<=$MaxMismatch){$s++;  print BB ">$target-$s\t$ar[0]\t$i\t-\t$mismatch2\n$seq2\n";}
        $i++;
      }
    
  }
  close FF;  close BB;
  #########################################
}
close AA;  

