$FileCount=$ARGV[0];    #  Output file of HTSeq-count
$FileGtf=$ARGV[1];      #  /home/yszhou/data/genome/GRCh38-ensembl/Homo_sapiens.GRCh38.96.chr.gtf   
#################
#example:  perl  OffTarget_gene.pl   $sample.HTSeq.out    Homo_sapiens.GRCh38.96.chr.gtf
#################
open(AA,"$FileGtf");
open(BB,">tmp.exon.loc");
open(BB2,">tmp.gene.list");
while($line=<AA>)
{  chomp $line;  $line=~s/\"//g; $line=~s/;//g;
   @bb=split/\s+/,$line; 
   if($bb[2] eq "exon"){  print BB $bb[0],"\t",$bb[3],"\t",$bb[4],"\t",$bb[9],"\n";
        }
   if($bb[2] eq "gene"){
           if(not exists($v{$bb[9]})){print BB2  $bb[9],"\t",$bb[13],"\n";}  $v{$bb[9]}=1;
        }
}close AA; close BB;  close BB2;
###########################################################################################################
`sort   -dk4,4  -k1,1n  -k2,2n  -k3,3n   tmp.exon.loc   >  tmp.exon.loc.sorted `;
###########################################################################################################
$s=0;
open(CC,"tmp.exon.loc.sorted ");
open(DD,">tmp.exon.loc.sorted.filtered");
while($line=<CC>)
{   chomp $line;
    @bb=split/\s+/,$line;
     if($s>0){
            if($bb[3] eq $pre)
                 {
                   if($bb[1]<$end)
                              {$end=$bb[2];}
                    else{  print  DD  $chr,"\t",$start,"\t",$end,"\t",$pre,"\n";  
                                  if(not exists($len{$pre})){$len{$pre}=$end-$start+1;}
                                    else{$len{$pre}=$end-$start+1+$len{$pre};}
                                $chr=$bb[0];  $start=$bb[1];  $end=$bb[2];   $pre=$bb[3];
                            }
                  }
                else{  print  DD  $chr,"\t",$start,"\t",$end,"\t",$pre,"\n";
                               if(not exists($len{$pre})){$len{$pre}=$end-$start+1;}
                                    else{$len{$pre}=$end-$start+1+$len{$pre};}
                            $chr=$bb[0];  $start=$bb[1];  $end=$bb[2];   $pre=$bb[3];
                         }
            }
          else{  $chr=$bb[0];  $start=$bb[1];  $end=$bb[2];   $pre=$bb[3];  }
       $s++;
}
print  DD  $chr,"\t",$start,"\t",$end,"\t",$pre,"\n";
 if(not exists($len{$pre})){$len{$pre}=$end-$start+1;}
 else{$len{$pre}=$end-$start+1+$len{$pre};}
close CC;  close DD;

open(FF,"tmp.gene.list");
open(FF2,">gene.len");
while($line=<FF>)
{  chomp $line;
   @bb=split/\s+/,$line;
   print FF2  $bb[0],"\t",$len{$bb[0]},"\t",$bb[1],"\n";
}close FF; close FF2;

`rm  tmp*`;
#####################################################################################################
open(FF,"gene.len");
while($line=<FF>)
{  chomp $line;
   @bb=split/\s+/,$line;
   $Len{$bb[0]}=$bb[1];
    $Name{$bb[0]}=$bb[2];
}
close FF;

$sum=0;
open(FF,"$FileCount");
while($line=<FF>)
{  chomp $line;
   @bb=split/\s+/,$line;
   if($bb[0]=~/^E/)
     {$sum=$sum+$bb[1];}
}
close FF;

open(FF,"$FileCount");
open(FF2,">$FileCount.fpkm");
while($line=<FF>)
{  chomp $line;
   @bb=split/\s+/,$line;
   if($bb[0]=~/^E/)
     {$fpkm=$bb[1]*1e9/($Len{$bb[0]}*$sum);
             print FF2  $bb[0],"\t",$fpkm,"\t",$Name{$bb[0]},"\n";
          }
}
close FF;  close FF2;
