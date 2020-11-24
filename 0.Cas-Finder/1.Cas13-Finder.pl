$file=$ARGV[0];
############################
# run " perl  1.Cas13-Finder.pl   $sample.pep.cas.fasta"
############################
$file2=$file;  $file2=~s/.fasta//;
open(FASTA,"$file");
open(FA,">$file2.fa");
$s=0;
while($line=<FASTA>)
{  chomp $line;
  if($line=~/^>/)
   { if($s>0){print FA "\n";}
     print FA $line,"\n";
     $s=1;
   }
    else{print FA $line;}
 }
 print FA "\n";
close FASTA; close FA;
###################################################
open(AA,"$file2.fa");
open(BB,">$file2.RxxxxH.fa");
while($line=<AA>)
{  chomp $line;
   $line2=<AA>;  chomp $line2;
   $ll=$line2;
   $len=length($ll);
   @bb=split//,$ll;
   $i=0;  
   $s1=0;  $s2=0;   $t1=0;  $t2=0;
   while($i<$len)
   {  if($bb[$i] eq "R" && $bb[$i+1] eq "H" && $bb[$i+5] eq "H")
        {$s1++;
            if($t1==0){$t1=$i;}
            else{$t2=$i;} 
        }
	  elsif($bb[$i] eq "R" &&  $bb[$i+1] eq "N" && $bb[$i+5] eq "H")
	    {$s2++;}
            elsif($bb[$i] eq "R" &&  $bb[$i+1] eq "Q" && $bb[$i+5] eq "H")
	    {$s2++;}  
      $i++;
   }
####################   
    if(not exists($h{$line2})){ 
      if(($s1+$s2)>=2){
        if($t1*2<$len && $t2*2>$len)
         {  print BB $line,"\_L",length($line2),"\n",$line2,"\n";  
            $h{$line2}=1;
         }
        }
      }
}
close AA; close BB;
