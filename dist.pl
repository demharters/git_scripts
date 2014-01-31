#!/usr/bin/perl -w
#
# this script gets the distribution 
#
# usage: 
# ./dist.pl dmin dmax resol data column > dist
#
 $dmin   = $ARGV[0];
 $dmax   = $ARGV[1];
 $resol  = $ARGV[2];
 $data   = $ARGV[3];
 $column = $ARGV[4];

 open(DATA,"$data");

 $hdim = ($dmax-$dmin)/$resol;
 for ($i=1;$i<=$hdim;$i++) {
     $dist[$i] = 0.0;
 }
 $sum_res = 0;
 while (<DATA>) {
       chomp;
       @a = split;       
       $w = $a[$column-1];
       if ($w > $dmax) {next;}
       $ic = int(($w-$dmin)/$resol) + 1;
       $sum_res += $w;
       $dist[$ic]++;
 }
# print "$sum_res\n";
 $sum = 0.0;
 for ($i=1;$i<=$hdim;$i++) {
     $sum += $dist[$i]*$resol;
     $sum[$i] = 0.0;
 }
 for ($i=1;$i<=$hdim;$i++) {
     $dist[$i] /= $sum;
 }
 for ($j=1;$j<=$hdim;$j++) {
     for ($i=1;$i<=$j;$i++) {
         $sum[$j] += $dist[$i]*$resol;
     }
 } 
 for ($i=1;$i<=$hdim;$i++) {
     $x = $dmin+$i*$resol;
     print "$x $dist[$i] $sum[$i]\n";
 }   
