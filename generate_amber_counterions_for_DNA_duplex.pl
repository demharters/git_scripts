#! /usr/bin/perl -w 
#
# This script generates counterions around 
# a double stranded DNA.
#
#
# usage:
# generate_amber_counterions_for_DNA_duplex.pl DNA_duplex.pdb
# 
# standard output is the .pdb with counter ions
#

$nline = 0;
$natom = 0;
$nlatm = 0;
while (<>) {
      chomp;
      $entry = $_;
      $nline++;
      $entry[$nline] = $entry;
      if (substr($entry,0,4) eq 'ATOM') {
         $nlatm++;
         $x_atm[$nlatm] = substr($entry,30,8); $x_atm[$nlatm] =~ tr/ //d;
         $y_atm[$nlatm] = substr($entry,38,8); $y_atm[$nlatm] =~ tr/ //d;
         $z_atm[$nlatm] = substr($entry,46,8); $z_atm[$nlatm] =~ tr/ //d;
         @list_entry = split;
         $atom = $list_entry[2];
         if ($atom eq 'C1\'') {
             $altLoc = substr($entry,16,1);
             if ($altLoc eq ' ' || $altLoc eq 'A') {
                $natom++;
                $atomSeq[$natom] = substr($entry,6,5);
                $resName[$natom] = substr($entry,17,3);
                $chainID[$natom] = substr($entry,21,1); $chainID[$natom] =~ tr/ //d;
                $resSeq[$natom] = substr($entry,22,4); $resSeq[$natom] =~ tr/ //d;
                $x[$natom] = substr($entry,30,8); $x[$natom] =~ tr/ //d;
                $y[$natom] = substr($entry,38,8); $y[$natom] =~ tr/ //d;
                $z[$natom] = substr($entry,46,8); $z[$natom] =~ tr/ //d;
                $entry_ = $entry;
             }
         }
      }
      if (substr($entry,0,3) eq 'END' || substr($entry,0,3) eq 'TER') {
         last;
      }
}

#$chainList = "ABCDEFGHIJKLMNOPQ";
$nres = $natom;
#print "nres $nres\n";

# Assume that the DNA double helix has 1..nres/2 nres/2+1..nres
# composition of molecules 

# find the central level

 $nlevel = $nres/2;
 if (($nlevel % 2)==0) {
      $clevel = $nlevel/2;
 } else {
      $clevel = $nlevel/2 + 1; 
 }   

# get the helical parts helix[level][1..2], $nlevel = $nres/2

 $ilevel = 0;
 for ($ilevel=1;$ilevel<=$nlevel;$ilevel++) {
      $helix[$ilevel][1] = $ilevel;
      $helix[$ilevel][2] = $nres-$ilevel+1;
      $iA = $helix[$ilevel][1];
      $iB = $helix[$ilevel][2];
      $entry_A = $entry_;
      $entry_B = $entry_;

      $atomSeq = $atomSeq[$iA];
      substr($entry_A,6,5) = sprintf("%5s",$atomSeq);
      $resName = $resName[$iA];
      substr($entry_A,17,3) = $resName;
      $chainID = $chainID[$iA];
      substr($entry_A,21,1) = $chainID;
      $resSeq  = $resSeq[$iA];
      substr($entry_A,22,4) = sprintf("%4s",$resSeq);
      $x = $x[$iA];
      substr($entry_A,30,8) = sprintf("%8s",$x);
      $y = $y[$iA];
      substr($entry_A,38,8) = sprintf("%8s",$y);
      $z = $z[$iA];
      substr($entry_A,46,8) = sprintf("%8s",$z);
      $xA[$ilevel] = $x;
      $yA[$ilevel] = $y;
      $zA[$ilevel] = $z;

#      print "$entry_A\n";

      $atomSeq = $atomSeq[$iB];
      substr($entry_B,6,5) = sprintf("%5s",$atomSeq);
      $resName = $resName[$iB];
      substr($entry_B,17,3) = $resName;
      $chainID = $chainID[$iB];
      substr($entry_B,21,1) = $chainID;
      $resSeq  = $resSeq[$iB];
      substr($entry_B,22,4) = sprintf("%4s",$resSeq);;
      $x = $x[$iB];
      substr($entry_B,30,8) = sprintf("%8s",$x);
      $y = $y[$iB];
      substr($entry_B,38,8) = sprintf("%8s",$y);
      $z = $z[$iB];
      substr($entry_B,46,8) = sprintf("%8s",$z);
      $xB[$ilevel] = $x;
      $yB[$ilevel] = $y;
      $zB[$ilevel] = $z;

#      print "$entry_B\n";
      $xA = $xA[$ilevel];
      $yA = $yA[$ilevel];
      $zA = $zA[$ilevel];

      $xB = $xB[$ilevel];
      $yB = $yB[$ilevel];
      $zB = $zB[$ilevel];

      $xc[$ilevel] = $xA + 0.5*($xB - $xA);
      $yc[$ilevel] = $yA + 0.5*($yB - $yA);
      $zc[$ilevel] = $zA + 0.5*($zB - $zA);


#      print "$helix[$ilevel][1] $helix[$ilevel][2]\n";
 }

#NOUSE
# locate chainID[$iB]
# for ($i=0;$i<=length($chainList)-1;$i++) {
#     if (substr($chainList,$i,1) eq $chainID[$nres]) {$chainLoc = $i;}
# } 
#NOUSE

# Find Coordinates of Dummy atom D10 at level 1(b)  0 (a,cp)
        
          $xD10 = $xc[1];
          $yD10 = $yc[1];
          $zD10 = $zc[1];

# define a vector

      $ax = $xB[$clevel] - $xA[$clevel];
      $ay = $yB[$clevel] - $yA[$clevel];
      $az = $zB[$clevel] - $zA[$clevel];

      $a2 = $ax*$ax + $ay*$ay + $az*$az;

      $ax /= sqrt($a2);
      $ay /= sqrt($a2);
      $az /= sqrt($a2);

# define b vector for central level

      $bx = $xc[$clevel+1] - $xc[$clevel];
      $by = $yc[$clevel+1] - $yc[$clevel];   
      $bz = $zc[$clevel+1] - $zc[$clevel];

      $b2 = $bx*$bx + $by*$by + $bz*$bz;        

      $bx /= sqrt($b2);
      $by /= sqrt($b2);
      $bz /= sqrt($b2);

# define b vector for each level 1...level-1
 
      for ($ilevel=1;$ilevel<=$nlevel-1;$ilevel++) {
         $bx[$ilevel] = $xc[$ilevel+1] - $xc[$ilevel];
         $by[$ilevel] = $yc[$ilevel+1] - $yc[$ilevel];   
         $bz[$ilevel] = $zc[$ilevel+1] - $zc[$ilevel];
     
         $b2[$ilevel] = $bx*$bx + $by*$by + $bz*$bz;        
     
         $bx[$ilevel] /= sqrt($b2[$ilevel]);
         $by[$ilevel] /= sqrt($b2[$ilevel]);
         $bz[$ilevel] /= sqrt($b2[$ilevel]);
      }

# define cp vector

      $cpx = $ay*$bz - $az*$by;
      $cpy = $az*$bx - $ax*$bz;
      $cpz = $ax*$by - $ay*$bx;

      $cp2 = $cpx*$cpx + $cpy*$cpy + $cpz*$cpz;   

      $cpx /= sqrt($cp2);
      $cpy /= sqrt($cp2);
      $cpz /= sqrt($cp2);


# write out the pdb
  for ($i=1;$i<=$nline;$i++) {
      if ( (substr($entry[$i],0,3) ne 'END') && 
           (substr($entry[$i],0,3) ne 'TER') ) { 
            print "$entry[$i]\n";
      }#endif#
  }#endfor#


# Find Coordinates for counter ions 
# Define two counter ions per nucleotide level
# Write out counter-ions

    $atomSeq =  100;
    $resSeq  =  100;
    $scaleD  =  15.0;
    $min_dist = 2; # 2.5
    for ($ilevel=1;$ilevel<=$nlevel;$ilevel++) {

        #along DNA axis#
	if ($ilevel==1) {               
	    $xx = 0.0;
	    $yy = 0.0;
	    $zz = 0.0;
            $bx = $bx[1];
            $by = $by[1];
            $bz = $bz[1];
        } else {
            $xx += sqrt($b2[$ilevel-1])*$bx[$ilevel-1];  
            $yy += sqrt($b2[$ilevel-1])*$by[$ilevel-1];  
            $zz += sqrt($b2[$ilevel-1])*$bz[$ilevel-1];  
            $bx  = $bx[$ilevel-1];
            $by  = $by[$ilevel-1];
            $bz  = $bz[$ilevel-1];
        }#endif#

        #towards the cylinder#

        # define a vector
         $ax = $xB[$ilevel] - $xA[$ilevel];
         $ay = $yB[$ilevel] - $yA[$ilevel];
         $az = $zB[$ilevel] - $zA[$ilevel];
     
         $a2 = $ax*$ax + $ay*$ay + $az*$az;
     
         $ax /= sqrt($a2);
         $ay /= sqrt($a2);
         $az /= sqrt($a2);

       # define cp vector
         $cpx = $ay*$bz - $az*$by;
         $cpy = $az*$bx - $ax*$bz;
         $cpz = $ax*$by - $ay*$bx;
    
         $cp2 = $cpx*$cpx + $cpy*$cpy + $cpz*$cpz;   
    
         $cpx /= sqrt($cp2);
         $cpy /= sqrt($cp2);
         $cpz /= sqrt($cp2);

        $dmin = $min_dist - 1;  
        while ($dmin < $min_dist) { 
	 $xx_rest = $scaleD*(rand()-0.5)*($ax+$cpx);
	 $yy_rest = $scaleD*(rand()-0.5)*($ay+$cpy);
	 $zz_rest = $scaleD*(rand()-0.5)*($az+$cpz);

          $xD00 = $xD10 + $xx + $xx_rest;
          $yD00 = $yD10 + $yy + $yy_rest;
          $zD00 = $zD10 + $zz + $zz_rest;

          $dmin = 1000.0;
          for ($iatm=1;$iatm<=$nlatm;$iatm++) {
               $dx = $xD00 - $x_atm[$iatm];
               $dy = $yD00 - $y_atm[$iatm];
               $dz = $zD00 - $z_atm[$iatm];
               $dist = sqrt($dx*$dx+$dy*$dy+$dz*$dz);
               if ($dist<$dmin) {$dmin = $dist;}
          }
        }
        $nlatm++;
        $x_atm[$nlatm] = $xD00;
        $y_atm[$nlatm] = $yD00;
        $z_atm[$nlatm] = $zD00;
#FROM HELIX
#         $xx_rest = $scaleD*(rand()-0.5);
#         $yy_rest = $scaleD*(rand()-0.5);
#         $zz_rest = $scaleD*(rand()-0.5);         
#          $xD00 = $xA[$ilevel] + $xx_rest;
#          $yD00 = $yA[$ilevel] + $yy_rest;
#          $zD00 = $zA[$ilevel] + $zz_rest;
#FROM HELIX 

          $entry_D00 = $entry_;
          $atomName = "NA";
          substr($entry_D00,13,3) = sprintf("%3s",$atomName);
          $atomSeq++;
          substr($entry_D00,6,5) = sprintf("%5s",$atomSeq);
          $resName = "CIP";
          substr($entry_D00,17,3) = $resName;
          $chainID = ' ';
          substr($entry_D00,21,1) = $chainID;
          $resSeq++;
          substr($entry_D00,22,4) = sprintf("%4s",$resSeq);
          substr($entry_D00,30,8) = sprintf("%8.4g",$xD00);
          substr($entry_D00,38,8) = sprintf("%8.4g",$yD00);
          substr($entry_D00,46,8) = sprintf("%8.4g",$zD00);
          print "$entry_D00\n";


         #other strand
        $dmin = $min_dist - 1;  
        while ($dmin < $min_dist) { 
	 $xx_rest = $scaleD*(rand()-0.5)*($ax+$cpx);
	 $yy_rest = $scaleD*(rand()-0.5)*($ay+$cpy);
	 $zz_rest = $scaleD*(rand()-0.5)*($az+$cpz);

         $xD00 = $xD10 + $xx + $xx_rest;
         $yD00 = $yD10 + $yy + $yy_rest;
         $zD00 = $zD10 + $zz + $zz_rest;

          $dmin = 1000.0;
          for ($iatm=1;$iatm<=$nlatm;$iatm++) {
               $dx = $xD00 - $x_atm[$iatm];
               $dy = $yD00 - $y_atm[$iatm];
               $dz = $zD00 - $z_atm[$iatm];
               $dist = sqrt($dx*$dx+$dy*$dy+$dz*$dz);
               if ($dist<$dmin) {$dmin = $dist;}
          }
        }
        $nlatm++;
        $x_atm[$nlatm] = $xD00;
        $y_atm[$nlatm] = $yD00;
        $z_atm[$nlatm] = $zD00;

#FROM HELIX
#         $xx_rest = $scaleD*(rand()-0.5);
#         $yy_rest = $scaleD*(rand()-0.5);
#         $zz_rest = $scaleD*(rand()-0.5);
#          $xD00 = $xB[$ilevel] + $xx_rest;
#          $yD00 = $yB[$ilevel] + $yy_rest;
#          $zD00 = $zB[$ilevel] + $zz_rest;
#FROM HELIX 

          $entry_D00 = $entry_;
          $atomName = "NA";
          substr($entry_D00,13,3) = sprintf("%3s",$atomName);
          $atomSeq++;
          substr($entry_D00,6,5) = sprintf("%5s",$atomSeq);
          $resName = "CIP";
          substr($entry_D00,17,3) = $resName;
          $chainID = ' ';
          substr($entry_D00,21,1) = $chainID;
          $resSeq++;
          substr($entry_D00,22,4) = sprintf("%4s",$resSeq);
          substr($entry_D00,30,8) = sprintf("%8.4g",$xD00);
          substr($entry_D00,38,8) = sprintf("%8.4g",$yD00);
          substr($entry_D00,46,8) = sprintf("%8.4g",$zD00);
          print "$entry_D00\n";
         
    }#endfor#      
