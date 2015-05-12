#! /usr/bin/perl -w 
#
# usage:
#  ./methyl2formyl_all_CME.pl 5Me-CYT-DNA.pdb > 5f-CYT-DNA.pdb
#
# This script changes 5MeC residues to 5fC in a DNA 
#
#                      H42  H41
#                        \  /
#                         N4
#                         |
#      H51             O  C4
#       |             || / \\
#  H52-C5M-        H5-C-C5   N3
#       |               ||   |
#      H53              C6   C2
#                      / \  / \\
#                     H6  N1   O2
#                           \
#                            \
#                             \
#         O1P    H5' H4'  O4'  \
#          |      |    \ /   \  \
#         -P-O5'-C5'---C4'    C1'
#          |      |     \     / \
#         O2P    H5''   C3'--C2' H1'
#                      / \   / \
#                   O3' H3' O2' H2''
#                    |       |
#                           H2'
# Copyright (c) 2007 Peter Minary
                                                #      "." atom type in AMBER
  $r_bond_C5_C5M        =  1.51;                # in A "CM"-"CT" distance in AMBER
  $r_bond_C5_C          =  1.444;               # in A "CM"-"C" distance in AMBER 
  $r_bond_C_O           =  1.229;               # in A "C"="O" distance in AMBER
  $r_bond_C_H5          =  1.08;                # in A "C"-"H4" distance in AMBER
  $theta_C5_C_O         =  125.3;               # in deg "CM"-"C"-"O" angle in AMBER
  $theta_C5_C_H5        =  115.0;               # in deg "CM"-"C"-"H4" angle in AMBER
  $pi                   =  3.14;     
  $phi_C6_C5_C_O        =  0.0;      # in rad
  $phi_C6_C5_C_H5       =  $pi;      # in rad 

  $nline  = 0;
  while (<>) {
        chomp;
        $entry = $_;
        if ( (substr($entry,0,3) eq 'TER') ||
             (substr($entry,0,3) eq 'END')  ) {
              $nline++;
              $entry[$nline]    = $entry;
        }
        if ( (substr($entry,0,4) eq 'ATOM')  ||
             (substr($entry,0,6) eq 'HETATM') ) { 
          $nline++;
          $entry[$nline]    = $entry;
          $atom_num[$nline] = substr($entry,6,5);  $atom_num[$nline]=~ tr/ //d;
          $atom[$nline]     = substr($entry,12,4); $atom[$nline]    =~ tr/ //d;
          $resName[$nline]  = substr($entry,17,3); $resName[$nline] =~ tr/ //d;
          $resSeq[$nline]   = substr($entry,22,4); $resSeq[$nline]  =~ tr/ //d;
#NOUSE
#         $chain_id[$nline] = substr($entry,21,1);
#NOUSE
          $x[$nline]  = substr($entry,30,8); $x[$nline] =~ tr/ //d;
          $y[$nline]  = substr($entry,38,8); $y[$nline] =~ tr/ //d;
          $z[$nline]  = substr($entry,46,8); $z[$nline] =~ tr/ //d;
#NOUSE
#         $char_first[$nline] = substr($atom,0,1);
#         $length[$nline]     = length($atom);
#NOUSE
        }
  }

  $resSeq_CME  = -100;
  $insert_line = 0;
  for ($iline=1;$iline<=$nline;$iline++) {
         $entry = $entry[$iline];
         $entry_new = $entry;

        if ( (substr($entry,0,3) eq 'TER') ||
             (substr($entry,0,3) eq 'END')  ) {
             next;
        }

         if (  ($resName[$iline] eq 'CME') && 
               ($resSeq_CME != $resSeq[$iline]) ) {
             $resSeq_CME = $resSeq[$iline];
             $natm_CME  = 0;
             for ($jline=$iline;$jline<=$nline;$jline++) {
                  $entry_jline = $entry[$jline];
                  if ( (substr($entry_jline,0,3) eq 'TER') ||
                       (substr($entry_jline,0,3) eq 'END')  ) {
                        next;
                  }
                  if ($resSeq[$jline] eq  $resSeq_CME) {
                      if ( substr($entry_jline,0,4) eq 'ATOM') { 
                           substr($entry_jline,17,3) = 'CFO';
                           $entry[$jline] = $entry_jline;
                      }
                      $natm_CME++;
                      if ($atom[$jline] eq 'C5') {
                          $xC5 = $x[$jline];
                          $yC5 = $y[$jline];
                          $zC5 = $z[$jline];
                      }
                      if ($atom[$jline] eq 'C5M') {
                          $xC5M = $x[$jline];
                          $yC5M = $y[$jline];
                          $zC5M = $z[$jline];
                      }
                      if ($atom[$jline] eq 'C6') {
                          $xC6 = $x[$jline];
                          $yC6 = $y[$jline];
                          $zC6 = $z[$jline];
                      }
                  } 
                  else {last;}
             }

             # build C 
               $ratio = $r_bond_C5_C/$r_bond_C5_C5M;
               $xC = $xC5 + $ratio*($xC5M - $xC5);
               $yC = $yC5 + $ratio*($yC5M - $yC5);
               $zC = $zC5 + $ratio*($zC5M - $zC5);
               
             # build O       

                 # use procedure applied within DC-REPSWA
                 # define r1=C6, r2=C5, r3=C, r, theta, phi
                 $r1x = $xC6; $r1y = $yC6; $r1z = $zC6;
                 $r2x = $xC5; $r2y = $yC5; $r2z = $zC5;
                 $r3x = $xC;  $r3y = $yC;  $r3z = $zC;

                 $r     = $r_bond_C_O;
                 $theta = $pi*(180.0-$theta_C5_C_O)/180.0;
               
                 # O
                 $phi   = $phi_C6_C5_C_O;
                 # build r4 based on r1, r2, r3, r, theta, phi
                 &build_r4;
                 # read out O
                 $xO = $r4x; $yO = $r4y; $zO = $r4z;                    

             # build H5       

                 # use procedure applied within DC-REPSWA
                 # define r1=C6, r2=C5, r3=C, r, theta, phi
                 $r1x = $xC6; $r1y = $yC6; $r1z = $zC6;
                 $r2x = $xC5; $r2y = $yC5; $r2z = $zC5;
                 $r3x = $xC;  $r3y = $yC;  $r3z = $zC;

                 $r     = $r_bond_C_H5;
                 $theta = $pi*(180.0-$theta_C5_C_H5)/180.0;
               
                 # O
                 $phi   = $phi_C6_C5_C_H5;
                 # build r4 based on r1, r2, r3, r, theta, phi
                 &build_r4;
                 # read out H5
                 $xH5 = $r4x; $yH5 = $r4y; $zH5 = $r4z;                    

                 $insert_line = $iline + $natm_CME;

	 }#endif-CME#

         $entry = $entry[$iline];
#         $entry_new = $entry;

         if ($insert_line == $iline) {

             $ilinem1   = $iline - 1;
             $entry_new = $entry[$ilinem1];

             $atom_new = ' C  ';
             substr($entry_new,12,4) = $atom_new;
             $x_field = sprintf("%8.5g",$xC);
             $y_field = sprintf("%8.5g",$yC);
             $z_field = sprintf("%8.5g",$zC);
             substr($entry_new,30,8) = $x_field;
             substr($entry_new,38,8) = $y_field;
             substr($entry_new,46,8) = $z_field;
             substr($entry_new,6,5)  = '    0';
#             substr($entry_new,70,15) = '   C';
             print "$entry_new\n";

             $atom_new = ' O  ';
             substr($entry_new,12,4) = $atom_new;
             $x_field = sprintf("%8.5g",$xO);
             $y_field = sprintf("%8.5g",$yO);
             $z_field = sprintf("%8.5g",$zO);
             substr($entry_new,30,8) = $x_field;
             substr($entry_new,38,8) = $y_field;
             substr($entry_new,46,8) = $z_field;
             substr($entry_new,6,5)  = '    0';
#             substr($entry_new,70,15) = '   O';
             print "$entry_new\n";

             $atom_new = 'H5  ';
             substr($entry_new,12,4) = $atom_new;
             $x_field = sprintf("%8.5g",$xH5);
             $y_field = sprintf("%8.5g",$yH5);
             $z_field = sprintf("%8.5g",$zH5);
             substr($entry_new,30,8) = $x_field;
             substr($entry_new,38,8) = $y_field;
             substr($entry_new,46,8) = $z_field;
             substr($entry_new,6,5)  = '    0';
#             substr($entry_new,70,15) = '   H';
             print "$entry_new\n";
             print "$entry\n";
         }
         else {
            $resName[$iline] = substr($entry[$iline],17,3); 
            $resName[$iline] =~ tr/ //d;
            if ($resName[$iline] eq 'CFO') { 
               if ((substr($entry,12,4) ne 'H51 ') &&
                   (substr($entry,12,4) ne ' H51') &&
                   (substr($entry,12,4) ne 'H52 ') &&
                   (substr($entry,12,4) ne ' H52') &&
                   (substr($entry,12,4) ne 'H53 ') &&
                   (substr($entry,12,4) ne ' H53') &&
                   (substr($entry,12,4) ne 'C5M ') &&
                   (substr($entry,12,4) ne ' C5M')) {print "$entry\n";}
            } else {
                    print "$entry\n"; 
            }
         }
         


#<STDIN>;
     
 }
 print "END\n";


 sub build_r4 
 {
     # builds r4 based on (r1, r2, r3, r, theta, phi)
     # defined in the main body of the program 

     # define e1 = (r3-r2)/||r3-r2||, e2 = (r1-r2)/||r1-r2||
     $e1x = $r3x - $r2x; $e1y = $r3y - $r2y; $e1z = $r3z - $r2z; 
     $e2x = $r1x - $r2x; $e2y = $r1y - $r2y; $e2z = $r1z - $r2z;

     $e1_2 = $e1x*$e1x + $e1y*$e1y + $e1z*$e1z;
     $e1   = sqrt($e1_2);
     $e1x  /= $e1; $e1y /= $e1; $e1z /= $e1;                 

     $e2_2 = $e2x*$e2x + $e2y*$e2y + $e2z*$e2z;
     $e2   = sqrt($e2_2);
     $e2x  /= $e2; $e2y /= $e2; $e2z /= $e2;                 

     # define ez = e1
     $ez_x = $e1x; $ez_y = $e1y; $ez_z = $e1z;
     
     # define ey = (e1 x e2) / ||e1 x e2||
     $ey_x  = $e1y*$e2z - $e1z*$e2y;
     $ey_y  = $e1z*$e2x - $e1x*$e2z;
     $ey_z  = $e1x*$e2y - $e1y*$e2x;

     $ey_2 = $ey_x*$ey_x + $ey_y*$ey_y + $ey_z*$ey_z;
     $ey   = sqrt($ey_2);
     $ey_x  /= $ey; $ey_y /= $ey; $ey_z /= $ey;                 
     
     # define ex = (ey x ez) / ||ey x ez||
     $ex_x  = $ey_y*$ez_z - $ey_z*$ez_y;
     $ex_y  = $ey_z*$ez_x - $ey_x*$ez_z;
     $ex_z  = $ey_x*$ez_y - $ey_y*$ez_x;

     $ex_2 = $ex_x*$ex_x + $ex_y*$ex_y + $ex_z*$ex_z;
     $ex   = sqrt($ex_2);
     $ex_x  /= $ex; $ex_y /= $ex; $ex_z /= $ex;                 

     #reconstruct in reference frame (ex, ey, ez)
     $rx    = $r*sin($theta)*cos($phi);
     $ry    = $r*sin($theta)*sin($phi);
     $rz    = $r*cos($theta);
     
     #rotate from reference frame to real frame
     $rrx   = $ex_x*$rx + $ey_x*$ry + $ez_x*$rz;
     $rry   = $ex_y*$rx + $ey_y*$ry + $ez_y*$rz;
     $rrz   = $ex_z*$rx + $ey_z*$ry + $ez_z*$rz;
     
     #place vector by adding it to $r3
     $r4x = $r3x + $rrx; 
     $r4y = $r3y + $rry; 
     $r4z = $r3z + $rrz;
 }
