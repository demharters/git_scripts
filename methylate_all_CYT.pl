#! /usr/bin/perl -w 
#
# usage:
#  ./methylate_CYT.pl DNA.pdb > 5MeC-DNA.pdb
#
# This script methylates all CYT residues in a DNA 
#
#                      H42  H41
#                        \  /
#                         N4
#                         |
#                  H51    C4
#                   |    / \\
#       H5-     H52-C5M-C5   N3
#                   |   ||   |
#                  H53  C6   C2
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

  $r_bond_C5_C5M        =  1.48;                # in A
  $r_bond_C5M_H         =  1.11;                # in A
  $theta_C5_C5M_H       =  110.10;              # in deg
  $phi_C6_C5_C5M_H51    =     3.14*120.0/180.0; # in rad
  $phi_C6_C5_C5M_H52  = 2.0*3.14*120.0/180.0; # in rad 
  $phi_C6_C5_C5M_H53  = 3.0*3.14*120.0/180.0; # in rad

  $pi                 =  3.14;     # in rad

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

  $resSeq_CYT  = -100;
  $insert_line = 0;
  for ($iline=1;$iline<=$nline;$iline++) {
         $entry = $entry[$iline];
         $entry_new = $entry;

        if ( (substr($entry,0,3) eq 'TER') ||
             (substr($entry,0,3) eq 'END')  ) {
             next;
        }

         if ( (($resName[$iline] eq 'CYT') ||
               ($resName[$iline] eq 'C'  )) &&  
              ($resSeq_CYT != $resSeq[$iline]) ) {
             $resSeq_CYT = $resSeq[$iline];
             $natm_CYT  = 0;
             for ($jline=$iline;$jline<=$nline;$jline++) {
                  $entry_jline = $entry[$jline];
                  if ( (substr($entry_jline,0,3) eq 'TER') ||
                       (substr($entry_jline,0,3) eq 'END')  ) {
                        next;
                  }
                  if ($resSeq[$jline] eq  $resSeq_CYT) {
                      if ( substr($entry_jline,0,4) eq 'ATOM') { 
                           substr($entry_jline,17,3) = 'CME';
                           $entry[$jline] = $entry_jline;
                      }
                      $natm_CYT++;
                      if ($atom[$jline] eq 'C5') {
                          $xC5 = $x[$jline];
                          $yC5 = $y[$jline];
                          $zC5 = $z[$jline];
                      }
                      if ($atom[$jline] eq 'H5') {
                          $xH5 = $x[$jline];
                          $yH5 = $y[$jline];
                          $zH5 = $z[$jline];
                      }
                      if ($atom[$jline] eq 'C6') {
                          $xC6 = $x[$jline];
                          $yC6 = $y[$jline];
                          $zC6 = $z[$jline];
                      }
                  } 
                  else {last;}
             }

             # build C5M
               $ratio = $r_bond_C5_C5M/$r_bond_C5M_H;
               $xC5M = $xC5 + $ratio*($xH5 - $xC5);
               $yC5M = $yC5 + $ratio*($yH5 - $yC5);
               $zC5M = $zC5 + $ratio*($zH5 - $zC5);
               
             # build H51, H52 H53       

                 # use procedure applied within DC-REPSWA
                 # define r1=C6, r2=C5, r3=C5M, r, theta, phi
                 $r1x = $xC6; $r1y = $yC6; $r1z = $zC6;
                 $r2x = $xC5; $r2y = $yC5; $r2z = $zC5;
                 $r3x = $xC5M; $r3y = $yC5M; $r3z = $zC5M;

                 $r     = $r_bond_C5M_H;
                 $theta = $pi*(180.0-$theta_C5_C5M_H)/180.0;
               
                 # H51
                 $phi   = $phi_C6_C5_C5M_H51;
                 # build r4 based on r1, r2, r3, r, theta, phi
                 &build_r4;
                 # read out H51
                 $xH51 = $r4x; $yH51 = $r4y; $zH51 = $r4z;                    

                 # H52
                 $phi   = $phi_C6_C5_C5M_H52;
                 # build r4 based on r1, r2, r3, r, theta, phi
                 &build_r4;
                 # read out H52
                 $xH52 = $r4x; $yH52 = $r4y; $zH52 = $r4z;                    

                 # H53
                 $phi   = $phi_C6_C5_C5M_H53;
                 # build r4 based on r1, r2, r3, r, theta, phi
                 &build_r4;
                 # read out H53
                 $xH53 = $r4x; $yH53 = $r4y; $zH53 = $r4z;                    

                 $insert_line = $iline + $natm_CYT - 1;

	 }#endif-CYT#

         $entry = $entry[$iline];
         $entry_new = $entry;

         if ($insert_line == $iline) {
             $atom_new = 'C5M ';
             substr($entry_new,12,4) = $atom_new;
             $x_field = sprintf("%8.5g",$xC5M);
             $y_field = sprintf("%8.5g",$yC5M);
             $z_field = sprintf("%8.5g",$zC5M);
             substr($entry_new,30,8) = $x_field;
             substr($entry_new,38,8) = $y_field;
             substr($entry_new,46,8) = $z_field;
             substr($entry_new,6,5)  = '    0';
#             substr($entry_new,70,15) = '   C';
             print "$entry_new\n";

             $atom_new = 'H51 ';
             substr($entry_new,12,4) = $atom_new;
             $x_field = sprintf("%8.5g",$xH51);
             $y_field = sprintf("%8.5g",$yH51);
             $z_field = sprintf("%8.5g",$zH51);
             substr($entry_new,30,8) = $x_field;
             substr($entry_new,38,8) = $y_field;
             substr($entry_new,46,8) = $z_field;
             substr($entry_new,6,5)  = '    0';
#             substr($entry_new,70,15) = '   H';
             print "$entry_new\n";

             $atom_new = 'H52 ';
             substr($entry_new,12,4) = $atom_new;
             $x_field = sprintf("%8.5g",$xH52);
             $y_field = sprintf("%8.5g",$yH52);
             $z_field = sprintf("%8.5g",$zH52);
             substr($entry_new,30,8) = $x_field;
             substr($entry_new,38,8) = $y_field;
             substr($entry_new,46,8) = $z_field;
             substr($entry_new,6,5)  = '    0';
#             substr($entry_new,70,15) = '   H';
             print "$entry_new\n";

             $atom_new = 'H53 ';
             substr($entry_new,12,4) = $atom_new;
             $x_field = sprintf("%8.5g",$xH53);
             $y_field = sprintf("%8.5g",$yH53);
             $z_field = sprintf("%8.5g",$zH53);
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
            if ($resName[$iline] eq 'CME') { 
               if ((substr($entry,12,4) ne 'H5  ') &&
                   (substr($entry,12,4) ne ' H5 ') && 
                   (substr($entry,12,4) ne '  H5')) {print "$entry\n";}
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
