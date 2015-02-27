#! /usr/bin/perl -w 
# usage:
#  ./hydroxylate_all_AME.pl 6MeA-RNA.pdb > 6MeOH-A-RNA.pdb 
#
# This script hydroxylates all 6-MeA residues in a DNA 
#
#
#
#
#                    
#    H51     HO-OH   H62                     
#     |         |    /                         
#           H53-C5M-N6 
#               |   |                       
#              H52  C6 
#                  // \ 
#                 N1  C5--N7\\
#                 |   ||     C8-H8 
#                 C2  C4--N9/
#                / \\ /     \
#               H2  N3       \
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
# Copyright (c) 2014 Peter Minary

  $r_bond_C5M_OH        =  1.41;                # in A "CT"-"OH" distance in AMBER
  $r_bond_C5M_H         =  1.11;                # in A
  $r_bond_OH_HO         =  0.96;                # in A "OH"-"HO" distance in AMBER
  $theta_C5M_OH_HO      =  108.5;               # in deg
  $phi_N6_C5M_OH_HO     = 2.0*3.14*120.0/180.0; # in rad   change 120, remove 2:

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

  $resSeq_AME  = -100;
  $insert_line = 0;
  for ($iline=1;$iline<=$nline;$iline++) {
         $entry = $entry[$iline];
         $entry_new = $entry;

        if ( (substr($entry,0,3) eq 'TER') ||
             (substr($entry,0,3) eq 'END')  ) {
             next;
        }

         if (  ($resName[$iline] eq 'AME') && 
               ($resSeq_AME != $resSeq[$iline]) ) {
             $resSeq_AME = $resSeq[$iline];
             $natm_AME  = 0;
             for ($jline=$iline;$jline<=$nline;$jline++) {
                  $entry_jline = $entry[$jline];
                  if ( (substr($entry_jline,0,3) eq 'TER') ||
                       (substr($entry_jline,0,3) eq 'END')  ) {
                        next;
                  }
                  if ($resSeq[$jline] eq  $resSeq_AME) {
                      if ( substr($entry_jline,0,4) eq 'ATOM') { 
                           substr($entry_jline,17,3) = 'AHM';
                           $entry[$jline] = $entry_jline;
                      }
                      $natm_AME++;
                      if ($atom[$jline] eq 'N6') {
                          $xN6 = $x[$jline];
                          $yN6 = $y[$jline];
                          $zN6 = $z[$jline];
                      }
                      if ($atom[$jline] eq 'C5M') {
                          $xC5M = $x[$jline];
                          $yC5M = $y[$jline];
                          $zC5M = $z[$jline];
                      }
                      if ($atom[$jline] eq 'H51') {
                          $xH51 = $x[$jline];
                          $yH51 = $y[$jline];
                          $zH51 = $z[$jline];
                      }
                  } 
                  else {last;}
             }

             # build OH 
               $ratio = $r_bond_C5M_OH/$r_bond_C5M_H;
               $xOH = $xC5M + $ratio*($xH51 - $xC5M);
               $yOH = $yC5M + $ratio*($yH51 - $yC5M);
               $zOH = $zC5M + $ratio*($zH51 - $zC5M);
               
             # build HO       

                 # use procedure applied within DC-REPSWA
                 # define r1=N6, r2=C5M, r3=OH, r, theta, phi
                 $r1x = $xN6; $r1y = $yN6; $r1z = $zN6;
                 $r2x = $xC5M; $r2y = $yC5M; $r2z = $zC5M;
                 $r3x = $xOH; $r3y = $yOH; $r3z = $zOH;

                 $r     = $r_bond_OH_HO;
                 $theta = $pi*(180.0-$theta_C5M_OH_HO)/180.0;
               
                 # HO
                 $phi   = $phi_N6_C5M_OH_HO;
                 # build r4 based on r1, r2, r3, r, theta, phi
                 &build_r4;
                 # read out HO
                 $xHO = $r4x; $yHO = $r4y; $zHO = $r4z;                    

                 $insert_line = $iline + $natm_AME - 1;

	 }#endif-AME#

         $entry = $entry[$iline];
         $entry_new = $entry;

         if ($insert_line == $iline) {
             $atom_new = 'OH  ';
             substr($entry_new,12,4) = $atom_new;
             $x_field = sprintf("%8.5g",$xOH);
             $y_field = sprintf("%8.5g",$yOH);
             $z_field = sprintf("%8.5g",$zOH);
             substr($entry_new,30,8) = $x_field;
             substr($entry_new,38,8) = $y_field;
             substr($entry_new,46,8) = $z_field;
             substr($entry_new,6,5)  = '    0';
             substr($entry_new,70,15) = '   O';
             print "$entry_new\n";

             $atom_new = 'HO  ';
             substr($entry_new,12,4) = $atom_new;
             $x_field = sprintf("%8.5g",$xHO);
             $y_field = sprintf("%8.5g",$yHO);
             $z_field = sprintf("%8.5g",$zHO);
             substr($entry_new,30,8) = $x_field;
             substr($entry_new,38,8) = $y_field;
             substr($entry_new,46,8) = $z_field;
             substr($entry_new,6,5)  = '    0';
             substr($entry_new,70,15) = '   H';
             print "$entry_new\n";
             print "$entry\n";
         }
         else {
            $resName[$iline] = substr($entry[$iline],17,3); 
            $resName[$iline] =~ tr/ //d;
            if ($resName[$iline] eq 'AHM') { 
               if ((substr($entry,12,4) ne 'H51 ') &&
                   (substr($entry,12,4) ne ' H51')) {print "$entry\n";}
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
