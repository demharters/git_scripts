#! /usr/bin/perl -w 
# 
# This script protonates all CYS residues that are missing atom HG1
#
#     |
#  HN-N
#     |   HB1  
#     |   | 
#  HA-CA--CB--SG
#     |   |    \
#     |   HB2   HG1
#   O=C
#     |
#
# Equilibrium values (charmm 22)
# Bond(SG-HG1)       =   Bond(S-HS)         =  1.325 A   
# Bend(CB-SG-HG1)    =   Bend(CT2-S-HS)     =  95 deg 
# Tors(CA-CB-SG-HG1) =   Tors(CT1-CT2-S-HS) =  180.0 deg     
#
#
#
  $r_bond_SG_HG1    =  1.325;    # in A
  $theta_CB_SG_HG1  =  95.0;     # in deg
  $phi_CA_CB_SG_HG1 =  180.0;    # in deg
  $pi               =  3.14;     # in rad


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

  $resSeq_CYS  = -100;
  $insert_line = 0;
  for ($iline=1;$iline<=$nline;$iline++) {
         $entry = $entry[$iline];
         $entry_new = $entry;


        if ( (substr($entry,0,3) eq 'TER') || 
             (substr($entry,0,3) eq 'END')  ) {
             print "$entry\n";
             next;
        } 


         if ( ($resName[$iline] eq 'CYS') &&  
              ($resSeq_CYS != $resSeq[$iline]) ) {
             $resSeq_CYS = $resSeq[$iline];
             $natm_CYS  = 0;
             $HG1_found = 0;
             for ($jline=$iline;$jline<=$nline;$jline++) {
                  if ($resSeq[$jline] eq  $resSeq_CYS) {
                      $natm_CYS++;
                      if ($atom[$jline] eq 'HG1') {$HG1_found=1;}
                      if ($atom[$jline] eq 'CA') {
                          $xCA = $x[$jline];
                          $yCA = $y[$jline];
                          $zCA = $z[$jline];
                      }
                      if ($atom[$jline] eq 'CB') {
                          $xCB = $x[$jline];
                          $yCB = $y[$jline];
                          $zCB = $z[$jline];
                      }
                      if ($atom[$jline] eq 'SG') {
                          $xSG = $x[$jline];
                          $ySG = $y[$jline];
                          $zSG = $z[$jline];
                      }
                  } 
                  else {last;}
             }
             
             if ($HG1_found==0) {
#DEBUG
#                 print "natm_CYS: $natm_CYS\n";
#DEBUG
                 # use procedure applied within DC-REPSWA

                 # define r1=CA, r2=CB, r3=SG, r, theta, phi
                 $r1x = $xCA; $r1y = $yCA; $r1z = $zCA;
                 $r2x = $xCB; $r2y = $yCB; $r2z = $zCB;
                 $r3x = $xSG; $r3y = $ySG; $r3z = $zSG;
                 
                 $r     = $r_bond_SG_HG1;
                 $theta = $pi*(180.0-$theta_CB_SG_HG1)/180.0;
                 $phi   = $pi*$phi_CA_CB_SG_HG1/180.0;

                 # build r4 based on r1, r2, r3, r, theta, phi
                 &build_r4;

                 # read out HG1
                 $xHG1 = $r4x; $yHG1 = $r4y; $zHG1 = $r4z;                    

                 $insert_line = $iline + $natm_CYS - 1;
             }
	 }

         if ($insert_line == $iline) {
             $atom_new = 'HG1 ';
             substr($entry_new,12,4) = $atom_new;
#DEBUG
#	      print " xHG1 : $xHG1\n";
#             print " yHG1 : $yHG1\n";
#             print " zHG1 : $zHG1\n";
#DEBUG
             $x_field = sprintf("%8.5g",$xHG1);
             $y_field = sprintf("%8.5g",$yHG1);
             $z_field = sprintf("%8.5g",$zHG1);
             substr($entry_new,30,8) = $x_field;
             substr($entry_new,38,8) = $y_field;
             substr($entry_new,46,8) = $z_field;
             substr($entry_new,6,5)  = '    0';

             print "$entry\n";
             print "$entry_new\n";
         }
         else { 
             print "$entry\n";
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
