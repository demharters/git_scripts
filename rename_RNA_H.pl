#! /usr/bin/perl -w
#
# This script renames H-atoms added to RNA, so that it is 
# compatible with amber99 topology
# The H atoms can be added using the pymol package for 
# example
#
# usage:
#       ./rename_RNA_H.pl A.pdb > A_.pdb
#
# A.pdb  :  pdb having some RNA with added H's
# A_.pdb :  RNA H's are renamed
#

  $pdb_file = $ARGV[0];

  open(PDB,"$pdb_file");

# read in .pdb

  $natom = 0;
  while (<PDB>) {
      chomp;
      $entry = $_;
      if (substr($entry,0,4) eq 'ATOM') {
         @list_entry = split;
         $atom = $list_entry[2];
         $altLoc = substr($entry,16,1);
         if ($altLoc eq ' ' || $altLoc eq 'A') {
             $natom++;
             $atom[$natom]       = substr($entry,12,4); $atom[$natom]    =~ tr/ //d;
             $resName[$natom]    = substr($entry,17,3); $resName[$natom] =~ tr/ //d;
             $chainID[$natom]    = substr($entry,21,1); $chainID[$natom] =~ tr/ //d;
             $resSeq[$natom]     = substr($entry,22,4); $resSeq[$natom]  =~ tr/ //d;
# NOUSE
#             $char_first[$natom] = substr($atom,0,1);
# NOUSE
             $x[$natom]     = substr($entry,30,8); $x[$natom] =~ tr/ //d;
             $y[$natom]     = substr($entry,38,8); $y[$natom] =~ tr/ //d;
             $z[$natom]     = substr($entry,46,8); $z[$natom] =~ tr/ //d;
             $entry[$natom] = $entry;
         }
#DEBUG
# print "$entry\n";
#DEBUG
      }
#NOUSE
#      if (substr($entry,0,3) eq 'END' || substr($entry,0,3) eq 'TER') {
#          last;
#      }
#NOUSE
  }

 
# 1) separate atoms into residues                           #

  $nres   = 1;
  $resSeq = $resSeq[1];
  $natom_res[$nres] = 0;
  for ($iatom=1;$iatom<=$natom;$iatom++) {
       if ($resSeq[$iatom] ne $resSeq) {
           $resSeq = $resSeq[$iatom];
           $nres++;
           $natom_res[$nres]++;
           $ind_atom[$nres][$natom_res[$nres]] = $iatom;  
       } else {
           $natom_res[$nres]++;
           $ind_atom[$nres][$natom_res[$nres]] = $iatom;  
       }
  }

# 2) rename H atoms of RNA residues                         #

  for ($ires=1;$ires<=$nres;$ires++) {
       $fst_atm  = $ind_atom[$ires][1];
       $res_Name = $resName[$fst_atm];

#DEBUG
#                      print("$res_Name\n");
#		      print("$fst_atm\n");
#		       print("Htype: $typ, Atom_type $typ_heavy\n");
#                       <STDIN>;
#DEBUG

       #GUA______________________________________________________
       if ( ($res_Name eq 'G')   ||
            ($res_Name eq 'GUA') ||
            ($res_Name eq 'DG') ) {
              
             $H21_found = 0;
             $H5p_found = 0;
             $H2p_found = 0;
	     for ($iatom=1;$iatom<=$natom_res[$ires];$iatom++) {
		  $iatm    = $ind_atom[$ires][$iatom];
	 	  $typ     = $atom[$iatm];
                  $fch_typ = substr($typ,0,1);
                  if ($fch_typ eq 'H') {
		       $dmin = 100.0;
                       for ($jatom=1;$jatom<=$natom_res[$ires];$jatom++) {
			    $jatm    = $ind_atom[$ires][$jatom];
                            $typj     = $atom[$jatm];
                            $fch_typj = substr($typj,0,1);
                            if ($fch_typj ne 'H') {
                                $dx = $x[$iatm] - $x[$jatm];        
                                $dy = $y[$iatm] - $y[$jatm];
                                $dz = $z[$iatm] - $z[$jatm];
                                $r = sqrt($dx*$dx+$dy*$dy+$dz*$dz);
                                if ($r < $dmin) {
                                    $dmin = $r;
                                    $ind_heavy = $jatm;
                                }#endif#
			   }#endif#
                       }#endif#
                       $typ_heavy = $atom[$ind_heavy];

#DEBUG
#                      print("$res_Name\n");
#		      print("$entry[$iatm]\n");
#		       print("Htype: $typ, Atom_type $typ_heavy\n");
#                       <STDIN>;
#DEBUG

                       if ($typ_heavy eq 'C8') {
                           $atom[$iatm] = 'H8';
                           $field = sprintf("%4s",$atom[$iatm]);
                           substr($entry[$iatm],12,4) = $field;
                       }
                       elsif ($typ_heavy eq 'C4\'') {
                           $atom[$iatm] = 'H4\'';
                           $field = sprintf("%4s",$atom[$iatm]);
                           substr($entry[$iatm],12,4) = $field;
                       }
                       elsif ($typ_heavy eq 'C1\'') {
                           $atom[$iatm] = 'H1\'';
                           $field = sprintf("%4s",$atom[$iatm]);
                           substr($entry[$iatm],12,4) = $field;
                       }
                       elsif ($typ_heavy eq 'N1') {
                           $atom[$iatm] = 'H1';
                           $field = sprintf("%4s",$atom[$iatm]);
                           substr($entry[$iatm],12,4) = $field;
                       }
                       elsif ($typ_heavy eq 'C3\'') {
                           $atom[$iatm] = 'H3\'';
                           $field = sprintf("%4s",$atom[$iatm]);
                           substr($entry[$iatm],12,4) = $field;
                       }
                       elsif ($typ_heavy eq 'N2') {
                           if ($H21_found==0) { 
                               $atom[$iatm] = 'H21';
                               $field = sprintf("%4s",$atom[$iatm]);
                               substr($entry[$iatm],12,4) = $field;
                               $H21_found = 1;
		           } else  { 
                               $atom[$iatm] = 'H22';
                               $field = sprintf("%4s",$atom[$iatm]);
                               substr($entry[$iatm],12,4) = $field;
		           }#endif#
		       }
                       elsif ($typ_heavy eq 'C2\'') {
                               $atom[$iatm] = 'H2\'\'';
                               $field = sprintf("%4s",$atom[$iatm]);
                               substr($entry[$iatm],12,4) = $field;
		       }
                       elsif ($typ_heavy eq 'O2\'') {
                               $atom[$iatm] = 'H2\'';
                               $field = sprintf("%4s",$atom[$iatm]);
                               substr($entry[$iatm],12,4) = $field;
                       }
                       elsif ($typ_heavy eq 'C5\'') {
                           if ($H5p_found==0) { 
                               $atom[$iatm] = 'H5\'';
                               $field = sprintf("%4s",$atom[$iatm]);
                               substr($entry[$iatm],12,4) = $field;
                               $H5p_found = 1;
		           } else  { 
                               $atom[$iatm] = 'H5\'\'';
                               $field = sprintf("%4s",$atom[$iatm]);
                               substr($entry[$iatm],12,4) = $field;
		           }#endif#
		       }
                       else {
			   $entry[$iatm] = 'excess';
                       }#endif#
                  }#endif#
             }#endif#
       }#endif#
       #GUA__________________________________________________END


       #ADE______________________________________________________
       if ( ($res_Name eq 'A')   ||
            ($res_Name eq 'ADE') ||
            ($res_Name eq 'DA') ) {
              
             $H61_found = 0;
             $H5p_found = 0;
             $H2p_found = 0;
	     for ($iatom=1;$iatom<=$natom_res[$ires];$iatom++) {
		  $iatm    = $ind_atom[$ires][$iatom];
	 	  $typ     = $atom[$iatm];
                  $fch_typ = substr($typ,0,1);
                  if ($fch_typ eq 'H') {
		       $dmin = 100.0;
                       for ($jatom=1;$jatom<=$natom_res[$ires];$jatom++) {
			    $jatm    = $ind_atom[$ires][$jatom];
                            $typj     = $atom[$jatm];
                            $fch_typj = substr($typj,0,1);
                            if ($fch_typj ne 'H') {
                                $dx = $x[$iatm] - $x[$jatm];        
                                $dy = $y[$iatm] - $y[$jatm];
                                $dz = $z[$iatm] - $z[$jatm];
                                $r = sqrt($dx*$dx+$dy*$dy+$dz*$dz);
                                if ($r < $dmin) {
                                    $dmin = $r;
                                    $ind_heavy = $jatm;
                                }#endif#
			   }#endif#
                       }#endif#
                       $typ_heavy = $atom[$ind_heavy];

#DEBUG
#                      print("$resName\n");
#		       print("Htype: $typ, Atom_type $typ_heavy\n");
#                       <STDIN>;
#DEBUG

                       if ($typ_heavy eq 'C8') {
                           $atom[$iatm] = 'H8';
                           $field = sprintf("%4s",$atom[$iatm]);
                           substr($entry[$iatm],12,4) = $field;
                       }
                       elsif ($typ_heavy eq 'C4\'') {
                           $atom[$iatm] = 'H4\'';
                           $field = sprintf("%4s",$atom[$iatm]);
                           substr($entry[$iatm],12,4) = $field;
                       }
                       elsif ($typ_heavy eq 'C1\'') {
                           $atom[$iatm] = 'H1\'';
                           $field = sprintf("%4s",$atom[$iatm]);
                           substr($entry[$iatm],12,4) = $field;
                       }
                       elsif ($typ_heavy eq 'C2') {
                           $atom[$iatm] = 'H2';
                           $field = sprintf("%4s",$atom[$iatm]);
                           substr($entry[$iatm],12,4) = $field;
                       }
                       elsif ($typ_heavy eq 'C3\'') {
                           $atom[$iatm] = 'H3\'';
                           $field = sprintf("%4s",$atom[$iatm]);
                           substr($entry[$iatm],12,4) = $field;
                       }
                       elsif ($typ_heavy eq 'N6') {
                           if ($H61_found==0) { 
                               $atom[$iatm] = 'H61';
                               $field = sprintf("%4s",$atom[$iatm]);
                               substr($entry[$iatm],12,4) = $field;
                               $H61_found = 1;
		           } else  { 
                               $atom[$iatm] = 'H62';
                               $field = sprintf("%4s",$atom[$iatm]);
                               substr($entry[$iatm],12,4) = $field;
		           }#endif#
		       }
                       elsif ($typ_heavy eq 'C2\'') {
                               $atom[$iatm] = 'H2\'\'';
                               $field = sprintf("%4s",$atom[$iatm]);
                               substr($entry[$iatm],12,4) = $field;
		       }
                       elsif ($typ_heavy eq 'O2\'') {
                               $atom[$iatm] = 'H2\'';
                               $field = sprintf("%4s",$atom[$iatm]);
                               substr($entry[$iatm],12,4) = $field;
                       }
                       elsif ($typ_heavy eq 'C5\'') {
                           if ($H5p_found==0) { 
                               $atom[$iatm] = 'H5\'';
                               $field = sprintf("%4s",$atom[$iatm]);
                               substr($entry[$iatm],12,4) = $field;
                               $H5p_found = 1;
		           } else  { 
                               $atom[$iatm] = 'H5\'\'';
                               $field = sprintf("%4s",$atom[$iatm]);
                               substr($entry[$iatm],12,4) = $field;
		           }#endif#
		       }
                       else {
			   $entry[$iatm] = 'excess';
                       }#endif#
                  }#endif#
             }#endif#
       }#endif#
       #ADE__________________________________________________END



       #CYT______________________________________________________
       if ( ($res_Name eq 'C')   ||
            ($res_Name eq 'CYT') ||
            ($res_Name eq 'DC') ) {
              
             $H41_found = 0;
             $H5p_found = 0;
             $H2p_found = 0;
	     for ($iatom=1;$iatom<=$natom_res[$ires];$iatom++) {
		  $iatm    = $ind_atom[$ires][$iatom];
	 	  $typ     = $atom[$iatm];
                  $fch_typ = substr($typ,0,1);
                  if ($fch_typ eq 'H') {
		       $dmin = 100.0;
                       for ($jatom=1;$jatom<=$natom_res[$ires];$jatom++) {
			    $jatm     = $ind_atom[$ires][$jatom];
                            $typj     = $atom[$jatm];
                            $fch_typj = substr($typj,0,1);
                            if ($fch_typj ne 'H') {
                                $dx = $x[$iatm] - $x[$jatm];        
                                $dy = $y[$iatm] - $y[$jatm];
                                $dz = $z[$iatm] - $z[$jatm];
                                $r = sqrt($dx*$dx+$dy*$dy+$dz*$dz);
                                if ($r < $dmin) {
                                    $dmin = $r;
                                    $ind_heavy = $jatm;
                                }#endif#
			   }#endif#
                       }#endif#
                       $typ_heavy = $atom[$ind_heavy];

#DEBUG
#                      print("$resName\n");
#		       print("Htype: $typ, Atom_type $typ_heavy\n");
#                       <STDIN>;
#DEBUG

                       if ($typ_heavy eq 'C5') {
                           $atom[$iatm] = 'H5';
                           $field = sprintf("%4s",$atom[$iatm]);
                           substr($entry[$iatm],12,4) = $field;
                       }
                       elsif ($typ_heavy eq 'C6') {
                           $atom[$iatm] = 'H6';
                           $field = sprintf("%4s",$atom[$iatm]);
                           substr($entry[$iatm],12,4) = $field;
                       }
                       elsif ($typ_heavy eq 'C4\'') {
                           $atom[$iatm] = 'H4\'';
                           $field = sprintf("%4s",$atom[$iatm]);
                           substr($entry[$iatm],12,4) = $field;
                       }
                       elsif ($typ_heavy eq 'C1\'') {
                           $atom[$iatm] = 'H1\'';
                           $field = sprintf("%4s",$atom[$iatm]);
                           substr($entry[$iatm],12,4) = $field;
                       }
                       elsif ($typ_heavy eq 'C2') {
                           $atom[$iatm] = 'H2';
                           $field = sprintf("%4s",$atom[$iatm]);
                           substr($entry[$iatm],12,4) = $field;
                       }
                       elsif ($typ_heavy eq 'C3\'') {
                           $atom[$iatm] = 'H3\'';
                           $field = sprintf("%4s",$atom[$iatm]);
                           substr($entry[$iatm],12,4) = $field;
                       }
                       elsif ($typ_heavy eq 'N4') {
                           if ($H41_found==0) { 
                               $atom[$iatm] = 'H41';
                               $field = sprintf("%4s",$atom[$iatm]);
                               substr($entry[$iatm],12,4) = $field;
                               $H41_found = 1;
		           } else  { 
                               $atom[$iatm] = 'H42';
                               $field = sprintf("%4s",$atom[$iatm]);
                               substr($entry[$iatm],12,4) = $field;
		           }#endif#
		       }
                       elsif ($typ_heavy eq 'C2\'') {
                               $atom[$iatm] = 'H2\'\'';
                               $field = sprintf("%4s",$atom[$iatm]);
                               substr($entry[$iatm],12,4) = $field;
		       }
                       elsif ($typ_heavy eq 'O2\'') {
                               $atom[$iatm] = 'H2\'';
                               $field = sprintf("%4s",$atom[$iatm]);
                               substr($entry[$iatm],12,4) = $field;
                       }
                       elsif ($typ_heavy eq 'C5\'') {
                           if ($H5p_found==0) { 
                               $atom[$iatm] = 'H5\'';
                               $field = sprintf("%4s",$atom[$iatm]);
                               substr($entry[$iatm],12,4) = $field;
                               $H5p_found = 1;
		           } else  { 
                               $atom[$iatm] = 'H5\'\'';
                               $field = sprintf("%4s",$atom[$iatm]);
                               substr($entry[$iatm],12,4) = $field;
		           }#endif#
		       }
                       else {
			   $entry[$iatm] = 'excess';
                       }#endif#
                  }#endif#
             }#endif#
       }#endif#
       #CYT__________________________________________________END



       #URA______________________________________________________
       if ( ($res_Name eq 'U')   ||
            ($res_Name eq 'URA') ||
            ($res_Name eq 'DU') ) {
              
             $H51_found = 0;
             $H52_found = 0;
             $H5p_found = 0;
             $H2p_found = 0;
	     for ($iatom=1;$iatom<=$natom_res[$ires];$iatom++) {
		  $iatm    = $ind_atom[$ires][$iatom];
	 	  $typ     = $atom[$iatm];
                  $fch_typ = substr($typ,0,1);
                  if ($fch_typ eq 'H') {
		       $dmin = 100.0;
                       for ($jatom=1;$jatom<=$natom_res[$ires];$jatom++) {
                            $jatm    = $ind_atom[$ires][$jatom];
                            $typj     = $atom[$jatm];
                            $fch_typj = substr($typj,0,1);
                            if ($fch_typj ne 'H') {
                                $dx = $x[$iatm] - $x[$jatm];        
                                $dy = $y[$iatm] - $y[$jatm];
                                $dz = $z[$iatm] - $z[$jatm];
                                $r = sqrt($dx*$dx+$dy*$dy+$dz*$dz);
                                if ($r < $dmin) {
                                    $dmin = $r;
                                    $ind_heavy = $jatm;
                                }#endif#
			   }#endif#
                       }#endif#
                       $typ_heavy = $atom[$ind_heavy];

#DEBUG
#                      print("$resName\n");
#		       print("Htype: $typ, Atom_type $typ_heavy\n");
#                       <STDIN>;
#DEBUG

                       if ($typ_heavy eq 'N3') {
                           $atom[$iatm] = 'H3';
                           $field = sprintf("%4s",$atom[$iatm]);
                           substr($entry[$iatm],12,4) = $field;
                       }
                       elsif ($typ_heavy eq 'C6') {
                           $atom[$iatm] = 'H6';
                           $field = sprintf("%4s",$atom[$iatm]);
                           substr($entry[$iatm],12,4) = $field;
                       }
                       elsif ($typ_heavy eq 'C4\'') {
                           $atom[$iatm] = 'H4\'';
                           $field = sprintf("%4s",$atom[$iatm]);
                           substr($entry[$iatm],12,4) = $field;
                       }
                       elsif ($typ_heavy eq 'C1\'') {
                           $atom[$iatm] = 'H1\'';
                           $field = sprintf("%4s",$atom[$iatm]);
                           substr($entry[$iatm],12,4) = $field;
                       }
                       elsif ($typ_heavy eq 'C2') {
                           $atom[$iatm] = 'H2';
                           $field = sprintf("%4s",$atom[$iatm]);
                           substr($entry[$iatm],12,4) = $field;
                       }
                       elsif ($typ_heavy eq 'C3\'') {
                           $atom[$iatm] = 'H3\'';
                           $field = sprintf("%4s",$atom[$iatm]);
                           substr($entry[$iatm],12,4) = $field;
                       }
                       elsif ($typ_heavy eq 'C5') {
                           $atom[$iatm] = 'H5';
                           $field = sprintf("%4s",$atom[$iatm]);
                           substr($entry[$iatm],12,4) = $field;
                       }
                       elsif ( ($typ_heavy eq 'C5M') ||
                               ($typ_heavy eq 'C7')  ) {
                           if ($H51_found==0) {
                               $atom[$iatm] = 'H51';
                               $field = sprintf("%4s",$atom[$iatm]);
                               substr($entry[$iatm],12,4) = $field;
                               $atom[$ind_heavy] = 'C5M';
                               $field = sprintf("%4s",$atom[$ind_heavy]);
                               substr($entry[$ind_heavy],12,4) = $field;
                               $H51_found = 1;
		           } else  { 
                              if ($H52_found==0) { 
                                $atom[$iatm] = 'H52';
                                $field = sprintf("%4s",$atom[$iatm]);
                                substr($entry[$iatm],12,4) = $field;
                                $H52_found = 1;
			      } else {
                                $atom[$iatm] = 'H53';
                                $field = sprintf("%4s",$atom[$iatm]);
                                substr($entry[$iatm],12,4) = $field;
			      }#endif#
		           }#endif#
		       }
                       elsif ($typ_heavy eq 'C2\'') {
                               $atom[$iatm] = 'H2\'\'';
                               $field = sprintf("%4s",$atom[$iatm]);
                               substr($entry[$iatm],12,4) = $field;
		       }
                       elsif ($typ_heavy eq 'O2\'') {
                               $atom[$iatm] = 'H2\'';
                               $field = sprintf("%4s",$atom[$iatm]);
                               substr($entry[$iatm],12,4) = $field;
		       }
                       elsif ($typ_heavy eq 'C5\'') {
                           if ($H5p_found==0) { 
                               $atom[$iatm] = 'H5\'';
                               $field = sprintf("%4s",$atom[$iatm]);
                               substr($entry[$iatm],12,4) = $field;
                               $H5p_found = 1;
		           } else  { 
                               $atom[$iatm] = 'H5\'\'';
                               $field = sprintf("%4s",$atom[$iatm]);
                               substr($entry[$iatm],12,4) = $field;
		           }#endif#
		       }
                       else {
			   $entry[$iatm] = 'excess';
                       }#endif#
                  }#endif#
             }#endif#
       }#endif#
       #URA__________________________________________________END


  }#endfor#
   
  
# 3) print                          #
 
  for ($iatm=1;$iatm<=$natom;$iatm++) {
      if ($entry[$iatm] ne 'excess') {
          print "$entry[$iatm]\n"; 
      }
  }#endfor#
