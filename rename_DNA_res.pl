#! /usr/bin/perl -w
#
# This script renames DNA residues, if they are 
# named by DT, DA, DC, DG, instead of T, A, C, G 
#
# usage:
#       ./rename_DNA_res.pl A.pdb > A_.pdb
#

  $pdb_file = $ARGV[0];
  $chainList   = "ABCDEFGHIJKLMNOPQRSTUWXYZ";
  $ind_chainList = 0;            

  $resSeq = 0;

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

             if ($resSeq != $resSeq[$natom]) { 
                 if (substr($resName[$natom],1,1) eq '5') {
                     $ind_chainList++;
                 }
                 $resSeq = $resSeq[$natom];
             }

             if ($chainID[$natom] eq '') {
		 $chainID[$natom] = substr($chainList,$ind_chainList-1,1);
                 $field = sprintf("%1s",$chainID[$natom]);
                 substr($entry[$natom],21,1) = $field;
             }
	    
             if ($resName[$natom] eq 'T3') {
		 $resName[$natom] = 'T'; $resName[$natom] =~ tr/ //d;
                 $field = sprintf("%3s",$resName[$natom]);
                 substr($entry[$natom],17,3) = $field;
	     }
             if ($resName[$natom] eq 'C3') {
		 $resName[$natom] = 'C'; $resName[$natom] =~ tr/ //d;
                 $field = sprintf("%3s",$resName[$natom]);
                 substr($entry[$natom],17,3) = $field;
	     }
             if ($resName[$natom] eq 'A3') {
		 $resName[$natom] = 'A'; $resName[$natom] =~ tr/ //d;
                 $field = sprintf("%3s",$resName[$natom]);
                 substr($entry[$natom],17,3) = $field;
	     }
             if ($resName[$natom] eq 'G3') {
		 $resName[$natom] = 'G'; $resName[$natom] =~ tr/ //d;
                 $field = sprintf("%3s",$resName[$natom]);
                 substr($entry[$natom],17,3) = $field;
	     }

             if ($resName[$natom] eq 'A5') {
		 $resName[$natom] = 'A'; $resName[$natom] =~ tr/ //d;
                 $field = sprintf("%3s",$resName[$natom]);
                 substr($entry[$natom],17,3) = $field;
	     }
             if ($resName[$natom] eq 'G5') {
		 $resName[$natom] = 'G'; $resName[$natom] =~ tr/ //d;
                 $field = sprintf("%3s",$resName[$natom]);
                 substr($entry[$natom],17,3) = $field;
	     }
             if ($resName[$natom] eq 'C5') {
		 $resName[$natom] = 'C'; $resName[$natom] =~ tr/ //d;
                 $field = sprintf("%3s",$resName[$natom]);
                 substr($entry[$natom],17,3) = $field;
	     }
             if ($resName[$natom] eq 'T5') {
		 $resName[$natom] = 'T'; $resName[$natom] =~ tr/ //d;
                 $field = sprintf("%3s",$resName[$natom]);
                 substr($entry[$natom],17,3) = $field;
	     }

             if ($resName[$natom] eq 'DT') {
		 $resName[$natom] = 'T';
                 $field = sprintf("%3s",$resName[$natom]);
                 substr($entry[$natom],17,3) = $field;
	     }

             if ($resName[$natom] eq 'DA') {
		 $resName[$natom] = 'A';
                 $field = sprintf("%3s",$resName[$natom]);
                 substr($entry[$natom],17,3) = $field;
	     }

             if ($resName[$natom] eq 'DC') {
		 $resName[$natom] = 'C';
                 $field = sprintf("%3s",$resName[$natom]);
                 substr($entry[$natom],17,3) = $field;
	     }

             if ($resName[$natom] eq 'DG') {
		 $resName[$natom] = 'G';
                 $field = sprintf("%3s",$resName[$natom]);
                 substr($entry[$natom],17,3) = $field;
	     }

             if ($resName[$natom] eq 'ADE') {
		 $resName[$natom] = 'A';
                 $field = sprintf("%3s",$resName[$natom]);
                 substr($entry[$natom],17,3) = $field;
	     }

             if ($resName[$natom] eq 'GUA') {
		 $resName[$natom] = 'G';
                 $field = sprintf("%3s",$resName[$natom]);
                 substr($entry[$natom],17,3) = $field;
	     }

             if ($resName[$natom] eq 'CYT') {
		 $resName[$natom] = 'C';
                 $field = sprintf("%3s",$resName[$natom]);
                 substr($entry[$natom],17,3) = $field;
	     }

             if ($resName[$natom] eq 'THY') {
		 $resName[$natom] = 'T';
                 $field = sprintf("%3s",$resName[$natom]);
                 substr($entry[$natom],17,3) = $field;
	     }
             
            if ($resName[$natom] eq 'URA') {
		 $resName[$natom] = 'U';
                 $field = sprintf("%3s",$resName[$natom]);
                 substr($entry[$natom],17,3) = $field;
	     }

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
  
# 3) print                          #
 
  for ($iatm=1;$iatm<=$natom;$iatm++) {
      if ($entry[$iatm] ne 'excess') {
          print "$entry[$iatm]\n"; 
      }
  }#endfor#
