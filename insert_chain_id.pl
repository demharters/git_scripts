#! /usr/bin/perl -w
#
# This script inserts a chain ID to a pdb between resSeq_min and resSeq_max 
#
# usage:
#       ./insert_chain_id.pl A.pdb ID resSeq_min resSeq_max > A_.pdb
#

  $pdb_file = $ARGV[0];
  $id       = $ARGV[1];
  $resSeq_min = $ARGV[2];
  $resSEq_max = $ARGV[3]; 
  open(PDB,"$pdb_file");

# read in .pdb

  $natom = 0;
  while (<PDB>) {
      chomp;
      $entry = $_;
      if (substr($entry,0,4) eq 'ATOM') {
         $altLoc = substr($entry,16,1);
         if ($altLoc eq ' ' || $altLoc eq 'A') {
             $natom++;
             $resSeq[$natom] = substr($entry,22,4); $resSeq[$natom] =~ tr/ //d;
             $entry[$natom] = $entry;
             if ( ($resSeq[$natom] >= $resSeq_min) && 
                  ($resSeq[$natom] <= $resSEq_max)  )  {
                   substr($entry[$natom],21,1) = $id;
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
