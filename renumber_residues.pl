#! /usr/bin/perl -w 
# 
# This script renumbers residues in a way present in a template pdb
# usage:
#
# ./renumber_residues.pl template.pdb target.pdb > target_renumbered.pdb 
#
#
 
  $template_file = $ARGV[0];
  $target_file   = $ARGV[1];

  open(TEMPLATE,"$template_file");
  open(TARGET,"$target_file");
 
  $resSeq = 1; 
  $nline  = 0;
  while (<TEMPLATE>) {
      chomp;
      $entry = $_;
      $ATOM_found    = (substr($entry,0,4) eq 'ATOM'); 
      $HETATM_found  = (substr($entry,0,6) eq 'HETATM');
      if ( $ATOM_found || $HETATM_found) {
         $nline++;
         $atom       = substr($entry,12,4); $atom =~ tr/ //d;
         $resName    = substr($entry,17,3); $resName =~ tr/ //d;
         $resSeq_now = substr($entry,22,4); $resSeq_now =~ tr/ //d;
         $entry_new = $entry;
         $entry[$nline]  = $entry;
         $resSeq[$nline] = $resSeq_now; 
      }#endif#
  }#endwhile#

  $nline = 0;
  while (<TARGET>) {
      chomp;
      $entry = $_;
      $ATOM_found    = (substr($entry,0,4) eq 'ATOM'); 
      $HETATM_found  = (substr($entry,0,6) eq 'HETATM');
      if ( $ATOM_found || $HETATM_found) {
           $nline++;
           $entry_new = $entry;
           $resSeq_field = sprintf("%4s",$resSeq[$nline]);
           substr($entry_new,22,4) = $resSeq_field;
           print "$entry_new\n";
      } else {
           print "$entry\n";
      }#endif#
  }#endwhile#  
