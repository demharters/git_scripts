#! /usr/bin/perl -w 
# 
# This script removes all H atoms from a pdb file 
#
#
  while (<>) {
      chomp;
      $entry = $_;
      if (substr($entry,0,4) eq 'ATOM') {
         $valid = 1;
         $atom       = substr($entry,12,4); $atom =~ tr/ //d;

         $char_first = substr($atom,0,1);
         $char_secon = substr($atom,1,1);

         if (($char_first =~ /\d/) && ($char_secon eq 'H')) {
             $valid = 0;
         }
         if ($char_first eq 'H') {
             $valid = 0;
         }
         if ($valid == 1) {
             print "$entry\n";
         }
      }
      else {
        print "$entry\n";
      }
 }
 print "END\n"
