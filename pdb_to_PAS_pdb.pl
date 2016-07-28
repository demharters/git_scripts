#! /usr/bin/perl -w 
# 
# This script renames atoms in the pdb in order to create 
# a pdb file that is compatible with the Palo ALto Sampler
#
#
  $resSeq = 1; 
  $nline  = 0;
  $chainList   = "ABCDEFGHIJKLMNOPQRSTUWXYZ";
  $chainoffset = 0;
  $chainpos_prv = 0;
  while (<>) {
      chomp;
      $entry = $_;
      $TER_found  = (substr($entry,0,3) eq 'TER');
      if ($TER_found) {print "$entry\n";}
      $ATOM_found    = (substr($entry,0,4) eq 'ATOM'); 
      $HETATM_found  = (substr($entry,0,6) eq 'HETATM');
      if ( $ATOM_found || $HETATM_found) {
         $nline++;
#NOUSE
#         @list_entry = split;
#         $atom = $list_entry[2];
#NOUSE
         $atom       = substr($entry,12,4); $atom =~ tr/ //d;
         $resName    = substr($entry,17,3); $resName =~ tr/ //d;
         $resSeq_now = substr($entry,22,4); $resSeq_now =~ tr/ //d;
         $chain_id   = substr($entry,21,1); 
         for ($i=1;$i<=length($chainList);$i++) {
	      $ii = $i-1;
              if (substr($chainList,$ii,1) eq $chain_id) {$chainpos = $i;}
         }

         $char_first = substr($atom,0,1);
         $length     = length($atom);         
#DEBUG
# print "$resName $atom\n";
# <STDIN>;
#DEBUG
         $entry_new = $entry;

         $modified  = 0;
# nucleic acids
         if ( ( ($resName eq 'A') || ($resName eq 'G') ||
                ($resName eq 'C') || ($resName eq 'T') ||
                ($resName eq 'U') ) ) {
                if ($atom eq 'OP1')  {
                    $atom_new = ' O1P';
                    substr($entry_new,12,4) = $atom_new;
                    $atom = $atom_new;
	        }
                if ($atom eq 'OR')  {
                    $atom_new = ' O1P';
                    substr($entry_new,12,4) = $atom_new;
                    $atom = $atom_new;
	        }
                if ($atom eq 'OP2')  {
                    $atom_new = ' O2P';
                    substr($entry_new,12,4) = $atom_new;
                    $atom = $atom_new;
	        }
                if ($atom eq 'OL')  {
                    $atom_new = ' O2P';
                    substr($entry_new,12,4) = $atom_new;
                    $atom = $atom_new;
	        }
                if ($atom eq 'O5*')  {
                    $atom_new = ' O5\'';
                    substr($entry_new,12,4) = $atom_new;
                    $atom = $atom_new;
	        }
                if ($atom eq 'C5*')  {
                    $atom_new = ' C5\'';
                    substr($entry_new,12,4) = $atom_new;
                    $atom = $atom_new;
	        }
                if ($atom eq 'C4*')  {
                    $atom_new = ' C4\'';
                    substr($entry_new,12,4) = $atom_new;
                    $atom = $atom_new;
	        }
                if ($atom eq 'O1\'')  {
                    $atom_new = ' O4\'';
                    substr($entry_new,12,4) = $atom_new;
                    $atom = $atom_new;
	        }
                if ($atom eq 'O4*')  {
                    $atom_new = ' O4\'';
                    substr($entry_new,12,4) = $atom_new;
                    $atom = $atom_new;
	        }
                if ($atom eq 'C3*')  {
                    $atom_new = ' C3\'';
                    substr($entry_new,12,4) = $atom_new;
                    $atom = $atom_new;
	        }
                if ($atom eq 'O3*')  {
                    $atom_new = ' O3\'';
                    substr($entry_new,12,4) = $atom_new;
                    $atom = $atom_new;
	        }
                if ($atom eq 'C2*')  {
                    $atom_new = ' C2\'';
                    substr($entry_new,12,4) = $atom_new;
                    $atom = $atom_new;
	        }
                if ($atom eq 'C1*')  {
                    $atom_new = ' C1\'';
                    substr($entry_new,12,4) = $atom_new;
                    $atom = $atom_new;
	        }
                if ($atom eq 'H1*')  {
                    $atom_new = ' H1\'';
                    substr($entry_new,12,4) = $atom_new;
                    $atom = $atom_new;
	        }
                if ($atom eq 'H3*')  {
                    $atom_new = ' H3\'';
                    substr($entry_new,12,4) = $atom_new;
                    $atom = $atom_new;
	        }
                if ($atom eq 'H4*')  {
                    $atom_new = ' H4\'';
                    substr($entry_new,12,4) = $atom_new;
                    $atom = $atom_new;
	        }
                if ($atom eq '1H2*')  {
                    $atom_new = ' H2\'';
                    substr($entry_new,12,4) = $atom_new;
                    $atom = $atom_new;
                    $modified = 1;
	        }
                if ($atom eq '2H2*')  {
                    $atom_new = 'H2\'\'';
                    substr($entry_new,12,4) = $atom_new;
                    $atom = $atom_new;
                    $modified = 1;
	        }
                if ($atom eq '1H5*')  {
                    $atom_new = ' H5\'';
                    substr($entry_new,12,4) = $atom_new;
                    $atom = $atom_new;
                    $modified = 1;
	        }
                if ($atom eq '2H5*')  {
                    $atom_new = 'H5\'\'';
                    substr($entry_new,12,4) = $atom_new;
                    $atom = $atom_new;
                    $modified = 1;
	        }
         }

         if ($resName eq 'A') {
                if ($atom eq '1H6')  {
                    $atom_new = ' H61';
                    substr($entry_new,12,4) = $atom_new;
                    $atom = $atom_new;
                    $modified = 1; 
	        }
                if ($atom eq '2H6')  {
                    $atom_new = ' H62';
                    substr($entry_new,12,4) = $atom_new;
                    $atom = $atom_new;
                    $modified = 1;
	        }
         }


         if ($resName eq 'G') {
                if ($atom eq '1H2')  {
                    $atom_new = ' H21';
                    substr($entry_new,12,4) = $atom_new;
                    $atom = $atom_new;
                    $modified = 1;
	        }
                if ($atom eq '2H2')  {
                    $atom_new = ' H22';
                    substr($entry_new,12,4) = $atom_new;
                    $atom = $atom_new;
                    $modified = 1;
	        }
         }


         if ($resName eq 'T') {
                if ($atom eq '1H5M')  {
                    $atom_new = ' H51';
                    substr($entry_new,12,4) = $atom_new;
                    $atom = $atom_new;
                    $modified = 1;
	        }
                if ($atom eq '2H5M')  {
                    $atom_new = ' H52';
                    substr($entry_new,12,4) = $atom_new;
                    $atom = $atom_new;
                    $modified = 1;
	        }
                if ($atom eq '3H5M')  {
                    $atom_new = ' H53';
                    substr($entry_new,12,4) = $atom_new;
                    $atom = $atom_new;
                    $modified = 1;
	        }
         }


         if ($resName eq 'C') {
                if ($atom eq '1H4')  {
                    $atom_new = ' H41';
                    substr($entry_new,12,4) = $atom_new;
                    $atom = $atom_new;
                    $modified = 1;
	        }
                if ($atom eq '2H4')  {
                    $atom_new = ' H42';
                    substr($entry_new,12,4) = $atom_new;
                    $atom = $atom_new;
                    $modified = 1;
	        }
         }


# amino acids
         if (($char_first =~ /\d/) && ($modified==0)) {
             $atom_new = substr($atom,1,$length-1).$char_first;
             substr($entry_new,12,$length) = $atom_new;
             $atom = $atom_new;
         }
         if ($atom eq 'H') {
             $atom_new = '  HN';
             substr($entry_new,12,4) = $atom_new;
             $atom = $atom_new;
         }
# changing HA3 to HA1
         if (      ($resName eq 'GLY')  ) {
              
              if ($atom eq 'HA3') {
                  $atom_new = ' HA1';
                  substr($entry_new,12,4) = $atom_new;
                  $atom = $atom_new;
              }
         }
# changing HB3 to HB1
         if (      ($resName eq 'ASP') || 
                   ($resName eq 'LEU') ||
                   ($resName eq 'PRO') ||
                   ($resName eq 'ARG') ||
                   ($resName eq 'ASN') ||
                   ($resName eq 'CYS') ||
                   ($resName eq 'GLN') ||
                   ($resName eq 'HSD') ||
                   ($resName eq 'HSE') ||
                   ($resName eq 'HSP') ||
                   ($resName eq 'LYS') ||
                   ($resName eq 'MET') ||
                   ($resName eq 'PHE') ||
                   ($resName eq 'SER') ||
                   ($resName eq 'TRP') ||
                   ($resName eq 'TYR')   ) {
              
              if ($atom eq 'HB3') {
                  $atom_new = ' HB1';
                  substr($entry_new,12,4) = $atom_new;
                  $atom = $atom_new;
              }
         }
# changing HG3 to HG1
         if (      ($resName eq 'ARG') || 
                   ($resName eq 'GLN') ||
                   ($resName eq 'GLU') ||
                   ($resName eq 'LYS') ||
                   ($resName eq 'MET') ||
                   ($resName eq 'PRO')    ) {
              
              if ($atom eq 'HG3') {
                  $atom_new = ' HG1';
                  substr($entry_new,12,4) = $atom_new;
                  $atom = $atom_new;
              }
         }
# changing HD3 to HD1
         if (      ($resName eq 'ARG') || 
                   ($resName eq 'LYS') ||
                   ($resName eq 'PRO')    ) {
              
              if ($atom eq 'HD3') {
                  $atom_new = ' HD1';
                  substr($entry_new,12,4) = $atom_new;
                  $atom = $atom_new;
              }
         }
# changing HE3 to HE1
         if (      ($resName eq 'LYS')  ) {
              
              if ($atom eq 'HE3') {
                  $atom_new = ' HE1';
                  substr($entry_new,12,4) = $atom_new;
                  $atom = $atom_new;
              }
         }



         if (($resName eq 'ILE') && ($atom eq 'CD1')) {
             $atom_new = '  CD';
             substr($entry_new,12,4) = $atom_new;
             $atom = $atom_new;
         }
         if (($resName eq 'ILE') && ($atom eq 'HG13')) {
             $atom_new = 'HG11';
             substr($entry_new,12,4) = $atom_new;
             $atom = $atom_new;
         }
         if (($resName eq 'ILE') && ($atom eq 'HD11')) {
             $atom_new = ' HD1';
             substr($entry_new,12,4) = $atom_new;
             $atom = $atom_new;
         }
         if (($resName eq 'ILE') && ($atom eq 'HD12')) {
             $atom_new = ' HD2';
             substr($entry_new,12,4) = $atom_new;
             $atom = $atom_new;
         }
         if (($resName eq 'ILE') && ($atom eq 'HD13')) {
             $atom_new = ' HD3';
             substr($entry_new,12,4) = $atom_new;
             $atom = $atom_new;
         }
         if (($resName eq 'SER') && ($atom eq 'HG')) {
             $atom_new = ' HG1';
             substr($entry_new,12,4) = $atom_new;
             $atom = $atom_new;
         }
         if (($resName eq 'CYS') && (($atom eq 'HG') || ($atom eq 'H'))) {
             $atom_new = ' HG1';
             substr($entry_new,12,4) = $atom_new;
             $atom = $atom_new;
         }


         if ($nline==1) {$resSeq_prv = $resSeq_now;}
         if ($resSeq_prv ne $resSeq_now) {$resSeq++;}
         $resSeq_diff = $resSeq_now - $resSeq_prv;
         $resSeq_prv  = $resSeq_now;
         $resSeq_field = sprintf("%4s",$resSeq); 
         substr($entry_new,22,4) = $resSeq_field;

         $chainpos_diff = $chainpos - $chainpos_prv;
         if ( ($resSeq_diff > 1) && 
              ($chain_id ne ' ') &&
              ($chainpos_diff == 0) ) {
               $chainoffset++;
#DEBUG
#               print "$chainoffset\n";
#               <STDIN>;
#DEBUG
         }
         $chainpos     += $chainoffset;
         $chain_id_new  = substr($chainList,$chainpos-1,1);
         $chainpos_prv  = $chainpos;
#DEBUG
#         if ($resSeq_diff > 1) {
#             print "chainpos: $chainpos\n";
#             print "new molecule\n";
#             print "chainoffset: $chainoffset, chain_id_new: $chain_id_new\n";
#             <STDIN>;
#         }
#DEBUG
         if ($ATOM_found==1) {
             substr($entry_new,21,1) = $chain_id_new;
         }
         print "$entry_new\n";
#DEBUG
#         print "$entry\n";
#         print "$atom $atom_new $length $char_first\n";
#         <STDIN>;
#DEBUG
     }#ATOM,HETATM#
  }#while (<>)#
  print "END\n"
