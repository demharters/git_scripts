#! /usr/bin/perl -w

   $file = $ARGV[0];
   system("/home/minary/Perl_Scripts/MOSAICS_scripts/pdb_handlers/pdb_to_PAS_pdb.pl $file |
           /home/minary/Perl_Scripts/MOSAICS_scripts/pdb_handlers/protonate_CYS.pl |
           /home/minary/Perl_Scripts/MOSAICS_scripts/pdb_handlers/protonate_HIS.pl");
