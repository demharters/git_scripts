#! /usr/bin/perl -w
# 
# usage 
# 
# ./energy_vs_rmsd.pl rmsd.dat pot_energy.dat inter_energy.dat > energy_vs_rmsd.dat
#
# rmsd.dat : contains the RMSD from native in one column
#            length: # of steps
# pot_energy : contains the pot_energy for each replica
#              in N coulmns (first column is the replica of
#              interest temperature), length: # of steps
#              ! skiping first line title !
# inter_energy : same as pot_energy
#              ! skipping first line title !
# 
# STDOUT :
# 
#    Col1       Col2         Col3
#
#    rmsd    pot_energy   inter_energy
#

  $rmsd_file   = $ARGV[0];
  $epot_file   = $ARGV[1];

  open(RMSD,"$rmsd_file");
  open(EPOT,"$epot_file");

  $ndata_rmsd = 0;
  while (<RMSD>) {
        chomp;
        $ndata_rmsd++;
        @entry = split;
        $rmsd[$ndata_rmsd] = $entry[0]; 
  }

  $ndata_epot = 0;
  while (<EPOT>) {
        chomp;
        $ndata_epot++;
        @entry = split;
        $epot[$ndata_epot-1] = $entry[0]; 
  }
  $ndata_epot;


  if ($ndata_rmsd != $ndata_epot) {die "$rmsd_file and $epot_file has different # of lines"} 
  $ndata = $ndata_epot;

  for ($idata=1;$idata<=$ndata;$idata++) {
      $rmsd = $rmsd[$idata];
      $epot = $epot[$idata];
      print "$rmsd  $epot\n";
  }
