#! /usr/bin/perl -w

  $nfold = 0;
  while (defined($fold = glob("*.pdb"))) {
        $nfold++;
        $fold[$nfold] = $fold;
        $fold_pre = $fold."pre";
        $fold_3pt = $fold."_3pt";
        print "$fold\n";
        system("/home/scratch/mosaics/scripts/extract_3pt_pdb.pl < $fold > $fold_pre");
        system("/home/scratch/mosaics/scripts/correct_3pt_pdb.pl $fold_pre > $fold_3pt");
        system("/home/scratch/mosaics/scripts/score_3pt.pl $fold_3pt > score");
        open(SCORE,"score");
        while (<SCORE>) {$score_=$_;}
        if ($score_ eq "corrupted") {
            system("rm $fold_3pt");
        } 
        system("rm $fold_pre");
  }
