#! /bin/bash
./remove_H.pl $1 >| 1.pdb
./rename_DNA_res.pl 1.pdb >| 2.pdb
./pymol_hadd.py 2.pdb
./rename_DNA_H.pl 2_hadd.pdb >| 4.pdb
./pdb_to_PAS_pdb.pl 4.pdb >| 5.pdb
./removeFirstLast.py 5.pdb
echo "CBLC >AB"|cat - 5_trunc.pdb > /tmp/out && mv /tmp/out init_tmp.pdb
mosaics_epi.x 0step.input > output
mv sampled.pos.pdb init.pdb
sed -i 's/MODEL\ 1/CBLC\ >AB/g' init.pdb
clean

rm 1.pdb 2.pdb 2_hadd.pdb 4.pdb 5.pdb 5_trunc.pdb init_tmp.pdb
