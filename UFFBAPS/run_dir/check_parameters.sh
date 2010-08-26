#! /bin/tcsh -f
set root = /rnase_local/chresten/UFFBAPS/run_dir/

#echo $root
echo energy.cfg
diff $root energy.cfg
echo es.cfg
diff $root es.cfg
echo hbonds.cfg
diff $root hbonds.cfg
echo hp.cfg
diff $root hp.cfg
echo ligand_entropy.cfg
diff $root ligand_entropy.cfg
echo msc.cfg
diff $root msc.cfg
echo protein_entropy.cfg
diff $root protein_entropy.cfg
echo vdw_line.cfg
diff $root vdw_line.cfg
