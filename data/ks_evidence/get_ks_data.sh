wget https://research.bioinformatics.udel.edu/iptmnet_data/files/current/ptm.txt
grep 'PHOSPHORYLATION' ptm.txt | grep 'Homo sapiens (Human)' > iPTMnet_ptm.txt

awk -F'\t' '{if ($7 != "") print $3"\t"$4"\t"$6"\t"$7"\t"$8}' iPTMnet_ptm.txt > iPTMnet_ks.tsv

sed -i '1i substrate_acc\tsubstrate_gene\tsite\tkinase_acc\tkinase_gene' iPTMnet_ks.tsv

# base URL for substrate - https://research.bioinformatics.udel.edu/iptmnet/entry/{substrate_uniprot_acc}/#asSub

rm iPTMnet_ptm.txt
rm ptm.txt
