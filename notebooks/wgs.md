This notebook contains the code to reproduce the population structure and geneflow analyses

Subset file so it can be used as input for PCAngsd and NGSadmix

```bash
#Extract header from BEAGLE file
head -n 1 wholegenome_pruned.beagle > header
#Pick 10000 lines at random
shuf -n 10000 wholegenome_pruned.beagle > wholegenome_pruned_1million.beagle
#OR if file is too big for your computer to handle
cat wholegenome_pruned.beagle | awk 'BEGIN {srand()} !/^$/ { if (rand() <= .001) print $0}' > wholegenome_pruned_1K.beagle
cat wholegenome_pruned.beagle | awk 'BEGIN {srand()} !/^$/ { if (rand() <= .01) print $0}' > wholegenome_pruned_10K.beagle
cat wholegenome_pruned.beagle | awk 'BEGIN {srand()} !/^$/ { if (rand() <= .1) print $0}' > wholegenome_pruned_100K.beagle
#Paste header
cat header wholegenome_pruned_1million.beagle > wholegenome_pruned_1million.beagle
#Keep only samples that are common to WGS and candidate genes (samples2keep)
python subset_beagle.py --beagle wholegenome_pruned_1million.beagle --samples samples2keep --out wholegenome_pruned_1million_subset.beagle
#Gzip file because PCAngsd needs it in this format
gzip -k wholegenome_pruned_1million_subset.beagle
```



Run PCAngsd

```bash
#!/bin/sh
PATH_INPUT="/data/Users/Andrea/silvereye/wgs_beagle/raw/"
PATH_PCANGSD="/home/zoo/sjoh4959/sjoh4959/projects/0.0_silvereye-candidate-genes/src/pcangsd/"
PATH_OUT="~/sjoh4959/projects/0.0_silvereye-candidate-genes/analyses/pca_pop_str"

python "${PATH_PCANGSD}pcangsd.py" -beagle "${PATH_INPUT}wholegenome_pruned_1million_subset.beagle.gz" -out "${PATH_OUT}output" -threads 16 &
```



Run NGSAdmix

```bash
#!/bin/sh
PATH_INPUT="/data/Users/Andrea/silvereye/wgs_beagle/raw/"
PATH_NGSADMIX="/home/zoo/sjoh4959/sjoh4959/projects/0.0_silvereye-candidate-genes/src/"
PATH_OUT="~/sjoh4959/projects/0.0_silvereye-candidate-genes/analyses/NGSadmix/"


for i in {1..25};
	do nohup "${PATH_NGSADMIX}NGSadmix" -likes "${PATH_INPUT}wholegenome_pruned_1million_subset.beagle" -K $i -o "${PATH_OUT}ngsadmix_$i" -minMaf 0.05 -seed 1 &
done

#Explore substructuring by splitting the dataset in K1 and K2

#!/bin/sh
for i in {1..15};
	do nohup "${PATH_NGSADMIX}NGSadmix" -likes "${PATH_INPUT}wholegenome_pruned_1million_subset_melanesia.beagle" -K $i -o "${PATH_OUT}ngsadmix_mel_$i" -minMaf 0.05 -seed 1 &
done

#!/bin/sh
for i in {1..10};
	do nohup "${PATH_NGSADMIX}NGSadmix" -likes "${PATH_INPUT}wholegenome_pruned_1million_subset_oz.beagle" -K $i -o "${PATH_OUT}ngsadmix_oz_$i" -minMaf 0.05 -seed 1 &
done
```

Generate files for BA3-SNP analysis

```bash
PATH_INPUT="/data/Users/Andrea/silvereye/wgs_beagle/raw/"
PATH_SCRIPTS="/home/zoo/sjoh4959/sjoh4959/projects/0.0_silvereye-candidate-genes/src/"

#Remove those bases with low confidence in genotype call (i.e. 0.33)

awk '{ for (i=1; i<=NF-2; i++) if ($i == "0.333333" && $(i+1) == "0.333333" && $(i+2) == "0.333333") next } 1' "${PATH_INPUT}wholegenome_pruned.beagle} >  "${PATH_INPUT}wholegenome_pruned_nouncertainSNPs.beagle"

#Convert BEAGLE file to VCF file
python "${PATH_SCRIPTS}my_scripts/beagle_to_vcf.py" "${PATH_INPUT}wholegenome_pruned_nouncertainSNPs.beagle" "${PATH_INPUT}wholegenome_pruned_nouncertainSNPs.vcf"
#Zip file
bgzip "${PATH_INPUT}wholegenome_pruned_nouncertainSNPs.vcf"
#Tabix zipped file
tabix "${PATH_INPUT}wholegenome_pruned_nouncertainSNPs.vcf.gz"

#Convert BCF to BA3 file
"${PATH_SCRIPTS}/others/ugnix/scripts/bcf2ba3" "${PATH_INPUT}wholegenome_pruned_nouncertainSNPs.vcf.gz" > "${PATH_INPUT}wholegenome_77K.ba"
#Update default labels with the population labels 
python "${PATH_SCRIPTS}others/ugnix/scripts/my_scripts/update_labels_BA.py" "${PATH_INPUT}wholegenome_77K.ba" poplabel > "${PATH_INPUT}wholegenome_pruned_77K.ba"

#Filter out those populations that are in the wgs data but not in the candidate gene dataset
cat "${PATH_INPUT}wholegenome_77K.ba" | "${PATH_SCRIPTS}my_scripts/filter_bayesass_file.sh" Chatham Lord_Howe_Island Norfolk_Island New_Zealand Tasmania Heron_Island Mainland > "${PATH_INPUT}aus_wholegenome_77K.ba"
cat "${PATH_INPUT}wholegenome_77K.ba" | "${PATH_SCRIPTS}my_scripts/filter_bayesass_file.sh" Efate Santo Gaua Grand_Terre Lifou Malekula Mare Pentecost Ouvea  Ambrym Ambae Tanna > "${PATH_INPUT}mel_wholegenome_77K.ba"

#create 5 subsets of 10K each
count=5
for i in $(seq $count)
do 
python "${PATH_SCRIPTS}select_snps_ba.py" 77K/mel_wholegenome_77K.ba 10000 mel_wholegenome_10K_${i}.ba
done

#Run autotune
conda activate bayesass

for size in 10000; do
	for i in *ba
do
	cd $size
		nohup python "${PATH_SCRIPTS}others/BA3-SNPS-autotune/BA3-SNPS-autotune.py" -i $i -l $size &
cd ..
	done
done
```

Run BA3-SNP 10 times on each subset

```bash
PATH_BA3="/home/zoo/sjoh4959/sjoh4959/projects/0.0_silvereye-candidate-genes/src/"
PATH_INPUT="/data/Users/Andrea/silvereye/wgs_beagle/raw/"
PATH_OUT="~/sjoh4959/projects/0.0_silvereye-candidate-genes/analyses/NGSadmix/"

for i in {1..10};
        do nohup "${PATH_BA3}BA3-SNPS-Ubuntu64" --file "${PATH_INPUT}mel_wholegenome_10K_1.inp" -t -s $RANDOM -m0.1 -f0.0001 -a0.6 -o "${PATH_OUT}out_mel_$i.txt" -i 1000000 -b 100000 -n 100 -v --loci 10000 &
done

for i in {1..10};
        do nohup "${PATH_BA3}BA3-SNPS-Ubuntu64" --file "${PATH_INPUT}/aus_wholegenome_10K_1.inp" -t -s $RANDOM -m0.5 -f0.7 -a0.7 -o "${PATH_OUT}out_aus_$i.txt" -i 10000000 -b 1000000 -n 1000 -v --loci 10000 &
done
```

