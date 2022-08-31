This notebook contains the code to reproduce the population structure and geneflow analyses

Subset file so it can be used as input for PCAngsd and NGSadmix

```bash
#Extract header from BEAGLE file
head -n 1 wholegenome_pruned.beagle > header
#Pick 10000 lines at random
shuf -n 10000 wholegenome_pruned.beagle > wholegenome_pruned_1million.beagle
#OR if file is too big for your computer to handle
cat wholegenome_pruned.beagle | awk 'BEGIN {srand()} !/^$/ { if (rand() <= .01) print $0}' > wholegenome_pruned_1million.beagle
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



Run BA3-SNP

```bash
PATH_BA3="/home/zoo/sjoh4959/sjoh4959/projects/0.0_silvereye-candidate-genes/src/"
PATH_INPUT="/data/Users/Andrea/silvereye/wgs_beagle/raw/"
PATH_OUT="~/sjoh4959/projects/0.0_silvereye-candidate-genes/analyses/NGSadmix/"

for i in {1..10};
        do nohup "${PATH_NGSADMIX}BA3-SNPS-Ubuntu64" --file "${PATH_INPUT}mel.inp" -t -s $RANDOM -m0.1 -f0.0001 -a0.6 -o "${PATH_OUT}out_mel_$i.txt" -i 1000000 -b 100000 -n 100 -v --loci 10000 &
done

for i in {1..10};
        do nohup "${PATH_NGSADMIX}BA3-SNPS-Ubuntu64" --file "${PATH_INPUT}/aus.inp" -t -s $RANDOM -m0.5 -f0.7 -a0.7 -o "${PATH_OUT}out_aus_$i.txt" -i 10000000 -b 1000000 -n 1000 -v --loci 10000 &
done
```



