import csv
import re
import gzip
from collections import Counter
from tqdm import tqdm

#iterate line by line, extract the
#PATH_SAMPLES='/home/zoo/sjoh4959/sjoh4959/projects/0.0_silvereye-candidate-genes/data/samples2remove.csv'
#PATH_BEAGLE='/data/Users/Andrea/ZLat.v1_WholeGenome.beagle.gz'
#PATH_OUTPUT='/data/Users/Andrea/ZLat.v1_WholeGenome_subset.beagle'

def subset_beagle(PATH_BEAGLE, PATH_SAMPLES, PATH_OUTPUT):
    """This function removes individuals from a BEAGLE file

    Args:
        PATH_BEAGLE (PosixPath): Path to the BEAGLE file in compressed format (.gz)
        PATH_SAMPLES (PosixPath): Path to samples to remove
        PATH_OUTPUT (PosixPath): Path to output file. This path should include the name of the file. If it does not exist it will create it. If it does exist it will overwrite it.
    """

#Open list of samples to keep
    with open(PATH_SAMPLES, newline='') as f:
        reader = f.read().split('\n')
        samples2keep = list(reader)

#Open the BEAGLE file by line, keep the three first columns as they are 
#For the rest, search and keep those that coincide with the ones to keep
#Write line by line the three first cols + those samples to keep
    with open(PATH_OUTPUT, 'w') as out:
        with gzip.GzipFile(PATH_BEAGLE, 'r') as input:
            for line in tqdm(input, position=0, leave=True): 
                linestr = re.split(r'\t+', line.decode("utf-8"))
                first_cols = linestr[:3]
                if line == 0:
                    common_idx = [i for i, item in enumerate(linestr) if item in samples2keep]
                a = []
                for index in range(len(common_idx)):
                    a.append(linestr[common_idx[index]])
                out.write('\t'.join(map(str,first_cols))+'\t'+'\t'.join(map(str,a))+'\n')


