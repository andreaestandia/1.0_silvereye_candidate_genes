import sys

# Define a function to translate ANGSD alleles to VCF alleles
def translate_alleles(allele1, allele2):
    if allele1 == '0':
        ref = 'A'
    elif allele1 == '1':
        ref = 'C'
    elif allele1 == '2':
        ref = 'G'
    elif allele1 == '3':
        ref = 'T'
    else:
        raise ValueError(f"Unrecognized allele: {allele1}")
    
    if allele2 == '0':
        alt = 'A'
    elif allele2 == '1':
        alt = 'C'
    elif allele2 == '2':
        alt = 'G'
    elif allele2 == '3':
        alt = 'T'
    else:
        raise ValueError(f"Unrecognized allele: {allele2}")
    
    return ref, alt

# Define the input and output files
input_file = sys.argv[1]
output_file = sys.argv[2]

# Define the list of sample names
sample_names = []

# Open the input and output files
with open(input_file, 'r') as infile, open(output_file, 'w') as outfile:
    
    # Process the header lines
    for line in infile:
        if line.startswith('marker'):
            # Extract the sample names from the header line
            sample_names = line.strip().split()[3::3]
            break
    
    # Write the VCF header
    outfile.write('##fileformat=VCFv4.2\n')
    outfile.write('##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n')
    outfile.write('#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT')
    
    # Write the sample names to the header line
    for name in sample_names:
        outfile.write(f"\t{name}")
    
    # Process the data lines
    for line in infile:
        # Split the line into columns
        cols = line.strip().split()
        
        # Translate the alleles
        try:
            ref, alt = translate_alleles(cols[1], cols[2])
        except ValueError as e:
            print(str(e))
            continue
        
        # Write the VCF record
        chrom, pos = cols[0].split('_')
        outfile.write(f"\n{chrom}\t{pos}\t.\t{ref}\t{alt}\t.\t.\t.\tGT")
        
        # Parse the genotype probabilities for each sample
        for i in range(3, len(cols), 3):
            probs = [float(cols[i]), float(cols[i+1]), float(cols[i+2])]
            gt = '0/0' if probs[0] > probs[1] and probs[0] > probs[2] else \
                 '0/1' if probs[1] > probs[0] and probs[1] > probs[2] else \
                 '1/1'
            outfile.write(f"\t{gt}")
