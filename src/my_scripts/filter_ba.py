import pandas as pd
import argparse

# Define command line arguments
parser = argparse.ArgumentParser(description='Filter a tab-separated file by a list of sample names.')
parser.add_argument('full_file', type=str, help='Path to the full tab-separated file')
parser.add_argument('sample_names', type=str, help='Path to the file containing the list of sample names')
parser.add_argument('output_file', type=str, help='Path to the output file')

# Parse command line arguments
args = parser.parse_args()

# Load the full tab-separated file into a Pandas dataframe
df = pd.read_csv(args.full_file, sep='\t')

# Load the list of sample names into a Python list
with open(args.sample_names, 'r') as f:
    sample_list = [line.strip() for line in f]

# Filter the dataframe by the list of sample names
filtered_df = df[df.iloc[:, 0].isin(sample_list)]

# Save the filtered dataframe to a new tab-separated file
filtered_df.to_csv(args.output_file, sep='\t', index=False)

