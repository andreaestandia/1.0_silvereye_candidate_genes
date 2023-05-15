import pandas as pd
import numpy as np
import argparse

# Parse command line arguments
parser = argparse.ArgumentParser()
parser.add_argument("input_file", help="path to input file")
parser.add_argument("num_markers", type=int, help="number of markers to select")
parser.add_argument("output_file", help="path to output file")
args = parser.parse_args()

# Load the input file into a pandas dataframe
df = pd.read_csv(args.input_file, sep="\t", header=None)

# Select the unique markers
unique_markers = df[2].unique()

# Randomly select the specified number of markers
if len(unique_markers) < args.num_markers:
    raise ValueError("Number of unique markers in the input file is less than the specified number of markers to select.")
selected_markers = np.random.choice(unique_markers, size=args.num_markers, replace=False)

# Filter the dataframe to include only the selected markers with all individuals
selected_df = df[df[2].isin(selected_markers)]
#grouped_df = selected_df.groupby(2).filter(lambda x: len(x) == len(df[4].unique()))

# Write the selected dataframe to the specified output file
selected_df.to_csv(args.output_file, sep="\t", header=False, index=False)
print(f"{len(selected_df[2].unique())} SNPs with all individuals were selected and written to {args.output_file}.")

