import pandas as pd
import glob
import os
from collections import defaultdict
import sys

# Path to the base directory containing the files
base_dir = "results/phasing"

# Collect all the `.tsv` files in the structure
files = glob.glob(os.path.join(base_dir, "*/*.tsv"))

# Dictionary to track phase sets, contigs, and the samples they occur in
phase_set_contig_samples = defaultdict(lambda: defaultdict(set))

# Process each file
for file in files:
    # Extract the sample name from the file path
    sample_name = os.path.basename(file).replace(".tsv", "")
    
    # Read the file into a DataFrame
    df = pd.read_csv(file, sep="\t")
    
    # Filter rows where variants >= 2
    filtered = df[df["variants"] >= 3]
    
    # Add phase sets to the dictionary, tracking which sample and contig they belong to
    for _, row in filtered.iterrows():
        phase_set = row["phase_set"]
        contig = row["chromosome"]
        
        # Add the sample to the phase_set and contig combination
        phase_set_contig_samples[phase_set][contig].add(sample_name)

# Filter phase sets that occur in at least 3 samples on the same contig
result = {}
for phase_set, contig_samples in phase_set_contig_samples.items():
    for contig, samples in contig_samples.items():
        if len(samples) >= 10:  # At least 3 samples must share this phase_set on the same contig
            if phase_set not in result:
                result[phase_set] = []
            result[phase_set].append((contig, list(samples)))

# Save the result to a file
with open(sys.argv[1], "w") as out_file:
    out_file.write("phase_set\tcontig\tsamples\n")
    for phase_set, contig_samples in result.items():
        for contig, samples in contig_samples:
            out_file.write(f"{phase_set}\t{contig}\t{','.join(samples)}\n")