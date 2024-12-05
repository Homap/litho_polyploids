import subprocess
import re
import sys

# File name for the input VCF
vcf_file = sys.argv[1]

# BCFtools command to extract genotypes for each SNP
bcftools_cmd = [
    "singularity", "exec", "containers/bcftools_1.21.sif", "bcftools", 
    "query", "-f", "[%GT\t]\n", vcf_file
]

# Run the BCFtools command and capture output
try:
    process = subprocess.run(
        bcftools_cmd, 
        text=True,  # Ensure output is in string format
        stdout=subprocess.PIPE,  # Capture standard output
        stderr=subprocess.PIPE   # Capture errors
    )
    if process.returncode != 0:
        raise Exception(f"BCFtools error: {process.stderr}")
    
    # Output from BCFtools
    genotype_lines = process.stdout.strip().split("\n")

except Exception as e:
    print(f"Error running BCFtools: {e}")
    exit(1)

# Regular expression to split genotypes by "/" or "|"
allele_splitter = re.compile(r"[\/|]")

# Define groups and columns
diploid_groups = [
    [0, 1, 2, 3],    # Group 1: Columns 1–4
    [4, 5, 6, 7],    # Group 2: Columns 5–8
    [8, 9, 10, 11],  # Group 3: Columns 9–12
    [22, 23, 24, 25] # Group 4: Columns 23–26 (newly added)
]
tetraploid_columns = [12, 13, 14, 15, 16, 17]  # Columns 13–18
hexaploid_columns = [18, 19, 20, 21]  # Columns 19–22

print(" ".join(["BAM", "MIN", "MXN", "GLA", "KNG4x", "MXN4x", "MIN4x", "MBL4x", "BAT4x", "WVR4x", "NHL6x", "SAM6x", "BAT6x", "HOP6x"]))
# Process each line of genotypes from BCFtools output
for line in genotype_lines:
    genotypes = line.strip().split("\t")
    
    # Initialize allele counts for each group and column
    diploid_counts = []
    tetraploid_counts = []
    hexaploid_counts = []
    
    # Handle diploid groups
    for group in diploid_groups:
        group_counts = {"0": 0, "1": 0}
        for col in group:
            genotype = genotypes[col]
            alleles = allele_splitter.split(genotype)
            for allele in alleles:
                group_counts[allele] = group_counts.get(allele, 0) + 1
        diploid_counts.append(f"{group_counts.get('0', 0)},{group_counts.get('1', 0)}")
    
    # Handle tetraploid columns
    for col in tetraploid_columns:
        genotype = genotypes[col]
        alleles = allele_splitter.split(genotype)
        tetraploid_counts.append(f"{alleles.count('0')},{alleles.count('1')}")
    
    # Handle hexaploid columns
    for col in hexaploid_columns:
        genotype = genotypes[col]
        alleles = allele_splitter.split(genotype)
        hexaploid_counts.append(f"{alleles.count('0')},{alleles.count('1')}")
    
    # Combine results for the line (one SNP)
    print(" ".join(diploid_counts + tetraploid_counts + hexaploid_counts))
