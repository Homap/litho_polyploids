import pandas as pd
import subprocess
import os

# Read the phase set file (phase_set_in_10_samples_same_contig.tsv)
phase_set_file = "results/grampa/input_processing/phase_sets_in_10_samples_same_contig.txt"
phase_set_data = pd.read_csv(phase_set_file, sep="\t")

# Function to create VCF file per sample and phase set
def create_vcf_for_sample(phase_set, contig, sample_list):
    # Iterate over each sample in the list
    output_dir = f"results/grampa/input_processing/vcf/{contig}_{phase_set}"
    os.makedirs(output_dir, exist_ok=True)
    print(output_dir, "created")
    for sample in sample_list:
        # Include contig name in the output file name
        sample = sample.replace(".phased", "")
        sample_shortened = "_".join(sample.replace(".phased", "").split("_")[0:2])
        input_vcf = f"results/phasing/{sample}/{sample}.phased.vcf.gz"
        # index_file = f"{input_vcf}.tbi"  
        output_vcf = os.path.join(output_dir, f"{sample_shortened}_{contig}_{phase_set}.vcf.gz")
        
        # Use bcftools to extract the subset of VCF for the sample and phase_set
        # print("Running Bcftools")
        # if not os.path.exists(index_file):
        #     print(f"Index not found for {input_vcf}, creating index...")
        #     bcftools_index_command = f"singularity exec containers/bcftools_1.21.sif bcftools index -t -f {input_vcf}"
        #     subprocess.run(bcftools_index_command, shell=True, check=True)
        # else:
        #     print(f"Index already exists for {input_vcf}, skipping indexing.")

        bcftools_command = f"singularity exec containers/bcftools_1.21.sif bcftools view -s {sample_shortened} -r {contig} -i PS={phase_set} {input_vcf} -o {output_vcf} -Oz"
        
        # Run the command to generate the VCF
        subprocess.run(bcftools_command, shell=True, check=True)
        print(f"Created VCF for sample {sample}, contig {contig}, and phase set {phase_set}: {output_vcf}")

print("Running create_vcf_for_sample function")
# Loop through each phase set entry in the file
for _, row in phase_set_data.iterrows():
    phase_set = row['phase_set']
    contig = row['contig']
    samples = row['samples'].split(',')
    
    # Create VCF files for each sample and phase set
    create_vcf_for_sample(phase_set, contig, samples)