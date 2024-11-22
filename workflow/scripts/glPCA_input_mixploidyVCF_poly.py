import gzip
import sys
#-------------------------------------------------------------------------------------------------------------------------------------------------------------------
# Sample names
#-------------------------------------------------------------------------------------------------------------------------------------------------------------------
diploid_sample_names = ["P18758_101_S1_L001", "P18758_102_S2_L001", "P18758_103_S3_L001", "P18758_104_S4_L001", 
"P18758_131_S31_L001", "P18758_132_S32_L001", "P18758_133_S33_L002", "P18758_134_S34_L002", 
"P18758_161_S61_L002", "P18758_162_S62_L002", "P18758_163_S63_L002", "P18758_164_S64_L002",
"P18758_201_S89_L004", "P18758_202_S90_L004", "P18758_203_S91_L004", "P18758_204_S92_L004"]
tetraploid_sample_names = ["P18758_169_S79_L004", "P18758_170_S80_L004", "P18758_171_S81_L004", "P18758_172_S82_L004", "P18758_173_S83_L004", "P18758_174_S84_L004"]
hexaploid_sample_names = ["P18758_175_S85_L004", "P18758_176_S86_L004", "P18758_177_S87_L004", "P18758_178_S88_L004"]

# Function to extract genotype from genotype field
def extract_genotype(genotype_field):
    return genotype_field.split(':')[0]  # Get only the genotype part

# Initialize dictionary to store sample genotypes
sample_genotypes = {}
geno_chr = []
geno_position = []
geno_alleles = []

# Open the gzipped VCF file
# with gzip.open("filtered_no_missing.prunned.vcf.gz", 'rt') as vcf_file:
with gzip.open(sys.argv[1], 'rt') as vcf_file:
    for line in vcf_file:
        # Skip header lines
        if line.startswith('#'):
            # Parse the sample IDs from the header line
            if line.startswith('#CHROM'):
                headers = line.strip().split('\t')
                sample_ids = headers[9:]  # Samples start from the 10th column
                # Initialize the dictionary with sample IDs
                sample_genotypes = {sample: [] for sample in sample_ids}
            continue
        
        # Parse data lines
        fields = line.strip().split('\t')
        geno_chr.append(fields[0])
        geno_position.append(fields[1])
        alleles = fields[3].lower()+"/"+fields[4].lower()
        geno_alleles.append(alleles)
        genotypes = fields[9:]  # Genotype data starts from the 10th column
        
        # Store genotypes in the dictionary
        for sample, genotype_field in zip(sample_ids, genotypes):
            genotype = extract_genotype(genotype_field)
            sample_genotypes[sample].append(genotype)

# P18758_169	KNG11544A_4x	P18758_169_S79_L004	4
# P18758_170	MXN12714A_4x	P18758_170_S80_L004	4
# P18758_171	MIN11453D_4x	P18758_171_S81_L004	4
# P18758_172	MBL8857A_4x	P18758_172_S82_L004	4
# P18758_173	BAT13308A_4x	P18758_173_S83_L004	4
# P18758_174	WVR12731B_4x	P18758_174_S84_L004	4
# P18758_175	NHL12853B_6x	P18758_175_S85_L004	6
# P18758_176	SAM13220A_6x	P18758_176_S86_L004	6
# P18758_177	BAT13299B_6x	P18758_177_S87_L004	6
# P18758_178	HOP966A_6x	P18758_178_S88_L004	6

print(">>>> begin comments - do not remove this line <<<<")
print(">>>> end comments - do not remove this line <<<<")
print(">>position")
print(" ".join(geno_position[::50]))
print(">>allele")
print(" ".join(geno_alleles[::50]))
print(">>chromosome")
print(" ".join(geno_chr[::50]))
print(">>population")
print("NHL6x", "SAM6x", "BAT6x", "HOP6x", "KNG4x", "MXN4x", "MIN4x", "MBL4x", "BAT4x", "WVR4x",
"BAM", "BAM", "BAM", "BAM",
"MIN", "MIN", "MIN", "MIN", 
"MXN", "MXN", "MXN", "MXN", 
"GLA", "GLA", "GLA", "GLA")
print(">>ploidy")
print(" ".join(["6","6","6","6", "4", "4", "4", "4", "4", "4", "2", "2","2","2","2","2","2","2","2","2","2","2","2","2","2","2"]))

for sample in hexaploid_sample_names:
    ind_genotypes_hexa = []
    genotype_list = sample_genotypes[sample]  # This should be a list of genotypes
    # Iterate through each genotype for the sample
    for genotype_hexaploid in genotype_list:
        if genotype_hexaploid == "0/0/0/0/0/0":
            ind_genotypes_hexa.append("0")
        elif genotype_hexaploid == "0/0/0/0/0/1":
            ind_genotypes_hexa.append("1")
        elif genotype_hexaploid == "0/0/0/0/1/1":
            ind_genotypes_hexa.append("2")
        elif genotype_hexaploid == "0/0/0/1/1/1":
            ind_genotypes_hexa.append("3")
        elif genotype_hexaploid == "0/0/1/1/1/1":
            ind_genotypes_hexa.append("4")
        elif genotype_hexaploid == "0/1/1/1/1/1":
            ind_genotypes_hexa.append("5")
        elif genotype_hexaploid == "1/1/1/1/1/1":
            ind_genotypes_hexa.append("6")
        elif genotype_hexaploid == "./././././.":
            ind_genotypes_hexa.append("-")
        elif genotype_hexaploid == "./.":
            ind_genotypes_hexa.append("-")
        else:
            print(genotype_hexaploid, "Error")
    # Corrected string concatenation
    print(">" + "_".join(sample.split("_")[0:2]))
    # Print the reduced genotype list (every 50th element)
    print("".join(ind_genotypes_hexa))

for sample in tetraploid_sample_names:
    ind_genotypes_tetraploid = []
    genotype_list = sample_genotypes[sample]  # This should be a list of genotypes
    # Iterate through each genotype for the sample
    for genotype_tetraploid in genotype_list:
        if genotype_tetraploid == "0/0/0/0":
            ind_genotypes_tetraploid.append("0")
        elif genotype_tetraploid == "0/0/0/1":
            ind_genotypes_tetraploid.append("1")
        elif genotype_tetraploid == "0/0/1/1":
            ind_genotypes_tetraploid.append("2")
        elif genotype_tetraploid == "0/1/1/1":
            ind_genotypes_tetraploid.append("3")
        elif genotype_tetraploid == "1/1/1/1":
            ind_genotypes_tetraploid.append("4")
        elif genotype_tetraploid == "./.":
            ind_genotypes_tetraploid.append("-")
        elif genotype_tetraploid == "./././.":
            ind_genotypes_tetraploid.append("-")
        else:
            print(genotype_tetraploid, "Error")
    # Corrected string concatenation
    print(">" + "_".join(sample.split("_")[0:2]))
    # Print the reduced genotype list (every 50th element)
    print("".join(ind_genotypes_tetraploid))

for sample in diploid_sample_names:
    ind_genotypes_diploid = []
    genotype_list = sample_genotypes[sample]  # This should be a list of genotypes
    
    # Iterate through each genotype for the sample
    for genotype_diploid in genotype_list:
        if genotype_diploid == "0/0" or genotype_diploid == "0|0":
            ind_genotypes_diploid.append("0")
        elif genotype_diploid == "0/1" or genotype_diploid == "0|1":
            ind_genotypes_diploid.append("1")
        elif genotype_diploid == "1/1" or genotype_diploid == "1|1":
            ind_genotypes_diploid.append("2")
        elif genotype_diploid == "./." or genotype_diploid == ".|.":
            ind_genotypes_diploid.append("-")
        else:
            print(genotype_diploid, "Error")
    # Corrected string concatenation
    print(">" + "_".join(sample.split("_")[0:2]))
    # Print the reduced genotype list (every 50th element)
    print("".join(ind_genotypes_diploid))