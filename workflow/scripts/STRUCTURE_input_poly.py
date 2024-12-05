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
# ----------------------------------------------------------------------------------------------------------------------------- #
# Hexaploids
# ----------------------------------------------------------------------------------------------------------------------------- #
for sample in hexaploid_sample_names:
    ind_genotypes_hexaploid_allele1 = []
    ind_genotypes_hexaploid_allele2 = []
    ind_genotypes_hexaploid_allele3 = []
    ind_genotypes_hexaploid_allele4 = []
    ind_genotypes_hexaploid_allele5 = []
    ind_genotypes_hexaploid_allele6 = []
    genotype_list = sample_genotypes[sample]  # This should be a list of genotypes
    # Iterate through each genotype for the sample
    for genotype_hexaploid in genotype_list:
        if genotype_hexaploid == "0/0/0/0/0/0":
            ind_genotypes_hexaploid_allele1.append("0")
            ind_genotypes_hexaploid_allele2.append("0")
            ind_genotypes_hexaploid_allele3.append("0")
            ind_genotypes_hexaploid_allele4.append("0")
            ind_genotypes_hexaploid_allele5.append("0")
            ind_genotypes_hexaploid_allele6.append("0")
        elif genotype_hexaploid == "0/0/0/0/0/1":
            ind_genotypes_hexaploid_allele1.append("0")
            ind_genotypes_hexaploid_allele2.append("0")
            ind_genotypes_hexaploid_allele3.append("0")
            ind_genotypes_hexaploid_allele4.append("0")
            ind_genotypes_hexaploid_allele5.append("0")
            ind_genotypes_hexaploid_allele6.append("1")
        elif genotype_hexaploid == "0/0/0/0/1/1":
            ind_genotypes_hexaploid_allele1.append("0")
            ind_genotypes_hexaploid_allele2.append("0")
            ind_genotypes_hexaploid_allele3.append("0")
            ind_genotypes_hexaploid_allele4.append("0")
            ind_genotypes_hexaploid_allele5.append("1")
            ind_genotypes_hexaploid_allele6.append("1")
        elif genotype_hexaploid == "0/0/0/1/1/1":
            ind_genotypes_hexaploid_allele1.append("0")
            ind_genotypes_hexaploid_allele2.append("0")
            ind_genotypes_hexaploid_allele3.append("0")
            ind_genotypes_hexaploid_allele4.append("1")
            ind_genotypes_hexaploid_allele5.append("1")
            ind_genotypes_hexaploid_allele6.append("1")
        elif genotype_hexaploid == "0/0/1/1/1/1":
            ind_genotypes_hexaploid_allele1.append("0")
            ind_genotypes_hexaploid_allele2.append("0")
            ind_genotypes_hexaploid_allele3.append("1")
            ind_genotypes_hexaploid_allele4.append("1")
            ind_genotypes_hexaploid_allele5.append("1")
            ind_genotypes_hexaploid_allele6.append("1")
        elif genotype_hexaploid == "0/1/1/1/1/1":
            ind_genotypes_hexaploid_allele1.append("0")
            ind_genotypes_hexaploid_allele2.append("1")
            ind_genotypes_hexaploid_allele3.append("1")
            ind_genotypes_hexaploid_allele4.append("1")
            ind_genotypes_hexaploid_allele5.append("1")
            ind_genotypes_hexaploid_allele6.append("1")
        elif genotype_hexaploid == "1/1/1/1/1/1":
            ind_genotypes_hexaploid_allele1.append("1")
            ind_genotypes_hexaploid_allele2.append("1")
            ind_genotypes_hexaploid_allele3.append("1")
            ind_genotypes_hexaploid_allele4.append("1")
            ind_genotypes_hexaploid_allele5.append("1")
            ind_genotypes_hexaploid_allele6.append("1")
        elif genotype_hexaploid == "./././././.":
            ind_genotypes_hexaploid_allele1.append("-9")
            ind_genotypes_hexaploid_allele2.append("-9")
            ind_genotypes_hexaploid_allele3.append("-9")
            ind_genotypes_hexaploid_allele4.append("-9")
            ind_genotypes_hexaploid_allele5.append("-9")
            ind_genotypes_hexaploid_allele6.append("-9")
        else:
            print(genotype_hexaploid, "Error")
    # Corrected string concatenation
    print(sample+" "+" ".join(ind_genotypes_hexaploid_allele1[::15])) # First allele of all loci
    print(sample+" "+" ".join(ind_genotypes_hexaploid_allele2[::15])) # Second allele of all loci
    print(sample+" "+" ".join(ind_genotypes_hexaploid_allele3[::15])) # Third allele of all loci # Fake hexaploid by adding missing allele
    print(sample+" "+" ".join(ind_genotypes_hexaploid_allele4[::15])) # Fourth allele of all loci # Fake hexaploid by adding missing allele
    print(sample+" "+" ".join(ind_genotypes_hexaploid_allele5[::15])) # Fifth allele of all loci # Fake hexaploid by adding missing allele
    print(sample+" "+" ".join(ind_genotypes_hexaploid_allele6[::15])) # Sixth allele of all loci # Fake hexaploid by adding missing allele
# ----------------------------------------------------------------------------------------------------------------------------- #
# Tetraploids
# ----------------------------------------------------------------------------------------------------------------------------- #
for sample in tetraploid_sample_names:
    ind_genotypes_tetraploid_allele1 = []
    ind_genotypes_tetraploid_allele2 = []
    ind_genotypes_tetraploid_allele3 = []
    ind_genotypes_tetraploid_allele4 = []
    ind_genotypes_tetraploid_allele5 = []
    ind_genotypes_tetraploid_allele6 = []
    genotype_list = sample_genotypes[sample]  # This should be a list of genotypes
    # Iterate through each genotype for the sample
    for genotype_tetraploid in genotype_list:
        if genotype_tetraploid == "0/0/0/0":
            ind_genotypes_tetraploid_allele1.append("0")
            ind_genotypes_tetraploid_allele2.append("0")
            ind_genotypes_tetraploid_allele3.append("0")
            ind_genotypes_tetraploid_allele4.append("0")
            ind_genotypes_tetraploid_allele5.append("-9")
            ind_genotypes_tetraploid_allele6.append("-9")
        elif genotype_tetraploid == "0/0/0/1":
            ind_genotypes_tetraploid_allele1.append("0")
            ind_genotypes_tetraploid_allele2.append("0")
            ind_genotypes_tetraploid_allele3.append("0")
            ind_genotypes_tetraploid_allele4.append("1")
            ind_genotypes_tetraploid_allele5.append("-9")
            ind_genotypes_tetraploid_allele6.append("-9")
        elif genotype_tetraploid == "0/0/1/1":
            ind_genotypes_tetraploid_allele1.append("0")
            ind_genotypes_tetraploid_allele2.append("0")
            ind_genotypes_tetraploid_allele3.append("1")
            ind_genotypes_tetraploid_allele4.append("1")
            ind_genotypes_tetraploid_allele5.append("-9")
            ind_genotypes_tetraploid_allele6.append("-9")
        elif genotype_tetraploid == "0/1/1/1":
            ind_genotypes_tetraploid_allele1.append("0")
            ind_genotypes_tetraploid_allele2.append("1")
            ind_genotypes_tetraploid_allele3.append("1")
            ind_genotypes_tetraploid_allele4.append("1")
            ind_genotypes_tetraploid_allele5.append("-9")
            ind_genotypes_tetraploid_allele6.append("-9")
        elif genotype_tetraploid == "1/1/1/1":
            ind_genotypes_tetraploid_allele1.append("1")
            ind_genotypes_tetraploid_allele2.append("1")
            ind_genotypes_tetraploid_allele3.append("1")
            ind_genotypes_tetraploid_allele4.append("1")
            ind_genotypes_tetraploid_allele5.append("-9")
            ind_genotypes_tetraploid_allele6.append("-9")
        elif genotype_tetraploid == "./././.":
            ind_genotypes_tetraploid_allele1.append("-9")
            ind_genotypes_tetraploid_allele2.append("-9")
            ind_genotypes_tetraploid_allele3.append("-9")
            ind_genotypes_tetraploid_allele4.append("-9")
            ind_genotypes_tetraploid_allele5.append("-9")
            ind_genotypes_tetraploid_allele6.append("-9")
        else:
            print(genotype_tetraploid, "Error")
    # Corrected string concatenation
    # Corrected string concatenation
    print(sample+" "+" ".join(ind_genotypes_tetraploid_allele1[::15])) # First allele of all loci
    print(sample+" "+" ".join(ind_genotypes_tetraploid_allele2[::15])) # Second allele of all loci
    print(sample+" "+" ".join(ind_genotypes_tetraploid_allele3[::15])) # Third allele of all loci # Fake hexaploid by adding missing allele
    print(sample+" "+" ".join(ind_genotypes_tetraploid_allele4[::15])) # Fourth allele of all loci # Fake hexaploid by adding missing allele
    print(sample+" "+" ".join(ind_genotypes_tetraploid_allele5[::15])) # Fifth allele of all loci # Fake hexaploid by adding missing allele
    print(sample+" "+" ".join(ind_genotypes_tetraploid_allele6[::15])) # Sixth allele of all loci # Fake hexaploid by adding missing allele
# ----------------------------------------------------------------------------------------------------------------------------- #
# Diploids
# ----------------------------------------------------------------------------------------------------------------------------- #
for sample in diploid_sample_names:
    ind_genotypes_diploid_allele1 = []
    ind_genotypes_diploid_allele2 = []
    ind_genotypes_diploid_allele3 = []
    ind_genotypes_diploid_allele4 = []
    ind_genotypes_diploid_allele5 = []
    ind_genotypes_diploid_allele6 = []
    genotype_list = sample_genotypes[sample]  # This should be a list of genotypes
    
    # Iterate through each genotype for the sample
    for genotype_diploid in genotype_list:
        if genotype_diploid == "0/0" or genotype_diploid == "0|0":
            ind_genotypes_diploid_allele1.append("0")
            ind_genotypes_diploid_allele2.append("0")
            ind_genotypes_diploid_allele3.append("-9")
            ind_genotypes_diploid_allele4.append("-9")
            ind_genotypes_diploid_allele5.append("-9")
            ind_genotypes_diploid_allele6.append("-9")
        elif genotype_diploid == "0/1" or genotype_diploid == "0|1":
            ind_genotypes_diploid_allele1.append("0")
            ind_genotypes_diploid_allele2.append("1")
            ind_genotypes_diploid_allele3.append("-9")
            ind_genotypes_diploid_allele4.append("-9")
            ind_genotypes_diploid_allele5.append("-9")
            ind_genotypes_diploid_allele6.append("-9")
        elif genotype_diploid == "1/1" or genotype_diploid == "1|1":
            ind_genotypes_diploid_allele1.append("1")
            ind_genotypes_diploid_allele2.append("1")
            ind_genotypes_diploid_allele3.append("-9")
            ind_genotypes_diploid_allele4.append("-9")
            ind_genotypes_diploid_allele5.append("-9")
            ind_genotypes_diploid_allele6.append("-9")
        elif genotype_diploid == "./." or genotype_diploid == ".|.":
            ind_genotypes_diploid_allele1.append("-9")
            ind_genotypes_diploid_allele2.append("-9")
            ind_genotypes_diploid_allele3.append("-9")
            ind_genotypes_diploid_allele4.append("-9")
            ind_genotypes_diploid_allele5.append("-9")
            ind_genotypes_diploid_allele6.append("-9")
        else:
            print(genotype_diploid, "Error")
    # Corrected string concatenation
    print(sample+" "+" ".join(ind_genotypes_diploid_allele1[::15])) # First allele of all loci
    print(sample+" "+" ".join(ind_genotypes_diploid_allele2[::15])) # Second allele of all loci
    print(sample+" "+" ".join(ind_genotypes_diploid_allele3[::15])) # Third allele of all loci # Fake hexaploid by adding missing allele
    print(sample+" "+" ".join(ind_genotypes_diploid_allele4[::15])) # Fourth allele of all loci # Fake hexaploid by adding missing allele
    print(sample+" "+" ".join(ind_genotypes_diploid_allele5[::15])) # Fifth allele of all loci # Fake hexaploid by adding missing allele
    print(sample+" "+" ".join(ind_genotypes_diploid_allele6[::15])) # Sixth allele of all loci # Fake hexaploid by adding missing allele
