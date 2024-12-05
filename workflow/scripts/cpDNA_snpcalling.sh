# bowtie_align_run.sh
#!/bin/bash
export read_group=../../../data/readgroup.txt
export genomeindex=../../../data/genome/cpDNA/chrC_index
export datadir=../../../data/raw_data
export outdir=../../../result/mapping/cpDNA

cat $read_group | while read line
do
samplename=$(echo $line | cut -f1 -d " ")
echo $samplename
for fastq in ${datadir}/${samplename}/*R1_*.fastq.gz
do
fastq1="${datadir}/${samplename}/${samplename}_*_R1_*.fastq.gz"
fastq2="${datadir}/${samplename}/${samplename}_*_R2_*.fastq.gz"
echo $fastq1 $fastq2
sbatch bowtie_align.sh $genomeindex $fastq1 $fastq2 ${outdir}/${samplename}.bam
done
done
#--------------------------------------------------------------------------------------------------------------------------------------------------#
#!/bin/bash

export read_group=../../../data/readgroup.txt
export outdir=../../../result/mapping/cpDNA

cat $read_group | while read line
do
samplename=$(echo $line | cut -f1 -d " ")
echo $samplename
    sbatch markduplicate.sh ${outdir}/${samplename}.sorted.bam ${outdir}/${samplename}.sorted.dedup.bam ${outdir}/${samplename}.metrics.txt
done
#--------------------------------------------------------------------------------------------------------------------------------------------------#
#!/bin/bash

export read_group=../../../data/readgroup.txt
export genome=../../../data/genome/cpDNA/chrC.fasta
export bamdir=../../../result/mapping/cpDNA
export snpdir=../../../result/snpcalling/cpDNA

cat $read_group | while read line
do
samplename=$(echo $line | cut -f1 -d " ")
echo $samplename
    sbatch haplotype_caller.sh ${genome} ${bamdir}/${samplename}.sorted.dedup.RG.bam ${snpdir}/${samplename}.gvcf
done
#--------------------------------------------------------------------------------------------------------------------------------------------------#
#!/bin/bash

export read_group=../../../data/readgroup.txt
export genome=../../../data/genome/cpDNA/chrC.fasta
export snpdir=../../../result/snpcalling/cpDNA

cat $read_group | while read line
do
samplename=$(echo $line | cut -f1 -d " ")
echo $samplename
    sbatch genotypeGVCF.sh ${genome} ${snpdir}/${samplename}.gvcf ${snpdir}/${samplename}.vcf.gz
done
#--------------------------------------------------------------------------------------------------------------------------------------------------#
#!/bin/bash -l

#SBATCH -A snic2022-5-643
#SBATCH -p core
#SBATCH -n 1
#SBATCH -t 4:00:00
#SBATCH -J sort
#SBATCH --mail-user=homa.papoli_yazdi@biol.lu.se
#SBATCH --mail-type=FAIL

# Load modules
module load bioinfo-tools vcftools/0.1.16 htslib/1.17 GATK/4.3.0.0

export read_group=../../../data/readgroup.txt
export genome=../../../data/genome/cpDNA/chrC.fasta
export snpdir=../../../result/snpcalling/cpDNA

cat $read_group | while read line
do
samplename=$(echo $line | cut -f1 -d " ")
echo $samplename
    # Keep biallelic SNPs
    vcftools --gzvcf ${snpdir}/${samplename}.vcf.gz \
    --remove-indels --min-alleles 2 --max-alleles 2 \
    --recode --recode-INFO-all \
    --stdout | bgzip -c > ${snpdir}/${samplename}.biallelic.snp.vcf.gz
    # Index VCF
    java -Xmx6g -jar $GATK_HOME/gatk-package-4.3.0.0-local.jar IndexFeatureFile -I ${snpdir}/${samplename}.biallelic.snp.vcf.gz
    # Filter VCF
    java -Xmx6g -jar $GATK_HOME/gatk-package-4.3.0.0-local.jar VariantFiltration \
    -R $genome \
    -V ${snpdir}/${samplename}.biallelic.snp.vcf.gz \
    -O ${snpdir}/${samplename}.biallelic.snp.filterflags.vcf.gz \
    --cluster-window-size 10 --cluster-size 3 \
    --filter-expression "QUAL < 30.0 " --filter-name "VeryLowQual" \
    --filter-expression "QD < 2.0 " --filter-name "LowQD" \
    --filter-expression "FS > 60.0" --filter-name "StrandBiasFishers" \
    --filter-expression "MQ < 40.0" --filter-name "Mapping Quality" \
    --filter-expression "MQRankSum < -12.5" --filter-name "MQRankSum" \
    --filter-expression "ReadPosRankSum < -8.0" --filter-name "ReadPosRankSum" \
    --filter-expression "SOR > 3.0" --filter-name "SOR" 
    # Remove filtered sites
    vcftools --gzvcf ${snpdir}/${samplename}.biallelic.snp.filterflags.vcf.gz --remove-filtered-all --minDP 5 \
    --max-missing 1.0 --recode --recode-INFO-all --stdout | bgzip -c > ${snpdir}/${samplename}.biallelic.snp.filtered.vcf.gz
done
#--------------------------------------------------------------------------------------------------------------------------------------------------#
#!/bin/bash -l

#SBATCH -A snic2022-5-643
#SBATCH -p core
#SBATCH -n 1
#SBATCH -t 4:00:00
#SBATCH -J sort
#SBATCH --mail-user=homa.papoli_yazdi@biol.lu.se
#SBATCH --mail-type=FAIL

# Load modules
module load bioinfo-tools vcftools/0.1.16 htslib/1.17 GATK/4.3.0.0

export read_group=../../../data/readgroup.txt
export genome=../../../data/genome/cpDNA/chrC.fasta
export snpdir=../../../result/snpcalling/cpDNA

cat $read_group | while read line
do
samplename=$(echo $line | cut -f1 -d " ")
echo $samplename
    #java -Xmx6g -jar $GATK_HOME/gatk-package-4.3.0.0-local.jar IndexFeatureFile -I ${snpdir}/${samplename}.biallelic.snp.filtered.vcf.gz
    zgrep -v "#" ${snpdir}/${samplename}.biallelic.snp.filtered.vcf.gz | wc
    cat $genome | vcf-consensus ${snpdir}/${samplename}.biallelic.snp.filtered.vcf.gz | sed "s/chrC/$samplename/g" > ${snpdir}/${samplename}.postfiltersnps.fa
done
#--------------------------------------------------------------------------------------------------------------------------------------------------#
cat ../../../result/snpcalling/cpDNA/*.postfiltersnps.fa > ../../../result/snpcalling/cpDNA/cpDNA.postfiltersnps.fa
#--------------------------------------------------------------------------------------------------------------------------------------------------#
module load bioinfo-tools iqtree/2.2.2.6-omp-mpi
iqtree2 -s ../../../result/snpcalling/cpDNA/cpDNA.postfiltersnps_poly.fa -m GTR+F+I+G4 -bb 1000 -alrt 1000 -nt AUTO