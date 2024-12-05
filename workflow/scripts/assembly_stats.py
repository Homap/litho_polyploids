#!/usr/bin/python
import sys
import pandas as pd
from Bio import SeqIO
from Bio.SeqUtils import gc_fraction
from N50 import calculate_N50
from chunks import chunks

species = sys.argv[1]

with open(sys.argv[2], "r") as fasta_in, open(sys.argv[3], "w") as fasta_edited, open(sys.argv[4], "w") as stats_file:
        print("Reading Fasta Genome into Dictionary ...")
        record_dict = SeqIO.to_dict(SeqIO.parse(fasta_in, "fasta"))
        genome_seq_list=[]
        genome_seq_len=[]
        print("Writing Edited Fasta Genome into a new file ...")
        for key in record_dict.items():
                genome_seq_list.append(str(key[1].seq))
                genome_seq_len.append(len(str(key[1].seq)))
                fasta_edited.write(">contig_"+key[0].replace("|", "_").replace(" ", "_").replace("_arrow", "")+"\n")
                for line in chunks(str(key[1].seq), 80):
                        fasta_edited.write(line.upper() + '\n')
        print("Calculate Fasta Genome Statistics ...")
        genome_seq_string="".join(genome_seq_list)
        A=genome_seq_string.upper().count("A")
        T=genome_seq_string.upper().count("T")
        C=genome_seq_string.upper().count("C")
        G=genome_seq_string.upper().count("G")
        N=genome_seq_string.upper().count("N")
        num_scaffolds=len(genome_seq_len)
        GC_percent=gc_fraction(genome_seq_string)*100
        genome_seq_len_wN = len(genome_seq_string)
        genome_seq_len_woN = A+T+C+G
        N_percent=N/genome_seq_len_wN
        N50=calculate_N50(genome_seq_len)
        print("Writing Fastas Genome Statistics ...")
        stats_file.write(species+"\t"+str(num_scaffolds)+"\t"+str(genome_seq_len_wN)+"\t"+str(genome_seq_len_woN)+"\t"+str(round(GC_percent, 3))+"\t"+str(round(N_percent, 3))+"\t"+str(N50)+"\n")
        print("Done! And Well Done You!")
