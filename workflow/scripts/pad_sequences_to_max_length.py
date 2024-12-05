import sys
from Bio import SeqIO

def pad_sequences_to_max_length(input_fasta, output_fasta):
    # Read all sequences from the input FASTA file
    sequences = list(SeqIO.parse(input_fasta, "fasta"))
    
    # Find the length of the longest sequence
    max_length = max(len(record.seq) for record in sequences)
    
    # Create a list of padded sequences
    for record in sequences:
        # If the sequence is shorter than the longest sequence, pad it with '-'
        if len(record.seq) < max_length:
            record.seq = record.seq + "-" * (max_length - len(record.seq))
    
    # Write the padded sequences to the output FASTA file
    SeqIO.write(sequences, output_fasta, "fasta")

# Example usage
input_fasta = sys.argv[1] #"input_sequences.fasta"  # Replace with your input FASTA file path
output_fasta = sys.argv[2] #"output_sequences.fasta"  # Replace with your desired output FASTA file path

pad_sequences_to_max_length(input_fasta, output_fasta)