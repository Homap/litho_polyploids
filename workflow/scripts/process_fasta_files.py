import os
from collections import defaultdict

# List of populations to check
pops = ["BAM", "MIN", "MXN", "GLA", "KNG4x", "MXN4x", "MIN4x", 
        "MBL4x", "BAT4x", "WVR4x", "NHL6x", "SAM6x", "BAT6x", "HOP6x"]

def parse_fasta(file):
    """Parses a FASTA file and returns a dictionary of sequences by header."""
    sequences = defaultdict(str)
    current_header = None
    with open(file, 'r') as f:
        for line in f:
            line = line.strip()
            if line.startswith(">"):
                current_header = line[1:]
            elif current_header:
                sequences[current_header] += line
    print(f"Parsed {len(sequences)} sequences from {file}")
    return sequences

def find_unique_sequences(sequences, population_tags):
    """Finds unique sequences for each population and assigns them unique numbers."""
    unique_sequences = defaultdict(list)  # To store sequences for each population
    seen_sequences = defaultdict(set)  # To track seen sequences for each population
    
    # Iterate over all sequences
    for header, seq in sequences.items():
        # Identify which population the sequence belongs to
        matched_population = None
        for pop in population_tags:
            if pop in header:
                matched_population = pop
                break

        if matched_population:  # If a population match was found
            # Only add the sequence if it hasn't been seen before for this population
            if seq not in seen_sequences[matched_population]:
                unique_sequences[matched_population].append(seq)
                seen_sequences[matched_population].add(seq)
    
    # Assign unique numbers to each sequence per population
    numbered_sequences = defaultdict(dict)
    for pop, seqs in unique_sequences.items():
        for idx, seq in enumerate(seqs, start=1):
            numbered_sequences[pop][f"{idx}_{pop}"] = seq

    return numbered_sequences

def write_fasta(numbered_sequences, output_file):
    """Writes unique sequences for all populations to a single FASTA file."""
    total_written = 0
    with open(output_file, 'w') as f:
        for pop, seqs in numbered_sequences.items():
            for header, seq in seqs.items():
                f.write(f">{header}\n{seq}\n")
                total_written += 1
    return total_written

def process_fasta_files(input_dir, output_dir):
    """Processes each FASTA file and creates one output file per input file."""
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)

    # Iterate through all FASTA files in the input directory
    for file in os.listdir(input_dir):
        if file.endswith(".fasta"):
            print(f"Processing file: {file}")
            input_file = os.path.join(input_dir, file)
            sequences = parse_fasta(input_file)

            # Find unique sequences for all populations
            numbered_sequences = find_unique_sequences(sequences, pops)

            # Write to a single output file
            output_file = os.path.join(output_dir, f"{file.rsplit('.', 1)[0]}_unique.fasta")
            total_written = write_fasta(numbered_sequences, output_file)

            if total_written == 0:
                print(f"No unique sequences found in {file}. Output file is empty.")
            else:
                print(f"Written {total_written} unique sequences to {output_file}")

# Input and output settings
input_directory = "in_fasta"  # Directory containing your input FASTA files
output_directory = "out_fasta"  # Directory to store the output FASTA files

# Process the FASTA files
process_fasta_files(input_directory, output_directory)

print(f"Processing completed. Unique sequences have been written to {output_directory}")
