import sys
import os
from Bio import SeqIO
from itertools import islice

def chunks(data, SIZE=10000):
    it = iter(data)
    for i in range(0, len(data), SIZE):
        yield {k:data[k] for k in islice(it, SIZE)}

with open(sys.argv[1], "r") as fasta_in:
    record_dict = SeqIO.to_dict(SeqIO.parse(fasta_in, "fasta"))

num_intervals=int(sys.argv[2])
output_path=str(sys.argv[3])

counter=0
for item in chunks(record_dict, num_intervals):
    counter+=1
    print("Interval", counter)
    outfile_name = "interval_" + str(counter) + ".bed"
    outfile = open(os.path.join(output_path, outfile_name), "w")
    for key in item.keys():
        outfile.write(key+"\t"+"0"+"\t"+str(len(item[key]))+"\n")
