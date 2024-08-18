import sys
from Bio import SeqIO

def calculate_gc(seq):
    return float(seq.count("G") + seq.count("C")) / len(seq) * 100

def filter_reads_by_gc(input_fasta, output_fasta, gc_min, gc_max):
    with open(output_fasta, "w") as output_handle:
        for record in SeqIO.parse(input_fasta, "fastq"):
            gc_content = calculate_gc(record.seq)
            if gc_min <= gc_content <= gc_max:
                SeqIO.write(record, output_handle, "fastq")

if __name__ == "__main__":
    input_fastq = sys.argv[1]
    output_fastq_gc_32 = sys.argv[2]
    output_fastq_gc_41 = sys.argv[3]
    
    # Filter reads with GC content around 32%
    filter_reads_by_gc(input_fastq, output_fastq_gc_32, 31, 33)

    # Filter reads with GC content around 41%
    filter_reads_by_gc(input_fastq, output_fastq_gc_41, 40, 42)
