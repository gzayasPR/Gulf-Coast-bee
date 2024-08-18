from Bio import SeqIO
import matplotlib.pyplot as plt

def gc_content(sequence):
    gc_count = sequence.count('G') + sequence.count('C')
    return (gc_count / len(sequence)) * 100

def calculate_gc_content_windows(fasta_file, window_size):
    gc_contents = []

    for record in SeqIO.parse(fasta_file, "fasta"):
        sequence = record.seq
        for i in range(0, len(sequence) - window_size + 1, window_size):
            window_seq = sequence[i:i+window_size]
            gc_contents.append(gc_content(window_seq))

    return gc_contents

def plot_gc_content_distribution(gc_contents, output_file):
    plt.hist(gc_contents, bins=50, edgecolor='black')
    plt.title('GC Content Distribution')
    plt.xlabel('GC Content (%)')
    plt.ylabel('Frequency')
    plt.savefig(output_file)
    plt.close()

# Replace 'your_genome.fasta' with the path to your FASTA file
fasta_file  = '/90daydata/beenome100/hesperapis_oraria_genomics/Hesperapis_oraria_Pop_genomics/data/final_assembly/Hesperapis_oraria_2.curatedScaff.clean_sort_rename_nuc.fasta'
window_size = 1000  # Define the window size (e.g., 1000 bp)
output_file = 'gc_content_distribution.png'  # Define the output file name

gc_contents = calculate_gc_content_windows(fasta_file, window_size)
plot_gc_content_distribution(gc_contents, output_file)
