import pysam
import pandas as pd
from tqdm import tqdm

sorted_bam_file = snakemake.input.sorted_bam
gff_file = snakemake.input.gff
output_csv = snakemake.output.csv

# Load GFF file into a DataFrame
gff_df = pd.read_csv(gff_file, sep='\t', header=None, 
                     names=["seqname", "source", "feature", "start", "end", "score", "strand", "frame", "attribute"])

# Open the position sorted BAM file (index must be present)
bam = pysam.AlignmentFile(sorted_bam_file, "rb")

# Initialize a dictionary to count reads in each bin
bin_counts = {}

for seq in ['2L', '2R', '3L', '3R', '4', 'X', 'Y']:
    print(seq, gff_df[gff_df.seqname==seq].shape[0])
    # Iterate through the BAM file only once per chromosome
    current_bin_idx = 0
    current_bin = gff_df.iloc[current_bin_idx]
    bin_start = current_bin['start']
    bin_end = current_bin['end']
    bin_id = f"{seq}:{bin_start}-{bin_end}"
    bin_counts[bin_id] = 0
    
    for read in tqdm(bam.fetch(seq), total=bam.mapped):
        read_start = read.reference_start
        read_end = read.reference_end

        # Move to the next bin if the read is past the current bin
        while read_start > bin_end and current_bin_idx < len(gff_df) - 1:
            # Go to the next bin
            current_bin_idx += 1
            current_bin = gff_df.iloc[current_bin_idx]
            bin_start = current_bin['start']
            bin_end = current_bin['end']
            bin_id = f"{seq}:{bin_start}-{bin_end}"
            
            # Initialize the bin count if it is not in the dictionary
            if bin_id not in bin_counts:
                bin_counts[bin_id] = 0

        # Count the read if it is within the current bin
        if read_start <= bin_end and read_end >= bin_start:
            bin_counts[bin_id] += 1

# Close the BAM file
bam.close()

pd.Series(bin_counts).to_csv(output_csv, header=False)