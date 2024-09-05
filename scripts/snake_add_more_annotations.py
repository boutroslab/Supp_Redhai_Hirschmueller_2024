import pandas as pd
from gtfparse import read_gtf
from tqdm import tqdm

peakfile = snakemake.input.peakfile
annotated_file = snakemake.input.annotated_file
gtf = snakemake.input.gtf
output = snakemake.output[0]

# Basically peaks in signal vs peaks in input
df = pd.read_csv(peakfile, sep='\t', header=None)
df.columns = ["chrom", "chromStart", "chromEnd", "name", "score", "strand", "signalValue", "pValue", "qValue", "peak"]
df.set_index("name", inplace=True)

# add annotation from ChIPseeker
df_anno = pd.read_csv(annotated_file, index_col=0)
df_anno.set_index("V4", inplace=True)
df['annotation'] = df_anno['annotation']
df['geneStart'] = df_anno['geneStart']
df['geneEnd'] = df_anno['geneEnd']
df['geneLength'] = df_anno['geneLength']
df['geneId'] = df_anno['geneId']
df['transcriptId'] = df_anno['transcriptId']
df['distanceToTSS'] = df_anno['distanceToTSS']

# annotate genes with gtf
gtf = read_gtf(gtf)
gtf = pd.DataFrame(gtf, columns=gtf.columns)
gene_gtf = gtf[gtf.feature == "gene"]
for chrom in tqdm(gene_gtf.seqname.unique()):
    chrom_gene_gtf = gene_gtf[gene_gtf.seqname == chrom]
    chrom_df = df[df.chrom == chrom]
    for row, index in chrom_df.iterrows():
        chromStart = index["chromStart"]
        chromEnd = index["chromEnd"]
        gene = chrom_gene_gtf[(chrom_gene_gtf["start"] <= chromStart) & (chrom_gene_gtf["end"] >= chromEnd)]
        if gene.shape[0] > 0:
            df.at[row, "gene"] = gene["gene_name"].values[0]
            df.at[row, "gene_start"] = gene["start"].values[0]
            df.at[row, "gene_end"] = gene["end"].values[0]
        else:
            df.at[row, "gene"] = "intergenic"
df.to_csv(output)