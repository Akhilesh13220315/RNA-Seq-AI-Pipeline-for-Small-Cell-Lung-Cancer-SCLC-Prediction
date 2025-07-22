# RNA-Seq-AI-Pipeline-for-Small-Cell-Lung-Cancer-SCLC-Prediction
This repository outlines an end-to-end RNA-seq analysis pipeline for Small Cell Lung Cancer (SCLC) using public datasets. The pipeline includes data download, preprocessing, alignment, quantification, differential expression analysis, and functional annotation.

1. Dataset Selection
A widely used public dataset for SCLC RNA-seq is available from the Gene Expression Omnibus (GEO).

Example Dataset:

GEO Accession: 

https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE60052
Description: RNA-seq of 7 SCLC cell lines and 2 normal bronchial epithelial cell lines.
Direct FASTQ download:

Use the 

https://www.ncbi.nlm.nih.gov/Traces/study/?acc=GSE60052
 to get SRA/FASTQ files.
2. Pipeline Overview
Data Download
Quality Control
Trimming
Alignment
Quantification
Differential Expression Analysis
Functional Annotation
3. Step-by-Step Pipeline
Step 1: Data Download
bash



# Install SRA Toolkit if not already installed
conda install -c bioconda sra-tools

# Download SRA files (replace SRRxxxxxxx with actual run IDs from GSE60052)
prefetch SRRxxxxxxx

# Convert SRA to FASTQ
fastq-dump --split-files SRRxxxxxxx.sra
Step 2: Quality Control
bash



# Install FastQC
conda install -c bioconda fastqc

# Run FastQC
fastqc SRRxxxxxxx_1.fastq SRRxxxxxxx_2.fastq
Review the HTML reports for quality issues.
Step 3: Trimming (if needed)
bash



# Install Trimmomatic
conda install -c bioconda trimmomatic

# Run Trimmomatic for paired-end data
trimmomatic PE SRRxxxxxxx_1.fastq SRRxxxxxxx_2.fastq \
  SRRxxxxxxx_1_paired.fastq SRRxxxxxxx_1_unpaired.fastq \
  SRRxxxxxxx_2_paired.fastq SRRxxxxxxx_2_unpaired.fastq \
  ILLUMINACLIP:TruSeq3-PE.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36
Step 4: Alignment
Reference genome: 

https://www.gencodegenes.org/human/

bash



# Install STAR aligner
conda install -c bioconda star

# Generate genome index (only once)
STAR --runThreadN 8 --runMode genomeGenerate --genomeDir /path/to/genomeDir \
  --genomeFastaFiles GRCh38.primary_assembly.genome.fa \
  --sjdbGTFfile gencode.vXX.annotation.gtf

# Align reads
STAR --runThreadN 8 --genomeDir /path/to/genomeDir \
  --readFilesIn SRRxxxxxxx_1_paired.fastq SRRxxxxxxx_2_paired.fastq \
  --outFileNamePrefix SRRxxxxxxx_ --outSAMtype BAM SortedByCoordinate
Step 5: Quantification
bash



# Install featureCounts
conda install -c bioconda subread

# Run featureCounts
featureCounts -T 8 -a gencode.vXX.annotation.gtf -o counts.txt SRRxxxxxxx_Aligned.sortedByCoord.out.bam
Step 6: Differential Expression Analysis
Recommended tool: 

https://bioconductor.org/packages/release/bioc/html/DESeq2.html
 in R.

r



# In R
library(DESeq2)

# Read in count matrix and sample info
counts <- read.table("counts.txt", header=TRUE, row.names=1)
coldata <- data.frame(
  row.names = colnames(counts),
  condition = c(rep("SCLC", 7), rep("Normal", 2))
)

dds <- DESeqDataSetFromMatrix(countData = counts, colData = coldata, design = ~ condition)
dds <- DESeq(dds)
res <- results(dds)
write.csv(as.data.frame(res), file="DEG_results.csv")
Step 7: Functional Annotation
Gene Ontology/Pathway Analysis:
Use 

https://bioconductor.org/packages/release/bioc/html/clusterProfiler.html
 or 

https://maayanlab.cloud/Enrichr/
.

r



library(clusterProfiler)
library(org.Hs.eg.db)

# Get significant genes
sig_genes <- rownames(res)[which(res$padj < 0.05 & abs(res$log2FoldChange) > 1)]

# GO enrichment
ego <- enrichGO(gene = sig_genes, OrgDb = org.Hs.eg.db, keyType = "SYMBOL", ont = "BP")
dotplot(ego)
Summary Table
Step	Tool(s)	Command/Script Example
Data Download	SRA Toolkit	prefetch, fastq-dump
Quality Control	FastQC	fastqc
Trimming	Trimmomatic	trimmomatic
Alignment	STAR	STAR
Quantification	featureCounts	featureCounts
Differential Expr.	DESeq2 (R)	DESeq2 workflow
Functional Analysis	clusterProfiler	enrichGO, dotplot
References & Resources

https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE60052

https://github.com/ncbi/sra-tools

https://www.bioinformatics.babraham.ac.uk/projects/fastqc/

http://www.usadellab.org/cms/?page=trimmomatic

https://github.com/alexdobin/STAR

http://bioinf.wehi.edu.au/featureCounts/

https://bioconductor.org/packages/release/bioc/html/DESeq2.html

https://bioconductor.org/packages/release/bioc/html/clusterProfiler.html  + in proper fomrat no change
