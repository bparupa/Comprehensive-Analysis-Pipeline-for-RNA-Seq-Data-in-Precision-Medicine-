# Precision-medicine-assignment_2
Name: Bhavana Parupalli
Programming language: linux , R
Date: 11/13/2024

A
Set Up the Environment: Ensure that all the necessary modules are loaded using the module load command
command: module load sra-toolkit fastqc star/2.7.11a samtools/1.10 python
and created a directory named assignment_2

B
Data Download:
The script downloads the RNA-Seq data from the SRA database using the prefetch and fasterq-dump commands.
The accession IDs for the samples are stored in the ACCESSION_IDS array.
These files will be saved in the Raw_Data/ directory.
command:prefetch $SRA_ID -O $DATA_DIR
fasterq-dump $DATA_DIR/$SRA_ID -O $DATA_DIR

C
Quality Control with FastQC:
The FastQC tool runs on the raw FASTQ files to check for issues like adapter contamination, GC content, and sequence quality.
The reports are saved in the fastqc_reports/ directory.
command: fastqc $DATA_DIR/*.fastq -o $FASTQC_DIR

D
Trimming with Trim Galore:
Trim Galore is used to trim low-quality bases and adapter sequences from the reads.
The trimmed reads are saved in the N/u/bparupa/Quartz/assignment_2/trimmed_reads/ directory.
These file is saved as trimmed_reads
Command: trim_galore --quality 20 --length 30 -o $TRIM_DIR $DATA_DIR/*.fastq

E
Quality Control on Trimmed Reads:
A second round of FastQC is performed on the trimmed reads to ensure that trimming was effective.
The reports are saved in the fastqc_reports/trimmed directory.
Command:fastqc $TRIM_DIR/*.fq -o $FASTQC_DIR/trimmed

F
Reference Genome Download and hisat Indexing:
Download the Reference Genome and GTF Annotation
Create a Directory for Storing Genome Files: This step will create a directory for storing the downloaded Homo sapiens genome files and the annotation GTF file.
command:HISAT2_GENOME_DIR="/N/u/bparupa/Assignment-2/hisat2_genomes"
mkdir -p $HISAT2_GENOME_DIR
Download the Genome: Use wget to download the Homo sapiens (hg38) genome in .tar.gz format from the official HISAT2 genome repository.
command:wget -P $HISAT2_GENOME_DIR https://genome-idx.s3.amazonaws.com/hg38.tar.gz
Extract the Files: Once the download is complete, extract the contents of the .tar.gz archive to the designated directory.
command:tar -xzf $HISAT2_GENOME_DIR/hg38.tar.gz -C $HISAT2_GENOME_DIR

G
Build the HISAT2 Index
After downloading and extracting the reference genome, I will index it using HISAT2:
Create a Directory for the Index Files: Create a directory for storing the generated HISAT2 index files.
Command:HISAT2_INDEX_DIR="/N/u/bparupa/Quartz/Assignment-2/hisat2_index"
mkdir -p $HISAT2_INDEX_DIR

H
Output Files
After running the hisat2-build command, I had the following outputs:
HISAT2 Index Files: The files will be stored in the hisat2_index/ directory and will have the .ht2 extension 
Genome and Annotation Files: The genome FASTA (genome.fa) and GTF (genes.gtf) files will also be  available in the hisat2_genomes/ directory.

I
Quantification with HTSeq:
HTSeq is used to count the number of reads mapping to each gene.
The results are saved in the counts/ directory.
command:htseq-count -f bam -s no "$bam_file" $STAR_INDEX_DIR/Homo_sapiens/UCSC/hg38/Annotation/Genes/genes.gtf > $counts_DIR/${sample_name}_counts.txt
 Gene expression counts: For each gene, the number of reads (mapped to that gene) across all samples will appear
These files are stored in the counts/ directory. They are used for downstream analysis, including differential gene expression analysis.

J
Differential gene expression analysis
This script performs the following steps for RNA-Seq data analysis:
Differential Expression Analysis (DEG): Using DESeq2, it compares two conditions (Control vs SARS-CoV-2) and two time points (24h vs 72h).
Gene Ontology (GO) Enrichment Analysis: Using clusterProfiler, it identifies enriched biological processes associated with differentially expressed genes (DEGs).
Install the required packages:
command:install.packages("BiocManager")
BiocManager::install("DESeq2")
BiocManager::install("clusterProfiler")
BiocManager::install("org.Hs.eg.db")  # Human gene annotation database
install.packages("openxlsx")

K
Data Input:
Count Data:RNA-Seq counts matrix in counts1.txt (rows are genes, columns are samples).
           Read the count data into R:
  command:counts_data <- read.table("C:/Users/bhava/downloads/counts1.txt", header = TRUE, row.names = 1)
rawdata 
rawdata file (assignment_2_info.csv) that includes information about the samples (e.g., condition, time point).
    Read the metadata into R:
command:meta_data <- read.csv("C:/Users/bhava/downloads/assignment_2_info.csv")

L
Create DESeq2 Object: The DESeqDataSetFromMatrix function creates a DESeq2 object, specifying the count data and rawdata:
command:dds <- DESeqDataSetFromMatrix(countData = counts_data, 
                              colData = meta_data, 
                              design = ~ Condition + Timepoint)
Differential Expression Analysis: The DESeq function runs the differential expression analysis. Two contrasts are tested:
Control vs SARS-CoV-2:
command:res_condition <- results(dds, contrast = c("Condition", "Control", "SARS-CoV-2"))
sig_res_condition <- res_condition[which(res_condition$padj < 0.05), ]
write.csv(as.data.frame(sig_res_condition), "DEGs_Control_vs_SARS-CoV-2.csv1")
24h vs 72h:
command:res_timepoint <- results(dds, contrast = c("Timepoint", "24", "72"))
sig_res_timepoint <- res_timepoint[which(res_timepoint$padj < 0.05), ]
write.csv(as.data.frame(sig_res_timepoint), "DEGs_24h_vs_72h.csv2")

M
Gene Ontology (GO) Enrichment Analysis: Using clusterProfiler, GO term enrichment analysis is performed for the DEGs identified in the two comparisons.
GO Enrichment for Control vs SARS-CoV-2:
command:deg_genes_condition <- rownames(sig_res_condition)
go_condition <- enrichGO(gene = deg_genes_condition,
                         OrgDb = org.Hs.eg.db,
                         keyType = "ENSEMBL", 
                         ont = "BP",  # Biological Process GO terms
                         pvalueCutoff = 0.05)
write.csv(as.data.frame(go_condition), "GO_Enrichment_Control_vs_SARS-CoV-2.csv1")
GO Enrichment for 24h vs 72h:
command:deg_genes_timepoint <- rownames(sig_res_timepoint)
go_timepoint <- enrichGO(gene = deg_genes_timepoint,
                         OrgDb = org.Hs.eg.db,
                         keyType = "ENSEMBL", 
                         ont = "BP",  # Biological Process GO terms
                         pvalueCutoff = 0.05)
write.csv(as.data.frame(go_timepoint), "GO_Enrichment_24h_vs_72h.csv2")

Output Files:

  Differential Expression Results (DEGs):
        DEGs_Control_vs_SARS-CoV-2.csv1: Contains DEGs for the comparison between Control and SARS-CoV-2.
        DEGs_24h_vs_72h.csv2: Contains DEGs for the comparison between 24h and 72h.

  GO Enrichment Results:
        GO_Enrichment_Control_vs_SARS-CoV-2.csv1: GO term enrichment results for the Control vs SARS-CoV-2 DEGs.
        GO_Enrichment_24h_vs_72h.csv2: GO term enrichment results for the 24h vs 72h DEGs

