# This is part of suite of scripts that analyses short read, next generation sequencing data for the human BRAF, KRAS, NRAS and PIK3CA genes generated from an Illumina Genome Analyser using the Qiagen "GeneRead DNAseq Targeted Panels V2 Human Tumor Actionable Mutations Panel" (Product no. 181900, Cat. no. NGHS-201X).  Please also see:
# ShortRead_quality_assessment.R
# BRAF_NGS_analysis_with_GenomicAlignments_and_Rsamtools_in_use_July_2019.R
# KRAS_NGS_analysis_with_GenomicAlignments_and_Rsamtools_in_use_July_2019.R
# NRAS_NGS_analysis_with_GenomicAlignments_and_Rsamtools_in_use_July_2019.R
# PIK3CA_NGS_analysis_with_GenomicAlignments_and_Rsamtools_in_use_July_2019.R

# This script loads the required R libraries and asks the data analyst to enter the names and locations of the files to be analysed and the desired name and location of the reports.  When these have been entered the sample name is extracted to be used for naming the reports.   

# Required libraries
library(Rsamtools)
library(dplyr)
library(tidyr)
library(reshape2)
library(ggplot2)

# INPUT FILES
# Path to folder containing the bam and bam.bai files.
bam_folder <- readline('Please enter the full path to the bam file location: ')

# Assign name of bam file.
bam_filename <- readline('Please enter the full bam file name: ')

# Assign name of bam.bai file.
bam.bai_filename <- readline('Please enter the full bam.bai file name: ')

# Extract the sample name from the file name by removing the suffix
bam_samplename <- unlist(strsplit(bam_filename, '\\.'))[1]
full_path_bam_file <- paste(bam_folder, bam_samplename, sep = '/')


# REPORTS
# Path for saving reports.
screening_report_folder <- readline('Please enter the desired name for the screening report folder: ')
path_to_report <- paste(bam_folder, screening_report_folder, sep = '/')
dir.create(path_to_report)

# BRAF
# Create names for percentage read reports 
BRAF_pc_reads_f <- paste(bam_samplename, '_percent_BRAF_fwd_reads.csv', sep = '')
BRAF_pc_reads_r <- paste(bam_samplename, '_percent_BRAF_rev_reads.csv', sep = '')

# Create names for stacked bar charts.
BRAF_chart_f <- paste(bam_samplename, '_BRAF_fwd_chart.pdf', sep = '')
BRAF_chart_r <- paste(bam_samplename, '_BRAF_rev_chart.pdf', sep = '')

# Percentage mutation reports.  Only one is required as the mutation report includes forward and reverse data.
BRAF_mutations <- paste(bam_samplename, '_BRAF_mutations.txt', sep = '')

# KRAS exon 2
# Create names for percentage read reports 
KRAS_2_pc_reads_f <- paste(bam_samplename, '_percent_KRAS_2_fwd_reads.csv', sep = '')
KRAS_2_pc_reads_r <- paste(bam_samplename, '_percent_KRAS_2_rev_reads.csv', sep = '')

# Create names for stacked bar charts.
KRAS_2_chart_f <- paste(bam_samplename, '_KRAS_2_fwd_chart.pdf', sep = '')
KRAS_2_chart_r <- paste(bam_samplename, '_KRAS_2_rev_chart.pdf', sep = '')

# Percentage mutation reports.  Only one is required as the mutation report includes forward and reverse data.
KRAS_2_mutations <- paste(bam_samplename, '_KRAS_2_fwd_mutations.txt', sep = '')

# KRAS exon 3
# Create names for percentage read reports 
KRAS_3_pc_reads_f <- paste(bam_samplename, '_percent_KRAS_3_fwd_reads.csv', sep = '')
KRAS_3_pc_reads_r <- paste(bam_samplename, '_percent_KRAS_3_rev_reads.csv', sep = '')

# Create names for stacked bar charts.
KRAS_3_chart_f <- paste(bam_samplename, '_KRAS_3_fwd_chart.pdf', sep = '')
KRAS_3_chart_r <- paste(bam_samplename, '_KRAS_3_rev_chart.pdf', sep = '')

# Percentage mutation reports.  Only one is required as the mutation report includes forward and reverse data.
KRAS_3_mutations <- paste(bam_samplename, '_KRAS_3_fwd_mutations.txt', sep = '')

# KRAS exon 4
# Create names for percentage read reports 
KRAS_4_pc_reads_f <- paste(bam_samplename, '_percent_KRAS_4_fwd_reads.csv', sep = '')
KRAS_4_pc_reads_r <- paste(bam_samplename, '_percent_KRAS_4_rev_reads.csv', sep = '')

# Create names for stacked bar charts.
KRAS_4_chart_f <- paste(bam_samplename, '_KRAS_4_fwd_chart.pdf', sep = '')
KRAS_4_chart_r <- paste(bam_samplename, '_KRAS_4_rev_chart.pdf', sep = '')

# Percentage mutation reports.  Only one is required as the mutation report includes forward and reverse data.
KRAS_4_mutations <- paste(bam_samplename, '_KRAS_4_fwd_mutations.txt', sep = '')

# NRAS exon 2
# Create names for percentage read reports 
NRAS_2_pc_reads_f <- paste(bam_samplename, '_percent_NRAS_2_fwd_reads.csv', sep = '')
NRAS_2_pc_reads_r <- paste(bam_samplename, '_percent_NRAS_2_rev_reads.csv', sep = '')

# Create names for stacked bar charts.
NRAS_2_chart_f <- paste(bam_samplename, '_NRAS_2_fwd_chart.pdf', sep = '')
NRAS_2_chart_r <- paste(bam_samplename, '_NRAS_2_rev_chart.pdf', sep = '')

# Percentage mutation reports.  Only one is required as the mutation report includes forward and reverse data.
NRAS_2_mutations <- paste(bam_samplename, '_NRAS_2_fwd_mutations.txt', sep = '')

# NRAS exon 3
# Create names for percentage read reports 
NRAS_3_pc_reads_f <- paste(bam_samplename, '_percent_NRAS_3_fwd_reads.csv', sep = '')
NRAS_3_pc_reads_r <- paste(bam_samplename, '_percent_NRAS_3_rev_reads.csv', sep = '')

# Create names for stacked bar charts.
NRAS_3_chart_f <- paste(bam_samplename, '_NRAS_3_fwd_chart.pdf', sep = '')
NRAS_3_chart_r <- paste(bam_samplename, '_NRAS_3_rev_chart.pdf', sep = '')

# Percentage mutation reports.  Only one is required as the mutation report includes forward and reverse data.
NRAS_3_mutations <- paste(bam_samplename, '_NRAS_3_fwd_mutations.txt', sep = '')

# PIK3CA exon 9
# Create names for percentage read reports 
PIK3CA_9_pc_reads_f <- paste(bam_samplename, '_percent_PIK3CA_9_fwd_reads.csv', sep = '')
PIK3CA_9_pc_reads_r <- paste(bam_samplename, '_percent_PIK3CA_9_rev_reads.csv', sep = '')

# Create names for stacked bar charts.
PIK3CA_9_chart_f <- paste(bam_samplename, '_PIK3CA_9_fwd_chart.pdf', sep = '')
PIK3CA_9_chart_r <- paste(bam_samplename, '_PIK3CA_9_rev_chart.pdf', sep = '')

# Percentage mutation reports.  Only one is required as the mutation report includes forward and reverse data.
PIK3CA_9_mutations <- paste(bam_samplename, '_PIK3CA_9_fwd_mutations.txt', sep = '')

# PIK3CA exon 21
# Create names for percentage read reports 
PIK3CA_21_pc_reads_f <- paste(bam_samplename, '_percent_PIK3CA_21_fwd_reads.csv', sep = '')
PIK3CA_21_pc_reads_r <- paste(bam_samplename, '_percent_PIK3CA_21_rev_reads.csv', sep = '')

# Create names for stacked bar charts.
PIK3CA_21_chart_f <- paste(bam_samplename, '_PIK3CA_21_fwd_chart.pdf', sep = '')
PIK3CA_21_chart_r <- paste(bam_samplename, '_PIK3CA_21_rev_chart.pdf', sep = '')

# Percentage mutation reports.  Only one is required as the mutation report includes forward and reverse data.
PIK3CA_21_mutations <- paste(bam_samplename, '_PIK3CA_21_fwd_mutations.txt', sep = '')

# Read in the bam and bam.bai files.  
sample_files <- BamFile(bam_filename, index = bam.bai_filename)
