# Quality assessment of short Illumina NGS reads.
# This assessment is meant to COMPLEMENT, rather than replace, quality assessment available from the Illumina Genome Analyzer.
# This code works on compressed fastq files, i.e. files named "[name].fastq.gz".  It accepts one file at a time.  Please get in touch if you require a version that accepts multiple files.


# Required libraries:
library(ShortRead)
library(htmltools)

# Direct the software to and read in ONE compressed fastq file:
sample_folder <- readline('Please enter the full path to the file location: ')
sample_filename <- readline('Please enter the full file name: ')
full_path_sample <- paste(sample_folder, sample_filename, sep = '/')

# Specify a path and location for saving a report:
report_folder <- readline('Please enter the desired name for the report folder: ')
full_path_report <- paste(sample_folder, report_folder, sep = '/')
dir.create(full_path_report)
report_name <- paste(report_folder, 'txt', sep = '.')

# Use ShortRead::readFastq to generate an object that stores the length (i.e. number of reads in the fastq file) and width (i.e. number of sequencing cycles) to a .txt file.
fq_reads <- readFastq(full_path_sample)

# Write a report containing the required data from the "readFastq" function
write.table(paste('Sample number', sample_filename, '\n\nThere were', length(fq_reads), 'short reads for sample number', sample_filename), paste(full_path_report, report_name, sep = '/'), quote = F, row.names = F, col.names = F, append = T)
write.table(paste('\nThe reads were', width(fq_reads[1]), 'bp long.'), paste(full_path_report, report_name, sep = '/'), quote = F, row.names = F, col.names = F, append = T)

# Create an object which is a quality summary for the one file.  This object is class FastqQA, which is unique to the ShortRead package.  This will generate and store an html file in the report folder.

readqa <- qa(full_path_sample, type = 'fastq')

report(readqa, dest = (full_path_report))
