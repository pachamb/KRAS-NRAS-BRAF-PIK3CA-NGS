# Author: Phil Chambers
# Date of last edit: 24/06/2019

# Purpose: 
# This analyses short read, next generation sequencing data for the human PIK3CA gene generated from an Illumina Genome Analyser using the Qiagen "GeneRead DNAseq Targeted Panels V2 Human Tumor Actionable Mutations Panel" (Product no. 181900, Cat. no. NGHS-201X).  It should also work with data generated from an Illumina Genome Analyser using and any PCR primers that amplify the regions of interest. Gene alignment and numbering is for human genome build hg19.
# It generates three outputs:
# 1. csv files with the percentage reads for each nucleotide at each position.
# 2. pdf files of stacked bar charts showing the percentage of each nucleotide at each position.
# 3. A txt file containing any mutations detected.


# PIK3CA EXONS 6 AND 21

# bam file scanning parameters.  These are for PIK3CA exons 9 and 21.
pik3ca.all.param <- ScanBamParam(which = GRanges("chr3", IRanges(start = c(178935998, 178951882), end = c(178936122, 178952152))))

# Pileup parameters.
pik3ca.all.pup.param <- PileupParam(max_depth = 10000, min_base_quality = 20, min_nucleotide_depth = 10, distinguish_strands = T)

# Perform the pileup operation
pik3ca.all.pup <- pileup(sample_files, scanBamParam = pik3ca.all.param, pileupParam = pik3ca.all.pup.param)

# Remove column called 'which_label' and keep the rest
pik3ca.all.pup.2 <- select(pik3ca.all.pup, -which_label)

# Split the data into 2 data.frames for the 2 exons.
pik3ca.9.pup <- filter(pik3ca.all.pup.2, pos >= 178935998, pos <= 178936122)
pik3ca.21.pup <- filter(pik3ca.all.pup.2, pos >= 178951882, pos <= 178952152)

# PIK3CA EXON 9/CODONS 542, 545+546

# Make two dataframes one for the '+' strand and one for the '-' strand.
pik3ca.9.pup.f <- filter(pik3ca.9.pup, strand == '+')
pik3ca.9.pup.r <- filter(pik3ca.9.pup, strand == '-')

# If there are >=2 nucleotides at one position the pileup function generates >=2 rows for that position.  This is particularly useless for mutation screening where the whole aim is to detect positions where there is more than 1 nucleotide...
# spread (tidyr) generates a dataframe that looks very like a consensus matrix.  spread (tidyr) needs to be told the dataframe it's working on, a 'key', ie name of the column whose values will be used as column headings and a 'value', ie name of the column whose values will populate the cells. 'fill' fills in the missing values.  Most importantly, it combines data for one position in to 1 column.
pik3ca.9.pup.f.s <- spread(pik3ca.9.pup.f, pos, count, fill = 0)
pik3ca.9.pup.r.s <- spread(pik3ca.9.pup.r, pos, count, fill = 0)

# Sum the counts in columns of the above dataframes.
pik3ca.9.pup.f.s.sum <- apply(pik3ca.9.pup.f.s[,4:128], 2, sum)
pik3ca.9.pup.r.s.sum <- apply(pik3ca.9.pup.r.s[,4:128], 2, sum)

# rbind results of apply operation to 'consensus matrix'.  This will lose the seqnames, strand and nucleotide columns.
pik3ca.9.pup.2.f.s.2 <- rbind(pik3ca.9.pup.f.s[,4:128], pik3ca.9.pup.f.s.sum)
pik3ca.9.pup.2.r.s.2 <- rbind(pik3ca.9.pup.r.s[,4:128], pik3ca.9.pup.r.s.sum)

# Calculate the proportion of counts for each nucleotide. 
pik3ca.9.a_pc.reads.f <- pik3ca.9.pup.2.f.s.2[1,]/pik3ca.9.pup.2.f.s.2[nrow(pik3ca.9.pup.2.f.s.2),]
pik3ca.9.c_pc.reads.f <- pik3ca.9.pup.2.f.s.2[2,]/pik3ca.9.pup.2.f.s.2[nrow(pik3ca.9.pup.2.f.s.2),]
pik3ca.9.g_pc.reads.f <- pik3ca.9.pup.2.f.s.2[3,]/pik3ca.9.pup.2.f.s.2[nrow(pik3ca.9.pup.2.f.s.2),]
pik3ca.9.t_pc.reads.f <- pik3ca.9.pup.2.f.s.2[4,]/pik3ca.9.pup.2.f.s.2[nrow(pik3ca.9.pup.2.f.s.2),]
pik3ca.9.a_pc.reads.r <- pik3ca.9.pup.2.r.s.2[1,]/pik3ca.9.pup.2.r.s.2[nrow(pik3ca.9.pup.2.r.s.2),]
pik3ca.9.c_pc.reads.r <- pik3ca.9.pup.2.r.s.2[2,]/pik3ca.9.pup.2.r.s.2[nrow(pik3ca.9.pup.2.r.s.2),]
pik3ca.9.g_pc.reads.r <- pik3ca.9.pup.2.r.s.2[3,]/pik3ca.9.pup.2.r.s.2[nrow(pik3ca.9.pup.2.r.s.2),]
pik3ca.9.t_pc.reads.r <- pik3ca.9.pup.2.r.s.2[4,]/pik3ca.9.pup.2.r.s.2[nrow(pik3ca.9.pup.2.r.s.2),]


# Make a single dataframe consisting of proportions for all nucleotides and convert in to percentages.
pik3ca.9.percent.reads.f <- rbind(pik3ca.9.a_pc.reads.f, pik3ca.9.c_pc.reads.f, pik3ca.9.g_pc.reads.f, pik3ca.9.t_pc.reads.f)
pik3ca.9.percent.reads.r <- rbind(pik3ca.9.a_pc.reads.r, pik3ca.9.c_pc.reads.r, pik3ca.9.g_pc.reads.r, pik3ca.9.t_pc.reads.r)

# The column names of 'pik3ca.9.percent.reads.f' and 'pik3ca.9.percent.reads.r' are gDNA co-ordinates.  This needs to be changed to cDNA positions.  This is done by making a dataframe of cDNA positions and assigning that as new column names.
# Make a dataframe of cDNA positions for PIK3CA exon 9. 
pik3ca.9.cDNA <- as.data.frame(1540:1664)
pik3ca.9.cDNA <- t(pik3ca.9.cDNA)

# Assign the the above column names.
colnames(pik3ca.9.percent.reads.f) <- pik3ca.9.cDNA
colnames(pik3ca.9.percent.reads.r) <- pik3ca.9.cDNA

# Use 'apply' to convert proportions for each nucleotide in to percentages.
pik3ca.9.percent.reads.f <- apply(pik3ca.9.percent.reads.f, 1:2, function(x) round(as.numeric(x*100),2))
pik3ca.9.percent.reads.r <- apply(pik3ca.9.percent.reads.r, 1:2, function(x) round(as.numeric(x*100),2))

# Make a dataframe of 'seqnames', 'strand' and 'nucleotide' and attach to the percent reads dataframe.
# Forward
pik3ca.9.percent.reads.f.info <- pik3ca.9.pup.f.s[1:4,1:3]
pik3ca.9.percent.reads.f.export <- cbind(pik3ca.9.percent.reads.f.info, pik3ca.9.percent.reads.f)
# As PIK3CA is in the reverse orientation the row names have been reversed too.  However, the 'Nucleotide' column is not reversed
row.names(pik3ca.9.percent.reads.f.export) <- c('T', 'G', 'C', 'A')
# Reverse
pik3ca.9.percent.reads.r.info <- pik3ca.9.pup.r.s[1:4,1:3]
pik3ca.9.percent.reads.r.export <- cbind(pik3ca.9.percent.reads.r.info, pik3ca.9.percent.reads.r)
# As PIK3CA is in the reverse orientation the row names have been reversed too.  However, the 'Nucleotide' column is not reversed
row.names(pik3ca.9.percent.reads.r.export) <- c('T', 'G', 'C', 'A')


# Output percentage reads to a .csv file.
write.csv(pik3ca.9.percent.reads.f.export, paste(path_to_report, PIK3CA_9_pc_reads_f, sep = '/') quote = F)
write.csv(pik3ca.9.percent.reads.r.export, paste(path_to_report, PIK3CA_9_pc_reads_r, sep = '/') quote = F)

# This section rearranges the data into a form suitable for making a stacked bar chart with ggplot2

# melt (reshape2) converts 'wide' data into 'long' data, which is needed for making a stacked bar chart.
# Forward
rownames(pik3ca.9.percent.reads.f) <- c('T', 'G', 'C', 'A')
pik3ca.9.percent.reads.f.melt <- melt(pik3ca.9.percent.reads.f, id.vars = colnames(pik3ca.9.percent.reads.f), value.name = 'Percent')
colnames(pik3ca.9.percent.reads.f.melt) <- c('Nucleotide', 'cDNA_Position', 'Percent')

# Re-order the columns.  This makes the mutation report make more sense.
pik3ca.9.percent.reads.f.melt <- pik3ca.9.percent.reads.f.melt[c('cDNA_Position', 'Nucleotide', 'Percent')]

# Reverse
rownames(pik3ca.9.percent.reads.r) <- c('T', 'G', 'C', 'A')
pik3ca.9.percent.reads.r.melt <- melt(pik3ca.9.percent.reads.r, id.vars = colnames(pik3ca.9.percent.reads.r), value.name = 'Percent')
colnames(pik3ca.9.percent.reads.r.melt) <- c('Nucleotide', 'cDNA_Position', 'Percent')

# Re-order the columns.  This makes the mutation report make more sense.
pik3ca.9.percent.reads.r.melt <- pik3ca.9.percent.reads.r.melt[c('cDNA_Position', 'Nucleotide', 'Percent')]


# First, assign some colours to a palette.  See 'http://www.hypergurl.com/rgbcolorchart.html' for hexadecimal codes.
seq.data.colours <- c('#FF0000', '#FF9900', '#0000FF', '#33CC00')

# STACKED BARCHARTS: FORWARD
# Use 'melted' data to make the barcharts. Start by adding a 'c.' to the cDNA number.  Each cDNA number will appear in the column 4 times, making the number of factors in the column = (number of rows)/4.  This makes the data suitable for a stacked barchart.
pik3ca.9.percent.reads.f.melt.1 <- as.data.frame(paste('c.', pik3ca.9.percent.reads.f.melt$cDNA_Position, sep = ''))
pik3ca.9.percent.reads.f.melt.2 <- cbind(pik3ca.9.percent.reads.f.melt.1, pik3ca.9.percent.reads.f.melt)
colnames(pik3ca.9.percent.reads.f.melt.2) <- c('cDNA_Position', 'Position', 'Nucleotide', 'Percent')

# Make a stacked barchart.  The plot has been assigned to an object to prevent it being printed.
pik3ca.9.fwd.plot <- ggplot(pik3ca.9.percent.reads.f.melt.2, aes(x = cDNA_Position, y = Percent, fill = Nucleotide)) + geom_bar(stat = 'identity') + scale_fill_manual(values = seq.data.colours) + labs(title = paste('Sample number', bam_filename, '\n\nStacked bar chart for forward reads in PIK3CA exon 9', sep = ' '), x = 'cDNA position', y = 'Percentage allele') + theme(axis.text.x = element_text(angle = 90))
ggsave(paste(path_to_report, PIK3CA_9_chart_f, sep = '/'))

# STACKED BARCHARTS: REVERSE
# Use 'melted' data to make the barcharts. Start by adding a 'c.' to the cDNA number.  Each cDNA number will appear in the column 4 times, making the number of factors in the column = (number of rows)/4.  This makes the data suitable for a stacked barchart.
pik3ca.9.percent.reads.r.melt.1 <- as.data.frame(paste('c.', pik3ca.9.percent.reads.r.melt$cDNA_Position, sep = ''))
pik3ca.9.percent.reads.r.melt.2 <- cbind(pik3ca.9.percent.reads.r.melt.1, pik3ca.9.percent.reads.r.melt)
colnames(pik3ca.9.percent.reads.r.melt.2) <- c('cDNA_Position', 'Position', 'Nucleotide', 'Percent')

# Make a stacked barchart.
pik3ca.9.rev.plot <- ggplot(pik3ca.9.percent.reads.r.melt.2, aes(x = cDNA_Position, y = Percent, fill = Nucleotide)) + geom_bar(stat = 'identity') + scale_fill_manual(values = seq.data.colours) + labs(title = paste('Sample number', bam_filename, '\n\nStacked bar chart for reverse reads in PIK3CA exon 9', sep = ' '), x = 'cDNA position', y = 'Percentage allele') + theme(axis.text.x = element_text(angle = 90))
ggsave(paste(path_to_report, PIK3CA_9_chart_r, sep = '/'))


# Mutation report: FORWARD 
# Add a 'c.' and the nucleotide to the melted data.  This means that the number of rows in the column = number of factors.  This makes the data suitable for a properly annotated mutation report. 
pik3ca.9.percent.reads.f.melt.3 <- as.data.frame(paste('c.', pik3ca.9.percent.reads.f.melt$cDNA_Position, pik3ca.9.percent.reads.f.melt$Nucleotide, sep = ''))
pik3ca.9.percent.reads.f.melt.4 <- cbind(pik3ca.9.percent.reads.f.melt.3, pik3ca.9.percent.reads.f.melt)
colnames(pik3ca.9.percent.reads.f.melt.4) <- c('cDNA_Position', 'Position', 'Nucleotide', 'Percent')

# Filter out the codons of interest plus one either side.  Filter by the 'Position' column.
pik3ca.9.percent.reads.f.melt.541_546 <- filter(pik3ca.9.percent.reads.f.melt.4, Position >= 31, Position <= 42)

# Mutation report: REVERSE
# Add a 'c.' and the nucleotide to the melted data.  This means that the number of rows in the column = number of factors.  This makes the data suitable for a properly annotated mutation report.  
pik3ca.9.percent.reads.r.melt.3 <- as.data.frame(paste('c.', pik3ca.9.percent.reads.r.melt$cDNA_Position, pik3ca.9.percent.reads.r.melt$Nucleotide, sep = ''))
pik3ca.9.percent.reads.r.melt.4 <- cbind(pik3ca.9.percent.reads.r.melt.3, pik3ca.9.percent.reads.r.melt)
colnames(pik3ca.9.percent.reads.r.melt.4) <- c('cDNA_Position', 'Position', 'Nucleotide', 'Percent')

# Filter out the codons of interest plus one either side.  Filter by the 'Position' column.
pik3ca.9.percent.reads.r.melt.541_546 <- filter(pik3ca.9.percent.reads.r.melt.4, Position >= 31, Position <= 42)

# Return positions that have nucleotides >=5% and <=95%.  This should detect nucleotide changes greater than 5%.
pik3ca.9.muts_5_f <- filter(pik3ca.9.percent.reads.f.melt, Percent >= 5.00, Percent <= 95.00); pik3ca.9.muts_5_f
pik3ca.9.muts_5_r <- filter(pik3ca.9.percent.reads.r.melt, Percent >= 5.00, Percent <= 95.00); pik3ca.9.muts_5_r

# Export the data to a mutation report.
write.table(paste('Sample number', bam_samplename, '\n\n', sep = ' '), paste(path_to_report, PIK3CA_9_mutations, sep = '/'), row.names = F, col.names = F, quote = F, append = T)
write.table(pik3ca.9.muts_5_f, paste(path_to_report, PIK3CA_9_mutations, sep = '/'), row.names = F, quote = F, append = T)
write.table('\n', paste(path_to_report, PIK3CA_9_mutations, sep = '/'), row.names = F, col.names = F, quote = F, append = T)
write.table(pik3ca.9.muts_5_r, paste(path_to_report, PIK3CA_9_mutations, sep = '/'), row.names = F, quote = F, append = T)

#=========================================================================================================

# PIK3CA EXON 21/CODONS 542, 545+546

# Make two dataframes one for the '+' strand and one for the '-' strand.
pik3ca.21.pup.f <- filter(pik3ca.21.pup, strand == '+')
pik3ca.21.pup.r <- filter(pik3ca.21.pup, strand == '-')

# If there are >=2 nucleotides at one position the pileup function generates >=2 rows for that position.  This is particularly useless for mutation screening where the whole aim is to detect positions where there is more than 1 nucleotide...
# spread (tidyr) generates a dataframe that looks very like a consensus matrix.  spread (tidyr) needs to be told the dataframe it's working on, a 'key', ie name of the column whose values will be used as column headings and a 'value', ie name of the column whose values will populate the cells. 'fill' fills in the missing values.  Most importantly, it combines data for one position in to 1 column.
pik3ca.21.pup.f.s <- spread(pik3ca.21.pup.f, pos, count, fill = 0)
pik3ca.21.pup.r.s <- spread(pik3ca.21.pup.r, pos, count, fill = 0)

# Sum the counts in columns of the above dataframes.
pik3ca.21.pup.f.s.sum <- apply(pik3ca.21.pup.f.s[,4:274], 2, sum)
pik3ca.21.pup.r.s.sum <- apply(pik3ca.21.pup.r.s[,4:274], 2, sum)

# rbind results of apply operation to 'consensus matrix'.  This will lose the seqnames, strand and nucleotide columns.
pik3ca.21.pup.2.f.s.2 <- rbind(pik3ca.21.pup.f.s[,4:274], pik3ca.21.pup.f.s.sum)
pik3ca.21.pup.2.r.s.2 <- rbind(pik3ca.21.pup.r.s[,4:274], pik3ca.21.pup.r.s.sum)

# Calculate the proportion of counts for each nucleotide. 
pik3ca.21.a_pc.reads.f <- pik3ca.21.pup.2.f.s.2[1,]/pik3ca.21.pup.2.f.s.2[nrow(pik3ca.21.pup.2.f.s.2),]
pik3ca.21.c_pc.reads.f <- pik3ca.21.pup.2.f.s.2[2,]/pik3ca.21.pup.2.f.s.2[nrow(pik3ca.21.pup.2.f.s.2),]
pik3ca.21.g_pc.reads.f <- pik3ca.21.pup.2.f.s.2[3,]/pik3ca.21.pup.2.f.s.2[nrow(pik3ca.21.pup.2.f.s.2),]
pik3ca.21.t_pc.reads.f <- pik3ca.21.pup.2.f.s.2[4,]/pik3ca.21.pup.2.f.s.2[nrow(pik3ca.21.pup.2.f.s.2),]
pik3ca.21.a_pc.reads.r <- pik3ca.21.pup.2.r.s.2[1,]/pik3ca.21.pup.2.r.s.2[nrow(pik3ca.21.pup.2.r.s.2),]
pik3ca.21.c_pc.reads.r <- pik3ca.21.pup.2.r.s.2[2,]/pik3ca.21.pup.2.r.s.2[nrow(pik3ca.21.pup.2.r.s.2),]
pik3ca.21.g_pc.reads.r <- pik3ca.21.pup.2.r.s.2[3,]/pik3ca.21.pup.2.r.s.2[nrow(pik3ca.21.pup.2.r.s.2),]
pik3ca.21.t_pc.reads.r <- pik3ca.21.pup.2.r.s.2[4,]/pik3ca.21.pup.2.r.s.2[nrow(pik3ca.21.pup.2.r.s.2),]


# Make a single dataframe consisting of proportions for all nucleotides and convert in to percentages.
pik3ca.21.percent.reads.f <- rbind(pik3ca.21.a_pc.reads.f, pik3ca.21.c_pc.reads.f, pik3ca.21.g_pc.reads.f, pik3ca.21.t_pc.reads.f)
pik3ca.21.percent.reads.r <- rbind(pik3ca.21.a_pc.reads.r, pik3ca.21.c_pc.reads.r, pik3ca.21.g_pc.reads.r, pik3ca.21.t_pc.reads.r)

# The column names of 'pik3ca.21.percent.reads.f' and 'pik3ca.21.percent.reads.r' are gDNA co-ordinates.  This needs to be changed to cDNA positions.  This is done by making a dataframe of cDNA positions and assigning that as new column names.
# Make a dataframe of cDNA positions for PIK3CA exon 21. 
pik3ca.21.cDNA <- as.data.frame(2937:3207)
pik3ca.21.cDNA <- t(pik3ca.21.cDNA)

# Assign the the above column names.
colnames(pik3ca.21.percent.reads.f) <- pik3ca.21.cDNA
colnames(pik3ca.21.percent.reads.r) <- pik3ca.21.cDNA

# Use 'apply' to convert proportions for each nucleotide in to percentages.
pik3ca.21.percent.reads.f <- apply(pik3ca.21.percent.reads.f, 1:2, function(x) round(as.numeric(x*100),2))
pik3ca.21.percent.reads.r <- apply(pik3ca.21.percent.reads.r, 1:2, function(x) round(as.numeric(x*100),2))

# Make a dataframe of 'seqnames', 'strand' and 'nucleotide' and attach to the percent reads dataframe.
# Forward
pik3ca.21.percent.reads.f.info <- pik3ca.21.pup.f.s[1:4,1:3]
pik3ca.21.percent.reads.f.export <- cbind(pik3ca.21.percent.reads.f.info, pik3ca.21.percent.reads.f)
# As PIK3CA is in the reverse orientation the row names have been reversed too.  However, the 'Nucleotide' column is not reversed
row.names(pik3ca.21.percent.reads.f.export) <- c('T', 'G', 'C', 'A')
# Reverse
pik3ca.21.percent.reads.r.info <- pik3ca.21.pup.r.s[1:4,1:3]
pik3ca.21.percent.reads.r.export <- cbind(pik3ca.21.percent.reads.r.info, pik3ca.21.percent.reads.r)
# As PIK3CA is in the reverse orientation the row names have been reversed too.  However, the 'Nucleotide' column is not reversed
row.names(pik3ca.21.percent.reads.r.export) <- c('T', 'G', 'C', 'A')


# Output percentage reads to a .csv file.
write.csv(pik3ca.21.percent.reads.f.export, paste(path_to_report, PIK3CA_21_pc_reads_f, sep = '/') quote = F)
write.csv(pik3ca.21.percent.reads.r.export, paste(path_to_report, PIK3CA_21_pc_reads_r, sep = '/') quote = F)

# This section rearranges the data into a form suitable for making a stacked bar chart with ggplot2

# melt (reshape2) converts 'wide' data into 'long' data, which is needed for making a stacked bar chart.
# Forward
rownames(pik3ca.21.percent.reads.f) <- c('T', 'G', 'C', 'A')
pik3ca.21.percent.reads.f.melt <- melt(pik3ca.21.percent.reads.f, id.vars = colnames(pik3ca.21.percent.reads.f), value.name = 'Percent')
colnames(pik3ca.21.percent.reads.f.melt) <- c('Nucleotide', 'cDNA_Position', 'Percent')

# Re-order the columns.  This makes the mutation report make more sense.
pik3ca.21.percent.reads.f.melt <- pik3ca.21.percent.reads.f.melt[c('cDNA_Position', 'Nucleotide', 'Percent')]

# Reverse
rownames(pik3ca.21.percent.reads.r) <- c('T', 'G', 'C', 'A')
pik3ca.21.percent.reads.r.melt <- melt(pik3ca.21.percent.reads.r, id.vars = colnames(pik3ca.21.percent.reads.r), value.name = 'Percent')
colnames(pik3ca.21.percent.reads.r.melt) <- c('Nucleotide', 'cDNA_Position', 'Percent')

# Re-order the columns.  This makes the mutation report make more sense.
pik3ca.21.percent.reads.r.melt <- pik3ca.21.percent.reads.r.melt[c('cDNA_Position', 'Nucleotide', 'Percent')]


# First, assign some colours to a palette.  See 'http://www.hypergurl.com/rgbcolorchart.html' for hexadecimal codes.
seq.data.colours <- c('#FF0000', '#FF9900', '#0000FF', '#33CC00')

# STACKED BARCHARTS: FORWARD
# Use 'melted' data to make the barcharts. Start by adding a 'c.' to the cDNA number.  Each cDNA number will appear in the column 4 times, making the number of factors in the column = (number of rows)/4.  This makes the data suitable for a stacked barchart.
pik3ca.21.percent.reads.f.melt.1 <- as.data.frame(paste('c.', pik3ca.21.percent.reads.f.melt$cDNA_Position, sep = ''))
pik3ca.21.percent.reads.f.melt.2 <- cbind(pik3ca.21.percent.reads.f.melt.1, pik3ca.21.percent.reads.f.melt)
colnames(pik3ca.21.percent.reads.f.melt.2) <- c('cDNA_Position', 'Position', 'Nucleotide', 'Percent')

# Make a stacked barchart.  The plot has been assigned to an object to prevent it being printed.
pik3ca.21.fwd.plot <- ggplot(pik3ca.21.percent.reads.f.melt.2, aes(x = cDNA_Position, y = Percent, fill = Nucleotide)) + geom_bar(stat = 'identity') + scale_fill_manual(values = seq.data.colours) + labs(title = paste('Sample number', bam_filename, '\n\nStacked bar chart for forward reads in PIK3CA exon 21', sep = ' '), x = 'cDNA position', y = 'Percentage allele') + theme(axis.text.x = element_text(angle = 90))
ggsave(paste(path_to_report, PIK3CA_21_chart_f, sep = '/'))

# STACKED BARCHARTS: REVERSE
# Use 'melted' data to make the barcharts. Start by adding a 'c.' to the cDNA number.  Each cDNA number will appear in the column 4 times, making the number of factors in the column = (number of rows)/4.  This makes the data suitable for a stacked barchart.
pik3ca.21.percent.reads.r.melt.1 <- as.data.frame(paste('c.', pik3ca.21.percent.reads.r.melt$cDNA_Position, sep = ''))
pik3ca.21.percent.reads.r.melt.2 <- cbind(pik3ca.21.percent.reads.r.melt.1, pik3ca.21.percent.reads.r.melt)
colnames(pik3ca.21.percent.reads.r.melt.2) <- c('cDNA_Position', 'Position', 'Nucleotide', 'Percent')

# Make a stacked barchart.
pik3ca.21.rev.plot <- ggplot(pik3ca.21.percent.reads.r.melt.2, aes(x = cDNA_Position, y = Percent, fill = Nucleotide)) + geom_bar(stat = 'identity') + scale_fill_manual(values = seq.data.colours) + labs(title = paste('Sample number', bam_filename, '\n\nStacked bar chart for reverse reads in PIK3CA exon 21', sep = ' '), x = 'cDNA position', y = 'Percentage allele') + theme(axis.text.x = element_text(angle = 90))
ggsave(paste(path_to_report, PIK3CA_21_chart_r, sep = '/'))


# Mutation report: FORWARD 
# Add a 'c.' and the nucleotide to the melted data.  This means that the number of rows in the column = number of factors.  This makes the data suitable for a properly annotated mutation report. 
pik3ca.21.percent.reads.f.melt.3 <- as.data.frame(paste('c.', pik3ca.21.percent.reads.f.melt$cDNA_Position, pik3ca.21.percent.reads.f.melt$Nucleotide, sep = ''))
pik3ca.21.percent.reads.f.melt.4 <- cbind(pik3ca.21.percent.reads.f.melt.3, pik3ca.21.percent.reads.f.melt)
colnames(pik3ca.21.percent.reads.f.melt.4) <- c('cDNA_Position', 'Position', 'Nucleotide', 'Percent')

# Filter out the codons of interest plus one either side.  Filter by the 'Position' column.
pik3ca.21.percent.reads.f.melt.541_546 <- filter(pik3ca.21.percent.reads.f.melt.4, Position >= 3136, Position <= 3144)

# Mutation report: REVERSE
# Add a 'c.' and the nucleotide to the melted data.  This means that the number of rows in the column = number of factors.  This makes the data suitable for a properly annotated mutation report.  
pik3ca.21.percent.reads.r.melt.3 <- as.data.frame(paste('c.', pik3ca.21.percent.reads.r.melt$cDNA_Position, pik3ca.21.percent.reads.r.melt$Nucleotide, sep = ''))
pik3ca.21.percent.reads.r.melt.4 <- cbind(pik3ca.21.percent.reads.r.melt.3, pik3ca.21.percent.reads.r.melt)
colnames(pik3ca.21.percent.reads.r.melt.4) <- c('cDNA_Position', 'Position', 'Nucleotide', 'Percent')

# Filter out the codons of interest plus one either side.  Filter by the 'Position' column.
pik3ca.21.percent.reads.r.melt.541_546 <- filter(pik3ca.21.percent.reads.r.melt.4, Position >= 3136, Position <= 3144)

# Return positions that have nucleotides >=5% and <=95%.  This should detect nucleotide changes greater than 5%.
pik3ca.21.muts_5_f <- filter(pik3ca.21.percent.reads.f.melt, Percent >= 5.00, Percent <= 95.00); pik3ca.21.muts_5_f
pik3ca.21.muts_5_r <- filter(pik3ca.21.percent.reads.r.melt, Percent >= 5.00, Percent <= 95.00); pik3ca.21.muts_5_r

# Export the data to a mutation report.
write.table(paste('Sample number', bam_samplename, '\n\n', sep = ' '), paste(path_to_report, PIK3CA_21_mutations, sep = '/'), row.names = F, col.names = F, quote = F, append = T)
write.table(pik3ca.21.muts_5_f, paste(path_to_report, PIK3CA_21_mutations, sep = '/'), row.names = F, quote = F, append = T)
write.table('\n', paste(path_to_report, PIK3CA_21_mutations, sep = '/'), row.names = F, col.names = F, quote = F, append = T)
write.table(pik3ca.21.muts_5_r, paste(path_to_report, PIK3CA_21_mutations, sep = '/'), row.names = F, quote = F, append = T)

