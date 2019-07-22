# Author: Phil Chambers
# Date of last edit: 24/06/2019

# Purpose: 
# This analyses short read, next generation sequencing data for the human KRAS gene generated from an Illumina Genome Analyser using the Qiagen "GeneRead DNAseq Targeted Panels V2 Human Tumor Actionable Mutations Panel" (Product no. 181900, Cat. no. NGHS-201X).  It should also work with data generated from an Illumina Genome Analyser using and any PCR primers that amplify the regions of interest. Gene alignment and numbering is for human genome build hg19.
# It generates three outputs:
# 1. csv files with the percentage reads for each nucleotide at each position.
# 2. pdf files of stacked bar charts showing the percentage of each nucleotide at each position.
# 3. A txt file containing any mutations detected.


# KRAS EXONS 2, 3 AND 4

# bam file scanning parameters.  These are for KRAS exons 2, 3 and 4.
kras.all.param <- ScanBamParam(which = GRanges('chr12', IRanges(start = c(25398208, 25380168, 25378548), end = c(25398318, 25380346, 25378707))))

# Pileup parameters.
kras.all.pup.param <- PileupParam(max_depth = 10000, min_base_quality = 20, min_nucleotide_depth = 10, distinguish_strands = T)

# Perform the pileup operation
kras.all.pup <- pileup(sample_files, scanBamParam = kras.all.param, pileupParam = kras.all.pup.param)

# Remove column called 'which_label' and keep the rest
kras.all.pup.2 <- select(kras.all.pup, -which_label)

# Split the data into 3 dataframes for the 3 exons.
kras.2.pup <- filter(kras.all.pup.2, pos >= 25398208, pos <= 25398318)
kras.3.pup <- filter(kras.all.pup.2, pos >= 25380168, pos <= 25380346)
kras.4.pup <- filter(kras.all.pup.2, pos >= 25378548, pos <= 25378707)

# KRAS EXON 2/CODONS 12+13

# Make two dataframes one for the '+' strand and one for the '-' strand.
kras.2.pup.f <- filter(kras.2.pup, strand == '+')
kras.2.pup.r <- filter(kras.2.pup, strand == '-')

# If there are >=2 nucleotides at one position the pileup function generates >=2 rows for that position.  This is particularly useless for mutation screening where the whole aim is to detect positions where there is more than 1 nucleotide...
# spread (tidyr) generates a dataframe that looks very like a consensus matrix.  spread (tidyr) needs to be told the dataframe it's working on, a 'key', ie name of the column whose values will be used as column headings and a 'value', ie name of the column whose values will populate the cells. 'fill' fills in the missing values.  Most importantly, it combines data for one position in to 1 column.
kras.2.pup.f.s <- spread(kras.2.pup.f, pos, count, fill = 0)
kras.2.pup.r.s <- spread(kras.2.pup.r, pos, count, fill = 0)

# Sum the counts in columns of the above dataframes.
kras.2.pup.f.s.sum <- apply(kras.2.pup.f.s[,4:114], 2, sum)
kras.2.pup.r.s.sum <- apply(kras.2.pup.r.s[,4:114], 2, sum)

# rbind results of apply operation to 'consensus matrix'.  This will lose the seqnames, strand and nucleotide columns.
kras.2.pup.2.f.s.2 <- rbind(kras.2.pup.f.s[,4:114], kras.2.pup.f.s.sum)
kras.2.pup.2.r.s.2 <- rbind(kras.2.pup.r.s[,4:114], kras.2.pup.r.s.sum)

# Calculate the proportion of counts for each nucleotide. 
kras.2.a_pc.reads.f <- kras.2.pup.2.f.s.2[1,]/kras.2.pup.2.f.s.2[nrow(kras.2.pup.2.f.s.2),]
kras.2.c_pc.reads.f <- kras.2.pup.2.f.s.2[2,]/kras.2.pup.2.f.s.2[nrow(kras.2.pup.2.f.s.2),]
kras.2.g_pc.reads.f <- kras.2.pup.2.f.s.2[3,]/kras.2.pup.2.f.s.2[nrow(kras.2.pup.2.f.s.2),]
kras.2.t_pc.reads.f <- kras.2.pup.2.f.s.2[4,]/kras.2.pup.2.f.s.2[nrow(kras.2.pup.2.f.s.2),]
kras.2.a_pc.reads.r <- kras.2.pup.2.r.s.2[1,]/kras.2.pup.2.r.s.2[nrow(kras.2.pup.2.r.s.2),]
kras.2.c_pc.reads.r <- kras.2.pup.2.r.s.2[2,]/kras.2.pup.2.r.s.2[nrow(kras.2.pup.2.r.s.2),]
kras.2.g_pc.reads.r <- kras.2.pup.2.r.s.2[3,]/kras.2.pup.2.r.s.2[nrow(kras.2.pup.2.r.s.2),]
kras.2.t_pc.reads.r <- kras.2.pup.2.r.s.2[4,]/kras.2.pup.2.r.s.2[nrow(kras.2.pup.2.r.s.2),]


# Make a single dataframe consisting of proportions for all nucleotides and convert in to percentages.
kras.2.percent.reads.f <- rbind(kras.2.a_pc.reads.f, kras.2.c_pc.reads.f, kras.2.g_pc.reads.f, kras.2.t_pc.reads.f)
kras.2.percent.reads.r <- rbind(kras.2.a_pc.reads.r, kras.2.c_pc.reads.r, kras.2.g_pc.reads.r, kras.2.t_pc.reads.r)

# The column names of 'kras.2.percent.reads.f' and 'kras.2.percent.reads.r' are gDNA co-ordinates.  This needs to be changed to cDNA positions.  This is done by making a dataframe of cDNA positions and assigning that as new column names.
# Make a dataframe of cDNA positions for KRAS exon 2. 
kras.2.cDNA <- as.data.frame(111:1)
kras.2.cDNA <- t(kras.2.cDNA)

# Assign the the above column names.
colnames(kras.2.percent.reads.f) <- kras.2.cDNA
colnames(kras.2.percent.reads.r) <- kras.2.cDNA

# Use 'apply' to convert proportions for each nucleotide in to percentages.
kras.2.percent.reads.f <- apply(kras.2.percent.reads.f, 1:2, function(x) round(as.numeric(x*100),2))
kras.2.percent.reads.r <- apply(kras.2.percent.reads.r, 1:2, function(x) round(as.numeric(x*100),2))

# Make a dataframe of 'seqnames', 'strand' and 'nucleotide' and attach to the percent reads dataframe.
# Forward
kras.2.percent.reads.f.info <- kras.2.pup.f.s[1:4,1:3]
kras.2.percent.reads.f.export <- cbind(kras.2.percent.reads.f.info, kras.2.percent.reads.f)
# As KRAS is in the reverse orientation the row names have been reversed too.  However, the 'Nucleotide' column is not reversed
row.names(kras.2.percent.reads.f.export) <- c('T', 'G', 'C', 'A')
# Reverse
kras.2.percent.reads.r.info <- kras.2.pup.r.s[1:4,1:3]
kras.2.percent.reads.r.export <- cbind(kras.2.percent.reads.r.info, kras.2.percent.reads.r)
# As KRAS is in the reverse orientation the row names have been reversed too.  However, the 'Nucleotide' column is not reversed
row.names(kras.2.percent.reads.r.export) <- c('T', 'G', 'C', 'A')


# Output percentage reads to a .csv file.
write.csv(kras.2.percent.reads.f.export, paste(path_to_report, KRAS_2_pc_reads_f, sep = '/') quote = F)
write.csv(kras.2.percent.reads.r.export, paste(path_to_report, KRAS_2_pc_reads_r, sep = '/') quote = F)

# This section rearranges the data into a form suitable for making a stacked bar chart with ggplot2

# melt (reshape2) converts 'wide' data into 'long' data, which is needed for making a stacked bar chart.
# Forward
rownames(kras.2.percent.reads.f) <- c('T', 'G', 'C', 'A')
kras.2.percent.reads.f.melt <- melt(kras.2.percent.reads.f, id.vars = colnames(kras.2.percent.reads.f), value.name = 'Percent')
colnames(kras.2.percent.reads.f.melt) <- c('Nucleotide', 'cDNA_Position', 'Percent')

# Re-order the columns.  This makes the mutation report make more sense.
kras.2.percent.reads.f.melt <- kras.2.percent.reads.f.melt[c('cDNA_Position', 'Nucleotide', 'Percent')]

# Reverse
rownames(kras.2.percent.reads.r) <- c('T', 'G', 'C', 'A')
kras.2.percent.reads.r.melt <- melt(kras.2.percent.reads.r, id.vars = colnames(kras.2.percent.reads.r), value.name = 'Percent')
colnames(kras.2.percent.reads.r.melt) <- c('Nucleotide', 'cDNA_Position', 'Percent')

# Re-order the columns.  This makes the mutation report make more sense.
kras.2.percent.reads.r.melt <- kras.2.percent.reads.r.melt[c('cDNA_Position', 'Nucleotide', 'Percent')]


# First, assign some colours to a palette.  See 'http://www.hypergurl.com/rgbcolorchart.html' for hexadecimal codes.
seq.data.colours <- c('#FF0000', '#FF9900', '#0000FF', '#33CC00')

# STACKED BARCHARTS: FORWARD
# Use 'melted' data to make the barcharts. Start by adding a 'c.' to the cDNA number.  Each cDNA number will appear in the column 4 times, making the number of factors in the column = (number of rows)/4.  This makes the data suitable for a stacked barchart.
kras.2.percent.reads.f.melt.1 <- as.data.frame(paste('c.', kras.2.percent.reads.f.melt$cDNA_Position, sep = ''))
kras.2.percent.reads.f.melt.2 <- cbind(kras.2.percent.reads.f.melt.1, kras.2.percent.reads.f.melt)
colnames(kras.2.percent.reads.f.melt.2) <- c('cDNA_Position', 'Position', 'Nucleotide', 'Percent')

# Make a stacked barchart.  The plot has been assigned to an object to prevent it being printed.
kras.2.fwd.plot <- ggplot(kras.2.percent.reads.f.melt.2, aes(x = cDNA_Position, y = Percent, fill = Nucleotide)) + geom_bar(stat = 'identity') + scale_fill_manual(values = seq.data.colours) + labs(title = paste('Sample number', bam_filename, '\n\nStacked bar chart for forward reads in KRAS exon 2', sep = ' '), x = 'cDNA position', y = 'Percentage allele') + theme(axis.text.x = element_text(angle = 90))
ggsave(paste(path_to_report, KRAS_2_chart_f, sep = '/'))

# STACKED BARCHARTS: REVERSE
# Use 'melted' data to make the barcharts. Start by adding a 'c.' to the cDNA number.  Each cDNA number will appear in the column 4 times, making the number of factors in the column = (number of rows)/4.  This makes the data suitable for a stacked barchart.
kras.2.percent.reads.r.melt.1 <- as.data.frame(paste('c.', kras.2.percent.reads.r.melt$cDNA_Position, sep = ''))
kras.2.percent.reads.r.melt.2 <- cbind(kras.2.percent.reads.r.melt.1, kras.2.percent.reads.r.melt)
colnames(kras.2.percent.reads.r.melt.2) <- c('cDNA_Position', 'Position', 'Nucleotide', 'Percent')

# Make a stacked barchart.
kras.2.rev.plot <- ggplot(kras.2.percent.reads.r.melt.2, aes(x = cDNA_Position, y = Percent, fill = Nucleotide)) + geom_bar(stat = 'identity') + scale_fill_manual(values = seq.data.colours) + labs(title = paste('Sample number', bam_filename, '\n\nStacked bar chart for reverse reads in KRAS exon 2', sep = ' '), x = 'cDNA position', y = 'Percentage allele') + theme(axis.text.x = element_text(angle = 90))
ggsave(paste(path_to_report, KRAS_2_chart_r, sep = '/'))


# Mutation report: FORWARD 
# Add a 'c.' and the nucleotide to the melted data.  This means that the number of rows in the column = number of factors.  This makes the data suitable for a properly annotated mutation report. 
kras.2.percent.reads.f.melt.3 <- as.data.frame(paste('c.', kras.2.percent.reads.f.melt$cDNA_Position, kras.2.percent.reads.f.melt$Nucleotide, sep = ''))
kras.2.percent.reads.f.melt.4 <- cbind(kras.2.percent.reads.f.melt.3, kras.2.percent.reads.f.melt)
colnames(kras.2.percent.reads.f.melt.4) <- c('cDNA_Position', 'Position', 'Nucleotide', 'Percent')

# Filter out the codons of interest plus one either side.  Filter by the 'Position' column.
kras.2.percent.reads.f.melt.11_14 <- filter(kras.2.percent.reads.f.melt.4, Position >= 31, Position <= 42)

# Mutation report: REVERSE
# Add a 'c.' and the nucleotide to the melted data.  This means that the number of rows in the column = number of factors.  This makes the data suitable for a properly annotated mutation report.  
kras.2.percent.reads.r.melt.3 <- as.data.frame(paste('c.', kras.2.percent.reads.r.melt$cDNA_Position, kras.2.percent.reads.r.melt$Nucleotide, sep = ''))
kras.2.percent.reads.r.melt.4 <- cbind(kras.2.percent.reads.r.melt.3, kras.2.percent.reads.r.melt)
colnames(kras.2.percent.reads.r.melt.4) <- c('cDNA_Position', 'Position', 'Nucleotide', 'Percent')

# Filter out the codons of interest plus one either side.  Filter by the 'Position' column.
kras.2.percent.reads.r.melt.11_14 <- filter(kras.2.percent.reads.r.melt.4, Position >= 31, Position <= 42)

# Return positions that have nucleotides >=5% and <=95%.  This should detect nucleotide changes greater than 5%.
kras.2.muts_5_f <- filter(kras.2.percent.reads.f.melt, Percent >= 5.00, Percent <= 95.00); kras.2.muts_5_f
kras.2.muts_5_r <- filter(kras.2.percent.reads.r.melt, Percent >= 5.00, Percent <= 95.00); kras.2.muts_5_r

# Export the data to a mutation report.
write.table(paste('Sample number', bam_samplename, '\n\n', sep = ' '), paste(path_to_report, KRAS_2_mutations, sep = '/'), row.names = F, col.names = F, quote = F, append = T)
write.table(kras.2.muts_5_f, paste(path_to_report, KRAS_2_mutations, sep = '/'), row.names = F, quote = F, append = T)
write.table('\n', paste(path_to_report, KRAS_2_mutations, sep = '/'), row.names = F, col.names = F, quote = F, append = T)
write.table(kras.2.muts_5_r, paste(path_to_report, KRAS_2_mutations, sep = '/'), row.names = F, quote = F, append = T)


#=========================================================================================================
# KRAS EXON 3/CODONS 59+61

# Make two dataframes one for the '+' strand and one for the '-' strand.
kras.3.pup.f <- filter(kras.3.pup, strand == '+')
kras.3.pup.r <- filter(kras.3.pup, strand == '-')

# If there are >=2 nucleotides at one position the pileup function generates >=2 rows for that position.  This is particularly useless for mutation screening where the whole aim is to detect positions where there is more than 1 nucleotide...
# spread (tidyr) generates a dataframe that looks very like a consensus matrix.  spread (tidyr) needs to be told the dataframe it's working on, a 'key', ie name of the column whose values will be used as column headings and a 'value', ie name of the column whose values will populate the cells. 'fill' fills in the missing values.  Most importantly, it combines data for one position in to 1 column.
kras.3.pup.f.s <- spread(kras.3.pup.f, pos, count, fill = 0)
kras.3.pup.r.s <- spread(kras.3.pup.r, pos, count, fill = 0)

# Sum the counts in columns of the above dataframes.
kras.3.pup.f.s.sum <- apply(kras.3.pup.f.s[,4:182], 2, sum)
kras.3.pup.r.s.sum <- apply(kras.3.pup.r.s[,4:182], 2, sum)

# rbind results of apply operation to 'consensus matrix'.  This will lose the seqnames, strand and nucleotide columns.
kras.3.pup.2.f.s.2 <- rbind(kras.3.pup.f.s[,4:182], kras.3.pup.f.s.sum)
kras.3.pup.2.r.s.2 <- rbind(kras.3.pup.r.s[,4:182], kras.3.pup.r.s.sum)

# Calculate the proportion of counts for each nucleotide. 
kras.3.a_pc.reads.f <- kras.3.pup.2.f.s.2[1,]/kras.3.pup.2.f.s.2[nrow(kras.3.pup.2.f.s.2),]
kras.3.c_pc.reads.f <- kras.3.pup.2.f.s.2[2,]/kras.3.pup.2.f.s.2[nrow(kras.3.pup.2.f.s.2),]
kras.3.g_pc.reads.f <- kras.3.pup.2.f.s.2[3,]/kras.3.pup.2.f.s.2[nrow(kras.3.pup.2.f.s.2),]
kras.3.t_pc.reads.f <- kras.3.pup.2.f.s.2[4,]/kras.3.pup.2.f.s.2[nrow(kras.3.pup.2.f.s.2),]
kras.3.a_pc.reads.r <- kras.3.pup.2.r.s.2[1,]/kras.3.pup.2.r.s.2[nrow(kras.3.pup.2.r.s.2),]
kras.3.c_pc.reads.r <- kras.3.pup.2.r.s.2[2,]/kras.3.pup.2.r.s.2[nrow(kras.3.pup.2.r.s.2),]
kras.3.g_pc.reads.r <- kras.3.pup.2.r.s.2[3,]/kras.3.pup.2.r.s.2[nrow(kras.3.pup.2.r.s.2),]
kras.3.t_pc.reads.r <- kras.3.pup.2.r.s.2[4,]/kras.3.pup.2.r.s.2[nrow(kras.3.pup.2.r.s.2),]


# Make a single dataframe consisting of proportions for all nucleotides and convert in to percentages.
kras.3.percent.reads.f <- rbind(kras.3.a_pc.reads.f, kras.3.c_pc.reads.f, kras.3.g_pc.reads.f, kras.3.t_pc.reads.f)
kras.3.percent.reads.r <- rbind(kras.3.a_pc.reads.r, kras.3.c_pc.reads.r, kras.3.g_pc.reads.r, kras.3.t_pc.reads.r)

# The column names of 'kras.3.percent.reads.f' and 'kras.3.percent.reads.r' are gDNA co-ordinates.  This needs to be changed to cDNA positions.  This is done by making a dataframe of cDNA positions and assigning that as new column names.
# Make a dataframe of cDNA positions for KRAS exon 3. 
kras.3.cDNA <- as.data.frame(290:112)
kras.3.cDNA <- t(kras.3.cDNA)

# Assign the the above column names.
colnames(kras.3.percent.reads.f) <- kras.3.cDNA
colnames(kras.3.percent.reads.r) <- kras.3.cDNA

# Use 'apply' to convert proportions for each nucleotide in to percentages.
kras.3.percent.reads.f <- apply(kras.3.percent.reads.f, 1:2, function(x) round(as.numeric(x*100),2))
kras.3.percent.reads.r <- apply(kras.3.percent.reads.r, 1:2, function(x) round(as.numeric(x*100),2))

# Make a dataframe of 'seqnames', 'strand' and 'nucleotide' and attach to the percent reads dataframe.
# Forward
kras.3.percent.reads.f.info <- kras.3.pup.f.s[1:4,1:3]
kras.3.percent.reads.f.export <- cbind(kras.3.percent.reads.f.info, kras.3.percent.reads.f)
# As KRAS is in the reverse orientation the row names have been reversed too.  However, the 'Nucleotide' column is not reversed
row.names(kras.3.percent.reads.f.export) <- c('T', 'G', 'C', 'A')
# Reverse
kras.3.percent.reads.r.info <- kras.3.pup.r.s[1:4,1:3]
kras.3.percent.reads.r.export <- cbind(kras.3.percent.reads.r.info, kras.3.percent.reads.r)
# As KRAS is in the reverse orientation the row names have been reversed too.  However, the 'Nucleotide' column is not reversed
row.names(kras.3.percent.reads.r.export) <- c('T', 'G', 'C', 'A')


# Output percentage reads to a .csv file.
write.csv(kras.3.percent.reads.f.export, paste(path_to_report, KRAS_3_pc_reads_f, sep = '/') quote = F)
write.csv(kras.3.percent.reads.r.export, paste(path_to_report, KRAS_3_pc_reads_r, sep = '/') quote = F)

# This section rearranges the data into a form suitable for making a stacked bar chart with ggplot2

# melt (reshape2) converts 'wide' data into 'long' data, which is needed for making a stacked bar chart.
# Forward
rownames(kras.3.percent.reads.f) <- c('T', 'G', 'C', 'A')
kras.3.percent.reads.f.melt <- melt(kras.3.percent.reads.f, id.vars = colnames(kras.3.percent.reads.f), value.name = 'Percent')
colnames(kras.3.percent.reads.f.melt) <- c('Nucleotide', 'cDNA_Position', 'Percent')

# Re-order the columns.  This makes the mutation report make more sense.
kras.3.percent.reads.f.melt <- kras.3.percent.reads.f.melt[c('cDNA_Position', 'Nucleotide', 'Percent')]

# Reverse
rownames(kras.3.percent.reads.r) <- c('T', 'G', 'C', 'A')
kras.3.percent.reads.r.melt <- melt(kras.3.percent.reads.r, id.vars = colnames(kras.3.percent.reads.r), value.name = 'Percent')
colnames(kras.3.percent.reads.r.melt) <- c('Nucleotide', 'cDNA_Position', 'Percent')

# Re-order the columns.  This makes the mutation report make more sense.
kras.3.percent.reads.r.melt <- kras.3.percent.reads.r.melt[c('cDNA_Position', 'Nucleotide', 'Percent')]


# First, assign some colours to a palette.  See 'http://www.hypergurl.com/rgbcolorchart.html' for hexadecimal codes.
seq.data.colours <- c('#FF0000', '#FF9900', '#0000FF', '#33CC00')

# STACKED BARCHARTS: FORWARD
# Use 'melted' data to make the barcharts. Start by adding a 'c.' to the cDNA number.  Each cDNA number will appear in the column 4 times, making the number of factors in the column = (number of rows)/4.  This makes the data suitable for a stacked barchart.
kras.3.percent.reads.f.melt.1 <- as.data.frame(paste('c.', kras.3.percent.reads.f.melt$cDNA_Position, sep = ''))
kras.3.percent.reads.f.melt.2 <- cbind(kras.3.percent.reads.f.melt.1, kras.3.percent.reads.f.melt)
colnames(kras.3.percent.reads.f.melt.2) <- c('cDNA_Position', 'Position', 'Nucleotide', 'Percent')

# Make a stacked barchart.  The plot has been assigned to an object to prevent it being printed.
kras.3.fwd.plot <- ggplot(kras.3.percent.reads.f.melt.2, aes(x = cDNA_Position, y = Percent, fill = Nucleotide)) + geom_bar(stat = 'identity') + scale_fill_manual(values = seq.data.colours) + labs(title = paste('Sample number', bam_filename, '\n\nStacked bar chart for forward reads in KRAS exon 3', sep = ' '), x = 'cDNA position', y = 'Percentage allele') + theme(axis.text.x = element_text(angle = 90))
ggsave(paste(path_to_report, KRAS_3_chart_f, sep = '/'))

# STACKED BARCHARTS: REVERSE
# Use 'melted' data to make the barcharts. Start by adding a 'c.' to the cDNA number.  Each cDNA number will appear in the column 4 times, making the number of factors in the column = (number of rows)/4.  This makes the data suitable for a stacked barchart.
kras.3.percent.reads.r.melt.1 <- as.data.frame(paste('c.', kras.3.percent.reads.r.melt$cDNA_Position, sep = ''))
kras.3.percent.reads.r.melt.2 <- cbind(kras.3.percent.reads.r.melt.1, kras.3.percent.reads.r.melt)
colnames(kras.3.percent.reads.r.melt.2) <- c('cDNA_Position', 'Position', 'Nucleotide', 'Percent')

# Make a stacked barchart.
kras.3.rev.plot <- ggplot(kras.3.percent.reads.r.melt.2, aes(x = cDNA_Position, y = Percent, fill = Nucleotide)) + geom_bar(stat = 'identity') + scale_fill_manual(values = seq.data.colours) + labs(title = paste('Sample number', bam_filename, '\n\nStacked bar chart for reverse reads in KRAS exon 3', sep = ' '), x = 'cDNA position', y = 'Percentage allele') + theme(axis.text.x = element_text(angle = 90))
ggsave(paste(path_to_report, KRAS_3_chart_r, sep = '/'))


# Mutation report: FORWARD 
# Add a 'c.' and the nucleotide to the melted data.  This means that the number of rows in the column = number of factors.  This makes the data suitable for a properly annotated mutation report. 
kras.3.percent.reads.f.melt.3 <- as.data.frame(paste('c.', kras.3.percent.reads.f.melt$cDNA_Position, kras.3.percent.reads.f.melt$Nucleotide, sep = ''))
kras.3.percent.reads.f.melt.4 <- cbind(kras.3.percent.reads.f.melt.3, kras.3.percent.reads.f.melt)
colnames(kras.3.percent.reads.f.melt.4) <- c('cDNA_Position', 'Position', 'Nucleotide', 'Percent')

# Filter out the codons of interest plus one either side.  Filter by the 'Position' column.
kras.3.percent.reads.f.melt.58_62 <- filter(kras.3.percent.reads.f.melt.4, Position >= 172, Position <= 186)

# Mutation report: REVERSE
# Add a 'c.' and the nucleotide to the melted data.  This means that the number of rows in the column = number of factors.  This makes the data suitable for a properly annotated mutation report.  
kras.3.percent.reads.r.melt.3 <- as.data.frame(paste('c.', kras.3.percent.reads.r.melt$cDNA_Position, kras.3.percent.reads.r.melt$Nucleotide, sep = ''))
kras.3.percent.reads.r.melt.4 <- cbind(kras.3.percent.reads.r.melt.3, kras.3.percent.reads.r.melt)
colnames(kras.3.percent.reads.r.melt.4) <- c('cDNA_Position', 'Position', 'Nucleotide', 'Percent')

# Filter out the codons of interest plus one either side.  Filter by the 'Position' column.
kras.3.percent.reads.r.melt.58_62 <- filter(kras.3.percent.reads.r.melt.4, Position >= 172, Position <= 186)

# Return positions that have nucleotides >=5% and <=95%.  This should detect nucleotide changes greater than 5%.
kras.3.muts_5_f <- filter(kras.3.percent.reads.f.melt, Percent >= 5.00, Percent <= 95.00); kras.3.muts_5_f
kras.3.muts_5_r <- filter(kras.3.percent.reads.r.melt, Percent >= 5.00, Percent <= 95.00); kras.3.muts_5_r

# Export the data to a mutation report.
write.table(paste('Sample number', bam_samplename, '\n\n', sep = ' '), paste(path_to_report, KRAS_3_mutations, sep = '/'), row.names = F, col.names = F, quote = F, append = T)
write.table(kras.3.muts_5_f, paste(path_to_report, KRAS_3_mutations, sep = '/'), row.names = F, quote = F, append = T)
write.table('\n', paste(path_to_report, KRAS_3_mutations, sep = '/'), row.names = F, col.names = F, quote = F, append = T)
write.table(kras.3.muts_5_r, paste(path_to_report, KRAS_3_mutations, sep = '/'), row.names = F, quote = F, append = T)


#=========================================================================================================
# KRAS EXON 4/CODONS 117+146

# Make two dataframes one for the '+' strand and one for the '-' strand.
kras.4.pup.f <- filter(kras.4.pup, strand == '+')
kras.4.pup.r <- filter(kras.4.pup, strand == '-')

# If there are >=2 nucleotides at one position the pileup function generates >=2 rows for that position.  This is particularly useless for mutation screening where the whole aim is to detect positions where there is more than 1 nucleotide...
# spread (tidyr) generates a dataframe that looks very like a consensus matrix.  spread (tidyr) needs to be told the dataframe it's working on, a 'key', ie name of the column whose values will be used as column headings and a 'value', ie name of the column whose values will populate the cells. 'fill' fills in the missing values.  Most importantly, it combines data for one position in to 1 column.
kras.4.pup.f.s <- spread(kras.4.pup.f, pos, count, fill = 0)
kras.4.pup.r.s <- spread(kras.4.pup.r, pos, count, fill = 0)

# Sum the counts in columns of the above dataframes.
kras.4.pup.f.s.sum <- apply(kras.4.pup.f.s[,4:163], 2, sum)
kras.4.pup.r.s.sum <- apply(kras.4.pup.r.s[,4:163], 2, sum)

# rbind results of apply operation to 'consensus matrix'.  This will lose the seqnames, strand and nucleotide columns.
kras.4.pup.2.f.s.2 <- rbind(kras.4.pup.f.s[,4:163], kras.4.pup.f.s.sum)
kras.4.pup.2.r.s.2 <- rbind(kras.4.pup.r.s[,4:163], kras.4.pup.r.s.sum)

# Calculate the proportion of counts for each nucleotide. 
kras.4.a_pc.reads.f <- kras.4.pup.2.f.s.2[1,]/kras.4.pup.2.f.s.2[nrow(kras.4.pup.2.f.s.2),]
kras.4.c_pc.reads.f <- kras.4.pup.2.f.s.2[2,]/kras.4.pup.2.f.s.2[nrow(kras.4.pup.2.f.s.2),]
kras.4.g_pc.reads.f <- kras.4.pup.2.f.s.2[3,]/kras.4.pup.2.f.s.2[nrow(kras.4.pup.2.f.s.2),]
kras.4.t_pc.reads.f <- kras.4.pup.2.f.s.2[4,]/kras.4.pup.2.f.s.2[nrow(kras.4.pup.2.f.s.2),]
kras.4.a_pc.reads.r <- kras.4.pup.2.r.s.2[1,]/kras.4.pup.2.r.s.2[nrow(kras.4.pup.2.r.s.2),]
kras.4.c_pc.reads.r <- kras.4.pup.2.r.s.2[2,]/kras.4.pup.2.r.s.2[nrow(kras.4.pup.2.r.s.2),]
kras.4.g_pc.reads.r <- kras.4.pup.2.r.s.2[3,]/kras.4.pup.2.r.s.2[nrow(kras.4.pup.2.r.s.2),]
kras.4.t_pc.reads.r <- kras.4.pup.2.r.s.2[4,]/kras.4.pup.2.r.s.2[nrow(kras.4.pup.2.r.s.2),]


# Make a single dataframe consisting of proportions for all nucleotides and convert in to percentages.
kras.4.percent.reads.f <- rbind(kras.4.a_pc.reads.f, kras.4.c_pc.reads.f, kras.4.g_pc.reads.f, kras.4.t_pc.reads.f)
kras.4.percent.reads.r <- rbind(kras.4.a_pc.reads.r, kras.4.c_pc.reads.r, kras.4.g_pc.reads.r, kras.4.t_pc.reads.r)

# The column names of 'kras.4.percent.reads.f' and 'kras.4.percent.reads.r' are gDNA co-ordinates.  This needs to be changed to cDNA positions.  This is done by making a dataframe of cDNA positions and assigning that as new column names.
# Make a dataframe of cDNA positions for KRAS exon 4. 
kras.4.cDNA <- as.data.frame(450:291)
kras.4.cDNA <- t(kras.4.cDNA)

# Assign the the above column names.
colnames(kras.4.percent.reads.f) <- kras.4.cDNA
colnames(kras.4.percent.reads.r) <- kras.4.cDNA

# Use 'apply' to convert proportions for each nucleotide in to percentages.
kras.4.percent.reads.f <- apply(kras.4.percent.reads.f, 1:2, function(x) round(as.numeric(x*100),2))
kras.4.percent.reads.r <- apply(kras.4.percent.reads.r, 1:2, function(x) round(as.numeric(x*100),2))

# Make a dataframe of 'seqnames', 'strand' and 'nucleotide' and attach to the percent reads dataframe.
# Forward
kras.4.percent.reads.f.info <- kras.4.pup.f.s[1:4,1:3]
kras.4.percent.reads.f.export <- cbind(kras.4.percent.reads.f.info, kras.4.percent.reads.f)
# As KRAS is in the reverse orientation the row names have been reversed too.  However, the 'Nucleotide' column is not reversed
row.names(kras.4.percent.reads.f.export) <- c('T', 'G', 'C', 'A')
# Reverse
kras.4.percent.reads.r.info <- kras.4.pup.r.s[1:4,1:3]
kras.4.percent.reads.r.export <- cbind(kras.4.percent.reads.r.info, kras.4.percent.reads.r)
# As KRAS is in the reverse orientation the row names have been reversed too.  However, the 'Nucleotide' column is not reversed
row.names(kras.4.percent.reads.r.export) <- c('T', 'G', 'C', 'A')


# Output percentage reads to a .csv file.
write.csv(kras.4.percent.reads.f.export, paste(path_to_report, KRAS_4_pc_reads_f, sep = '/') quote = F)
write.csv(kras.4.percent.reads.r.export, paste(path_to_report, KRAS_4_pc_reads_r, sep = '/') quote = F)

# This section rearranges the data into a form suitable for making a stacked bar chart with ggplot2

# melt (reshape2) converts 'wide' data into 'long' data, which is needed for making a stacked bar chart.
# Forward
rownames(kras.4.percent.reads.f) <- c('T', 'G', 'C', 'A')
kras.4.percent.reads.f.melt <- melt(kras.4.percent.reads.f, id.vars = colnames(kras.4.percent.reads.f), value.name = 'Percent')
colnames(kras.4.percent.reads.f.melt) <- c('Nucleotide', 'cDNA_Position', 'Percent')

# Re-order the columns.  This makes the mutation report make more sense.
kras.4.percent.reads.f.melt <- kras.4.percent.reads.f.melt[c('cDNA_Position', 'Nucleotide', 'Percent')]

# Reverse
rownames(kras.4.percent.reads.r) <- c('T', 'G', 'C', 'A')
kras.4.percent.reads.r.melt <- melt(kras.4.percent.reads.r, id.vars = colnames(kras.4.percent.reads.r), value.name = 'Percent')
colnames(kras.4.percent.reads.r.melt) <- c('Nucleotide', 'cDNA_Position', 'Percent')

# Re-order the columns.  This makes the mutation report make more sense.
kras.4.percent.reads.r.melt <- kras.4.percent.reads.r.melt[c('cDNA_Position', 'Nucleotide', 'Percent')]


# First, assign some colours to a palette.  See 'http://www.hypergurl.com/rgbcolorchart.html' for hexadecimal codes.
seq.data.colours <- c('#FF0000', '#FF9900', '#0000FF', '#33CC00')

# STACKED BARCHARTS: FORWARD
# Use 'melted' data to make the barcharts. Start by adding a 'c.' to the cDNA number.  Each cDNA number will appear in the column 4 times, making the number of factors in the column = (number of rows)/4.  This makes the data suitable for a stacked barchart.
kras.4.percent.reads.f.melt.1 <- as.data.frame(paste('c.', kras.4.percent.reads.f.melt$cDNA_Position, sep = ''))
kras.4.percent.reads.f.melt.2 <- cbind(kras.4.percent.reads.f.melt.1, kras.4.percent.reads.f.melt)
colnames(kras.4.percent.reads.f.melt.2) <- c('cDNA_Position', 'Position', 'Nucleotide', 'Percent')

# Make a stacked barchart.  The plot has been assigned to an object to prevent it being printed.
kras.4.fwd.plot <- ggplot(kras.4.percent.reads.f.melt.2, aes(x = cDNA_Position, y = Percent, fill = Nucleotide)) + geom_bar(stat = 'identity') + scale_fill_manual(values = seq.data.colours) + labs(title = paste('Sample number', bam_filename, '\n\nStacked bar chart for forward reads in KRAS exon 4', sep = ' '), x = 'cDNA position', y = 'Percentage allele') + theme(axis.text.x = element_text(angle = 90))
ggsave(paste(path_to_report, KRAS_4_chart_f, sep = '/'))

# STACKED BARCHARTS: REVERSE
# Use 'melted' data to make the barcharts. Start by adding a 'c.' to the cDNA number.  Each cDNA number will appear in the column 4 times, making the number of factors in the column = (number of rows)/4.  This makes the data suitable for a stacked barchart.
kras.4.percent.reads.r.melt.1 <- as.data.frame(paste('c.', kras.4.percent.reads.r.melt$cDNA_Position, sep = ''))
kras.4.percent.reads.r.melt.2 <- cbind(kras.4.percent.reads.r.melt.1, kras.4.percent.reads.r.melt)
colnames(kras.4.percent.reads.r.melt.2) <- c('cDNA_Position', 'Position', 'Nucleotide', 'Percent')

# Make a stacked barchart.
kras.4.rev.plot <- ggplot(kras.4.percent.reads.r.melt.2, aes(x = cDNA_Position, y = Percent, fill = Nucleotide)) + geom_bar(stat = 'identity') + scale_fill_manual(values = seq.data.colours) + labs(title = paste('Sample number', bam_filename, '\n\nStacked bar chart for reverse reads in KRAS exon 4', sep = ' '), x = 'cDNA position', y = 'Percentage allele') + theme(axis.text.x = element_text(angle = 90))
ggsave(paste(path_to_report, KRAS_4_chart_r, sep = '/'))


# Mutation report: FORWARD 
# Add a 'c.' and the nucleotide to the melted data.  This means that the number of rows in the column = number of factors.  This makes the data suitable for a properly annotated mutation report. 
kras.4.percent.reads.f.melt.3 <- as.data.frame(paste('c.', kras.4.percent.reads.f.melt$cDNA_Position, kras.4.percent.reads.f.melt$Nucleotide, sep = ''))
kras.4.percent.reads.f.melt.4 <- cbind(kras.4.percent.reads.f.melt.3, kras.4.percent.reads.f.melt)
colnames(kras.4.percent.reads.f.melt.4) <- c('cDNA_Position', 'Position', 'Nucleotide', 'Percent')

# Filter out the codons of interest plus one either side.  Filter by the 'Position' column.
kras.4.percent.reads.f.melt.116_147 <- filter(kras.4.percent.reads.f.melt.4, Position >= 346, Position <= 441)

# Mutation report: REVERSE
# Add a 'c.' and the nucleotide to the melted data.  This means that the number of rows in the column = number of factors.  This makes the data suitable for a properly annotated mutation report.  
kras.4.percent.reads.r.melt.3 <- as.data.frame(paste('c.', kras.4.percent.reads.r.melt$cDNA_Position, kras.4.percent.reads.r.melt$Nucleotide, sep = ''))
kras.4.percent.reads.r.melt.4 <- cbind(kras.4.percent.reads.r.melt.3, kras.4.percent.reads.r.melt)
colnames(kras.4.percent.reads.r.melt.4) <- c('cDNA_Position', 'Position', 'Nucleotide', 'Percent')

# Filter out the codons of interest plus one either side.  Filter by the 'Position' column.
kras.4.percent.reads.r.melt.116_147 <- filter(kras.4.percent.reads.r.melt.4, Position >= 346, Position <= 441)

# Return positions that have nucleotides >=5% and <=95%.  This should detect nucleotide changes greater than 5%.
kras.4.muts_5_f <- filter(kras.4.percent.reads.f.melt, Percent >= 5.00, Percent <= 95.00); kras.4.muts_5_f
kras.4.muts_5_r <- filter(kras.4.percent.reads.r.melt, Percent >= 5.00, Percent <= 95.00); kras.4.muts_5_r

# Export the data to a mutation report.
write.table(paste('Sample number', bam_samplename, '\n\n', sep = ' '), paste(path_to_report, KRAS_4_mutations, sep = '/'), row.names = F, col.names = F, quote = F, append = T)
write.table(kras.4.muts_5_f, paste(path_to_report, kras_4_mutations, sep = '/'), row.names = F, quote = F, append = T)
write.table('\n', paste(path_to_report, KRAS_4_mutations, sep = '/'), row.names = F, col.names = F, quote = F, append = T)
write.table(kras.4.muts_5_r, paste(path_to_report, KRAS_4_mutations, sep = '/'), row.names = F, quote = F, append = T)

