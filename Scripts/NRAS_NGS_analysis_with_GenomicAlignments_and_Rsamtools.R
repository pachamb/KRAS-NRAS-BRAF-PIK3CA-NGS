# Author: Phil Chambers
# Date of last edit: 24/06/2019

# Purpose: 
# This analyses short read, next generation sequencing data for the human NRAS gene generated from an Illumina Genome Analyser using the Qiagen "GeneRead DNAseq Targeted Panels V2 Human Tumor Actionable Mutations Panel" (Product no. 181900, Cat. no. NGHS-201X).  It should also work with data generated from an Illumina Genome Analyser using and any PCR primers that amplify the regions of interest. Gene alignment and numbering is for human genome build hg19.
# It generates three outputs:
# 1. csv files with the percentage reads for each nucleotide at each position.
# 2. pdf files of stacked bar charts showing the percentage of each nucleotide at each position.
# 3. A txt file containing any mutations detected.


# NRAS EXONS 2 AND 3

# bam file scanning parameters.  These are for NRAS exons 2, 3 and 4.
nras.all.param <- ScanBamParam(which = GRanges('chr12', IRanges(start = c(115258671, 115256421), end = c(115258781, 115256599))))

# Pileup parameters.
nras.all.pup.param <- PileupParam(max_depth = 10000, min_base_quality = 20, min_nucleotide_depth = 10, distinguish_strands = T)

# Perform the pileup operation
nras.all.pup <- pileup(sample_files, scanBamParam = nras.all.param, pileupParam = nras.all.pup.param)

# Remove column called 'which_label' and keep the rest
nras.all.pup.2 <- select(nras.all.pup, -which_label)

# Split the data into 3 dataframes for the 3 exons.
nras.2.pup <- filter(nras.all.pup.2, pos >= 115258671, pos <= 115258781)
nras.3.pup <- filter(nras.all.pup.2, pos >= 115256421, pos <= 115256599)

# NRAS EXON 2/CODONS 12+13

# Make two dataframes one for the '+' strand and one for the '-' strand.
nras.2.pup.f <- filter(nras.2.pup, strand == '+')
nras.2.pup.r <- filter(nras.2.pup, strand == '-')

# If there are >=2 nucleotides at one position the pileup function generates >=2 rows for that position.  This is particularly useless for mutation screening where the whole aim is to detect positions where there is more than 1 nucleotide...
# spread (tidyr) generates a dataframe that looks very like a consensus matrix.  spread (tidyr) needs to be told the dataframe it's working on, a 'key', ie name of the column whose values will be used as column headings and a 'value', ie name of the column whose values will populate the cells. 'fill' fills in the missing values.  Most importantly, it combines data for one position in to 1 column.
nras.2.pup.f.s <- spread(nras.2.pup.f, pos, count, fill = 0)
nras.2.pup.r.s <- spread(nras.2.pup.r, pos, count, fill = 0)

# Sum the counts in columns of the above dataframes.
nras.2.pup.f.s.sum <- apply(nras.2.pup.f.s[,4:114], 2, sum)
nras.2.pup.r.s.sum <- apply(nras.2.pup.r.s[,4:114], 2, sum)

# rbind results of apply operation to 'consensus matrix'.  This will lose the seqnames, strand and nucleotide columns.
nras.2.pup.2.f.s.2 <- rbind(nras.2.pup.f.s[,4:114], nras.2.pup.f.s.sum)
nras.2.pup.2.r.s.2 <- rbind(nras.2.pup.r.s[,4:114], nras.2.pup.r.s.sum)

# Calculate the proportion of counts for each nucleotide. 
nras.2.a_pc.reads.f <- nras.2.pup.2.f.s.2[1,]/nras.2.pup.2.f.s.2[nrow(nras.2.pup.2.f.s.2),]
nras.2.c_pc.reads.f <- nras.2.pup.2.f.s.2[2,]/nras.2.pup.2.f.s.2[nrow(nras.2.pup.2.f.s.2),]
nras.2.g_pc.reads.f <- nras.2.pup.2.f.s.2[3,]/nras.2.pup.2.f.s.2[nrow(nras.2.pup.2.f.s.2),]
nras.2.t_pc.reads.f <- nras.2.pup.2.f.s.2[4,]/nras.2.pup.2.f.s.2[nrow(nras.2.pup.2.f.s.2),]
nras.2.a_pc.reads.r <- nras.2.pup.2.r.s.2[1,]/nras.2.pup.2.r.s.2[nrow(nras.2.pup.2.r.s.2),]
nras.2.c_pc.reads.r <- nras.2.pup.2.r.s.2[2,]/nras.2.pup.2.r.s.2[nrow(nras.2.pup.2.r.s.2),]
nras.2.g_pc.reads.r <- nras.2.pup.2.r.s.2[3,]/nras.2.pup.2.r.s.2[nrow(nras.2.pup.2.r.s.2),]
nras.2.t_pc.reads.r <- nras.2.pup.2.r.s.2[4,]/nras.2.pup.2.r.s.2[nrow(nras.2.pup.2.r.s.2),]


# Make a single dataframe consisting of proportions for all nucleotides and convert in to percentages.
nras.2.percent.reads.f <- rbind(nras.2.a_pc.reads.f, nras.2.c_pc.reads.f, nras.2.g_pc.reads.f, nras.2.t_pc.reads.f)
nras.2.percent.reads.r <- rbind(nras.2.a_pc.reads.r, nras.2.c_pc.reads.r, nras.2.g_pc.reads.r, nras.2.t_pc.reads.r)

# The column names of 'nras.2.percent.reads.f' and 'nras.2.percent.reads.r' are gDNA co-ordinates.  This needs to be changed to cDNA positions.  This is done by making a dataframe of cDNA positions and assigning that as new column names.
# Make a dataframe of cDNA positions for NRAS exon 2. 
nras.2.cDNA <- as.data.frame(111:1)
nras.2.cDNA <- t(nras.2.cDNA)

# Assign the the above column names.
colnames(nras.2.percent.reads.f) <- nras.2.cDNA
colnames(nras.2.percent.reads.r) <- nras.2.cDNA

# Use 'apply' to convert proportions for each nucleotide in to percentages.
nras.2.percent.reads.f <- apply(nras.2.percent.reads.f, 1:2, function(x) round(as.numeric(x*100),2))
nras.2.percent.reads.r <- apply(nras.2.percent.reads.r, 1:2, function(x) round(as.numeric(x*100),2))

# Make a dataframe of 'seqnames', 'strand' and 'nucleotide' and attach to the percent reads dataframe.
# Forward
nras.2.percent.reads.f.info <- nras.2.pup.f.s[1:4,1:3]
nras.2.percent.reads.f.export <- cbind(nras.2.percent.reads.f.info, nras.2.percent.reads.f)
# As NRAS is in the reverse orientation the row names have been reversed too.  However, the 'Nucleotide' column is not reversed
row.names(nras.2.percent.reads.f.export) <- c('T', 'G', 'C', 'A')
# Reverse
nras.2.percent.reads.r.info <- nras.2.pup.r.s[1:4,1:3]
nras.2.percent.reads.r.export <- cbind(nras.2.percent.reads.r.info, nras.2.percent.reads.r)
# As NRAS is in the reverse orientation the row names have been reversed too.  However, the 'Nucleotide' column is not reversed
row.names(nras.2.percent.reads.r.export) <- c('T', 'G', 'C', 'A')


# Output percentage reads to a .csv file.
write.csv(nras.2.percent.reads.f.export, paste(path_to_report, NRAS_2_pc_reads_f, sep = '/') quote = F)
write.csv(nras.2.percent.reads.r.export, paste(path_to_report, NRAS_2_pc_reads_r, sep = '/') quote = F)

# This section rearranges the data into a form suitable for making a stacked bar chart with ggplot2

# melt (reshape2) converts 'wide' data into 'long' data, which is needed for making a stacked bar chart.
# Forward
rownames(nras.2.percent.reads.f) <- c('T', 'G', 'C', 'A')
nras.2.percent.reads.f.melt <- melt(nras.2.percent.reads.f, id.vars = colnames(nras.2.percent.reads.f), value.name = 'Percent')
colnames(nras.2.percent.reads.f.melt) <- c('Nucleotide', 'cDNA_Position', 'Percent')

# Re-order the columns.  This makes the mutation report make more sense.
nras.2.percent.reads.f.melt <- nras.2.percent.reads.f.melt[c('cDNA_Position', 'Nucleotide', 'Percent')]

# Reverse
rownames(nras.2.percent.reads.r) <- c('T', 'G', 'C', 'A')
nras.2.percent.reads.r.melt <- melt(nras.2.percent.reads.r, id.vars = colnames(nras.2.percent.reads.r), value.name = 'Percent')
colnames(nras.2.percent.reads.r.melt) <- c('Nucleotide', 'cDNA_Position', 'Percent')

# Re-order the columns.  This makes the mutation report make more sense.
nras.2.percent.reads.r.melt <- nras.2.percent.reads.r.melt[c('cDNA_Position', 'Nucleotide', 'Percent')]


# First, assign some colours to a palette.  See 'http://www.hypergurl.com/rgbcolorchart.html' for hexadecimal codes.
seq.data.colours <- c('#FF0000', '#FF9900', '#0000FF', '#33CC00')

# STACKED BARCHARTS: FORWARD
# Use 'melted' data to make the barcharts. Start by adding a 'c.' to the cDNA number.  Each cDNA number will appear in the column 4 times, making the number of factors in the column = (number of rows)/4.  This makes the data suitable for a stacked barchart.
nras.2.percent.reads.f.melt.1 <- as.data.frame(paste('c.', nras.2.percent.reads.f.melt$cDNA_Position, sep = ''))
nras.2.percent.reads.f.melt.2 <- cbind(nras.2.percent.reads.f.melt.1, nras.2.percent.reads.f.melt)
colnames(nras.2.percent.reads.f.melt.2) <- c('cDNA_Position', 'Position', 'Nucleotide', 'Percent')

# Make a stacked barchart.  The plot has been assigned to an object to prevent it being printed.
nras.2.fwd.plot <- ggplot(nras.2.percent.reads.f.melt.2, aes(x = cDNA_Position, y = Percent, fill = Nucleotide)) + geom_bar(stat = 'identity') + scale_fill_manual(values = seq.data.colours) + labs(title = paste('Sample number', bam_filename, '\n\nStacked bar chart for forward reads in NRAS exon 2', sep = ' '), x = 'cDNA position', y = 'Percentage allele') + theme(axis.text.x = element_text(angle = 90))
ggsave(paste(path_to_report, NRAS_2_chart_f, sep = '/'))

# STACKED BARCHARTS: REVERSE
# Use 'melted' data to make the barcharts. Start by adding a 'c.' to the cDNA number.  Each cDNA number will appear in the column 4 times, making the number of factors in the column = (number of rows)/4.  This makes the data suitable for a stacked barchart.
nras.2.percent.reads.r.melt.1 <- as.data.frame(paste('c.', nras.2.percent.reads.r.melt$cDNA_Position, sep = ''))
nras.2.percent.reads.r.melt.2 <- cbind(nras.2.percent.reads.r.melt.1, nras.2.percent.reads.r.melt)
colnames(nras.2.percent.reads.r.melt.2) <- c('cDNA_Position', 'Position', 'Nucleotide', 'Percent')

# Make a stacked barchart.
nras.2.rev.plot <- ggplot(nras.2.percent.reads.r.melt.2, aes(x = cDNA_Position, y = Percent, fill = Nucleotide)) + geom_bar(stat = 'identity') + scale_fill_manual(values = seq.data.colours) + labs(title = paste('Sample number', bam_filename, '\n\nStacked bar chart for reverse reads in NRAS exon 2', sep = ' '), x = 'cDNA position', y = 'Percentage allele') + theme(axis.text.x = element_text(angle = 90))
ggsave(paste(path_to_report, NRAS_2_chart_r, sep = '/'))


# Mutation report: FORWARD 
# Add a 'c.' and the nucleotide to the melted data.  This means that the number of rows in the column = number of factors.  This makes the data suitable for a properly annotated mutation report. 
nras.2.percent.reads.f.melt.3 <- as.data.frame(paste('c.', nras.2.percent.reads.f.melt$cDNA_Position, nras.2.percent.reads.f.melt$Nucleotide, sep = ''))
nras.2.percent.reads.f.melt.4 <- cbind(nras.2.percent.reads.f.melt.3, nras.2.percent.reads.f.melt)
colnames(nras.2.percent.reads.f.melt.4) <- c('cDNA_Position', 'Position', 'Nucleotide', 'Percent')

# Filter out the codons of interest plus one either side.  Filter by the 'Position' column.
nras.2.percent.reads.f.melt.11_14 <- filter(nras.2.percent.reads.f.melt.4, Position >= 31, Position <= 42)

# Mutation report: REVERSE
# Add a 'c.' and the nucleotide to the melted data.  This means that the number of rows in the column = number of factors.  This makes the data suitable for a properly annotated mutation report.  
nras.2.percent.reads.r.melt.3 <- as.data.frame(paste('c.', nras.2.percent.reads.r.melt$cDNA_Position, nras.2.percent.reads.r.melt$Nucleotide, sep = ''))
nras.2.percent.reads.r.melt.4 <- cbind(nras.2.percent.reads.r.melt.3, nras.2.percent.reads.r.melt)
colnames(nras.2.percent.reads.r.melt.4) <- c('cDNA_Position', 'Position', 'Nucleotide', 'Percent')

# Filter out the codons of interest plus one either side.  Filter by the 'Position' column.
nras.2.percent.reads.r.melt.11_14 <- filter(nras.2.percent.reads.r.melt.4, Position >= 31, Position <= 42)

# Return positions that have nucleotides >=5% and <=95%.  This should detect nucleotide changes greater than 5%.
nras.2.muts_5_f <- filter(nras.2.percent.reads.f.melt, Percent >= 5.00, Percent <= 95.00); nras.2.muts_5_f
nras.2.muts_5_r <- filter(nras.2.percent.reads.r.melt, Percent >= 5.00, Percent <= 95.00); nras.2.muts_5_r

# Export the data to a mutation report.
write.table(paste('Sample number', bam_samplename, '\n\n', sep = ' '), paste(path_to_report, NRAS_2_mutations, sep = '/'), row.names = F, col.names = F, quote = F, append = T)
write.table(nras.2.muts_5_f, paste(path_to_report, NRAS_2_mutations, sep = '/'), row.names = F, quote = F, append = T)
write.table('\n', paste(path_to_report, NRAS_2_mutations, sep = '/'), row.names = F, col.names = F, quote = F, append = T)
write.table(nras.2.muts_5_r, paste(path_to_report, NRAS_2_mutations, sep = '/'), row.names = F, quote = F, append = T)


#=========================================================================================================
# NRAS EXON 3/CODONS 59+61

# Make two dataframes one for the '+' strand and one for the '-' strand.
nras.3.pup.f <- filter(nras.3.pup, strand == '+')
nras.3.pup.r <- filter(nras.3.pup, strand == '-')

# If there are >=2 nucleotides at one position the pileup function generates >=2 rows for that position.  This is particularly useless for mutation screening where the whole aim is to detect positions where there is more than 1 nucleotide...
# spread (tidyr) generates a dataframe that looks very like a consensus matrix.  spread (tidyr) needs to be told the dataframe it's working on, a 'key', ie name of the column whose values will be used as column headings and a 'value', ie name of the column whose values will populate the cells. 'fill' fills in the missing values.  Most importantly, it combines data for one position in to 1 column.
nras.3.pup.f.s <- spread(nras.3.pup.f, pos, count, fill = 0)
nras.3.pup.r.s <- spread(nras.3.pup.r, pos, count, fill = 0)

# Sum the counts in columns of the above dataframes.
nras.3.pup.f.s.sum <- apply(nras.3.pup.f.s[,4:182], 2, sum)
nras.3.pup.r.s.sum <- apply(nras.3.pup.r.s[,4:182], 2, sum)

# rbind results of apply operation to 'consensus matrix'.  This will lose the seqnames, strand and nucleotide columns.
nras.3.pup.2.f.s.2 <- rbind(nras.3.pup.f.s[,4:182], nras.3.pup.f.s.sum)
nras.3.pup.2.r.s.2 <- rbind(nras.3.pup.r.s[,4:182], nras.3.pup.r.s.sum)

# Calculate the proportion of counts for each nucleotide. 
nras.3.a_pc.reads.f <- nras.3.pup.2.f.s.2[1,]/nras.3.pup.2.f.s.2[nrow(nras.3.pup.2.f.s.2),]
nras.3.c_pc.reads.f <- nras.3.pup.2.f.s.2[2,]/nras.3.pup.2.f.s.2[nrow(nras.3.pup.2.f.s.2),]
nras.3.g_pc.reads.f <- nras.3.pup.2.f.s.2[3,]/nras.3.pup.2.f.s.2[nrow(nras.3.pup.2.f.s.2),]
nras.3.t_pc.reads.f <- nras.3.pup.2.f.s.2[4,]/nras.3.pup.2.f.s.2[nrow(nras.3.pup.2.f.s.2),]
nras.3.a_pc.reads.r <- nras.3.pup.2.r.s.2[1,]/nras.3.pup.2.r.s.2[nrow(nras.3.pup.2.r.s.2),]
nras.3.c_pc.reads.r <- nras.3.pup.2.r.s.2[2,]/nras.3.pup.2.r.s.2[nrow(nras.3.pup.2.r.s.2),]
nras.3.g_pc.reads.r <- nras.3.pup.2.r.s.2[3,]/nras.3.pup.2.r.s.2[nrow(nras.3.pup.2.r.s.2),]
nras.3.t_pc.reads.r <- nras.3.pup.2.r.s.2[4,]/nras.3.pup.2.r.s.2[nrow(nras.3.pup.2.r.s.2),]


# Make a single dataframe consisting of proportions for all nucleotides and convert in to percentages.
nras.3.percent.reads.f <- rbind(nras.3.a_pc.reads.f, nras.3.c_pc.reads.f, nras.3.g_pc.reads.f, nras.3.t_pc.reads.f)
nras.3.percent.reads.r <- rbind(nras.3.a_pc.reads.r, nras.3.c_pc.reads.r, nras.3.g_pc.reads.r, nras.3.t_pc.reads.r)

# The column names of 'nras.3.percent.reads.f' and 'nras.3.percent.reads.r' are gDNA co-ordinates.  This needs to be changed to cDNA positions.  This is done by making a dataframe of cDNA positions and assigning that as new column names.
# Make a dataframe of cDNA positions for NRAS exon 3. 
nras.3.cDNA <- as.data.frame(290:112)
nras.3.cDNA <- t(nras.3.cDNA)

# Assign the the above column names.
colnames(nras.3.percent.reads.f) <- nras.3.cDNA
colnames(nras.3.percent.reads.r) <- nras.3.cDNA

# Use 'apply' to convert proportions for each nucleotide in to percentages.
nras.3.percent.reads.f <- apply(nras.3.percent.reads.f, 1:2, function(x) round(as.numeric(x*100),2))
nras.3.percent.reads.r <- apply(nras.3.percent.reads.r, 1:2, function(x) round(as.numeric(x*100),2))

# Make a dataframe of 'seqnames', 'strand' and 'nucleotide' and attach to the percent reads dataframe.
# Forward
nras.3.percent.reads.f.info <- nras.3.pup.f.s[1:4,1:3]
nras.3.percent.reads.f.export <- cbind(nras.3.percent.reads.f.info, nras.3.percent.reads.f)
# As NRAS is in the reverse orientation the row names have been reversed too.  However, the 'Nucleotide' column is not reversed
row.names(nras.3.percent.reads.f.export) <- c('T', 'G', 'C', 'A')
# Reverse
nras.3.percent.reads.r.info <- nras.3.pup.r.s[1:4,1:3]
nras.3.percent.reads.r.export <- cbind(nras.3.percent.reads.r.info, nras.3.percent.reads.r)
# As NRAS is in the reverse orientation the row names have been reversed too.  However, the 'Nucleotide' column is not reversed
row.names(nras.3.percent.reads.r.export) <- c('T', 'G', 'C', 'A')


# Output percentage reads to a .csv file.
write.csv(nras.3.percent.reads.f.export, paste(path_to_report, NRAS_3_pc_reads_f, sep = '/') quote = F)
write.csv(nras.3.percent.reads.r.export, paste(path_to_report, NRAS_3_pc_reads_r, sep = '/') quote = F)

# This section rearranges the data into a form suitable for making a stacked bar chart with ggplot2

# melt (reshape2) converts 'wide' data into 'long' data, which is needed for making a stacked bar chart.
# Forward
rownames(nras.3.percent.reads.f) <- c('T', 'G', 'C', 'A')
nras.3.percent.reads.f.melt <- melt(nras.3.percent.reads.f, id.vars = colnames(nras.3.percent.reads.f), value.name = 'Percent')
colnames(nras.3.percent.reads.f.melt) <- c('Nucleotide', 'cDNA_Position', 'Percent')

# Re-order the columns.  This makes the mutation report make more sense.
nras.3.percent.reads.f.melt <- nras.3.percent.reads.f.melt[c('cDNA_Position', 'Nucleotide', 'Percent')]

# Reverse
rownames(nras.3.percent.reads.r) <- c('T', 'G', 'C', 'A')
nras.3.percent.reads.r.melt <- melt(nras.3.percent.reads.r, id.vars = colnames(nras.3.percent.reads.r), value.name = 'Percent')
colnames(nras.3.percent.reads.r.melt) <- c('Nucleotide', 'cDNA_Position', 'Percent')

# Re-order the columns.  This makes the mutation report make more sense.
nras.3.percent.reads.r.melt <- nras.3.percent.reads.r.melt[c('cDNA_Position', 'Nucleotide', 'Percent')]


# First, assign some colours to a palette.  See 'http://www.hypergurl.com/rgbcolorchart.html' for hexadecimal codes.
seq.data.colours <- c('#FF0000', '#FF9900', '#0000FF', '#33CC00')

# STACKED BARCHARTS: FORWARD
# Use 'melted' data to make the barcharts. Start by adding a 'c.' to the cDNA number.  Each cDNA number will appear in the column 4 times, making the number of factors in the column = (number of rows)/4.  This makes the data suitable for a stacked barchart.
nras.3.percent.reads.f.melt.1 <- as.data.frame(paste('c.', nras.3.percent.reads.f.melt$cDNA_Position, sep = ''))
nras.3.percent.reads.f.melt.2 <- cbind(nras.3.percent.reads.f.melt.1, nras.3.percent.reads.f.melt)
colnames(nras.3.percent.reads.f.melt.2) <- c('cDNA_Position', 'Position', 'Nucleotide', 'Percent')

# Make a stacked barchart.  The plot has been assigned to an object to prevent it being printed.
nras.3.fwd.plot <- ggplot(nras.3.percent.reads.f.melt.2, aes(x = cDNA_Position, y = Percent, fill = Nucleotide)) + geom_bar(stat = 'identity') + scale_fill_manual(values = seq.data.colours) + labs(title = paste('Sample number', bam_filename, '\n\nStacked bar chart for forward reads in NRAS exon 3', sep = ' '), x = 'cDNA position', y = 'Percentage allele') + theme(axis.text.x = element_text(angle = 90))
ggsave(paste(path_to_report, NRAS_3_chart_f, sep = '/'))

# STACKED BARCHARTS: REVERSE
# Use 'melted' data to make the barcharts. Start by adding a 'c.' to the cDNA number.  Each cDNA number will appear in the column 4 times, making the number of factors in the column = (number of rows)/4.  This makes the data suitable for a stacked barchart.
nras.3.percent.reads.r.melt.1 <- as.data.frame(paste('c.', nras.3.percent.reads.r.melt$cDNA_Position, sep = ''))
nras.3.percent.reads.r.melt.2 <- cbind(nras.3.percent.reads.r.melt.1, nras.3.percent.reads.r.melt)
colnames(nras.3.percent.reads.r.melt.2) <- c('cDNA_Position', 'Position', 'Nucleotide', 'Percent')

# Make a stacked barchart.
nras.3.rev.plot <- ggplot(nras.3.percent.reads.r.melt.2, aes(x = cDNA_Position, y = Percent, fill = Nucleotide)) + geom_bar(stat = 'identity') + scale_fill_manual(values = seq.data.colours) + labs(title = paste('Sample number', bam_filename, '\n\nStacked bar chart for reverse reads in NRAS exon 3', sep = ' '), x = 'cDNA position', y = 'Percentage allele') + theme(axis.text.x = element_text(angle = 90))
ggsave(paste(path_to_report, NRAS_3_chart_r, sep = '/'))


# Mutation report: FORWARD 
# Add a 'c.' and the nucleotide to the melted data.  This means that the number of rows in the column = number of factors.  This makes the data suitable for a properly annotated mutation report. 
nras.3.percent.reads.f.melt.3 <- as.data.frame(paste('c.', nras.3.percent.reads.f.melt$cDNA_Position, nras.3.percent.reads.f.melt$Nucleotide, sep = ''))
nras.3.percent.reads.f.melt.4 <- cbind(nras.3.percent.reads.f.melt.3, nras.3.percent.reads.f.melt)
colnames(nras.3.percent.reads.f.melt.4) <- c('cDNA_Position', 'Position', 'Nucleotide', 'Percent')

# Filter out the codons of interest plus one either side.  Filter by the 'Position' column.
nras.3.percent.reads.f.melt.60_62 <- filter(nras.3.percent.reads.f.melt.4, Position >= 178, Position <= 186)

# Mutation report: REVERSE
# Add a 'c.' and the nucleotide to the melted data.  This means that the number of rows in the column = number of factors.  This makes the data suitable for a properly annotated mutation report.  
nras.3.percent.reads.r.melt.3 <- as.data.frame(paste('c.', nras.3.percent.reads.r.melt$cDNA_Position, nras.3.percent.reads.r.melt$Nucleotide, sep = ''))
nras.3.percent.reads.r.melt.4 <- cbind(nras.3.percent.reads.r.melt.3, nras.3.percent.reads.r.melt)
colnames(nras.3.percent.reads.r.melt.4) <- c('cDNA_Position', 'Position', 'Nucleotide', 'Percent')

# Filter out the codons of interest plus one either side.  Filter by the 'Position' column.
nras.3.percent.reads.r.melt.60_62 <- filter(nras.3.percent.reads.r.melt.4, Position >= 178, Position <= 186)

# Return positions that have nucleotides >=5% and <=95%.  This should detect nucleotide changes greater than 5%.
nras.3.muts_5_f <- filter(nras.3.percent.reads.f.melt, Percent >= 5.00, Percent <= 95.00); nras.3.muts_5_f
nras.3.muts_5_r <- filter(nras.3.percent.reads.r.melt, Percent >= 5.00, Percent <= 95.00); nras.3.muts_5_r

# Export the data to a mutation report.
write.table(paste('Sample number', bam_samplename, '\n\n', sep = ' '), paste(path_to_report, NRAS_3_mutations, sep = '/'), row.names = F, col.names = F, quote = F, append = T)
write.table(nras.3.muts_5_f, paste(path_to_report, NRAS_3_mutations, sep = '/'), row.names = F, quote = F, append = T)
write.table('\n', paste(path_to_report, NRAS_3_mutations, sep = '/'), row.names = F, col.names = F, quote = F, append = T)
write.table(nras.3.muts_5_r, paste(path_to_report, NRAS_3_mutations, sep = '/'), row.names = F, quote = F, append = T)
