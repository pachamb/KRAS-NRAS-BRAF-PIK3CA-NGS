
# This analyses short read, next generation sequencing data for the human BRAF gene generated from an Illumina Genome Analyser using the Qiagen "GeneRead DNAseq Targeted Panels V2 Human Tumor Actionable Mutations Panel" (Product no. 181900, Cat. no. NGHS-201X).  It should also work with data generated from an Illumina Genome Analyser using and any PCR primers that amplify the regions of interest. Gene alignment and numbering is for human genome build hg19.
# It generates three outputs:
# 1. csv files with the percentage reads for each nucleotide at each position.
# 2. pdf files of stacked bar charts showing the percentage of each nucleotide at each position.
# 3. A txt file containing any mutations detected.


# BRAF EXON 15

# bam file scanning parameters.  These are for BRAF exon 15.
braf.15.param <- ScanBamParam(which = GRanges('chr7', IRanges(start = 140453075, end = 140453193)))

# Pileup parameters.
braf.15.pup.param <- PileupParam(max_depth = 10000, min_base_quality = 20, min_nucleotide_depth = 10, distinguish_strands = T)

# Perform the pileup operation.
braf.15.pup <- pileup(sample_files, scanBamParam = braf.15.param, pileupParam = braf.15.pup.param)

# Remove column called 'which_label' and keep the rest.
braf.15.pup.2 <- select(braf.15.pup, -which_label)

# Make two dataframes one for the '+' strand and one for the '-' strand.
braf.15.pup.f <- filter(braf.15.pup, strand == '+')
braf.15.pup.r <- filter(braf.15.pup, strand == '-')

# If there are >=2 nucleotides at one position the pileup function generates >=2 rows for that position.  This is particularly useless for mutation screening where the whole aim is to detect positions where there is more than 1 nucleotide...
# spread (tidyr) generates a dataframe that looks very like a consensus matrix.  spread (tidyr) needs to be told the dataframe it's working on, a 'key', i.e. name of the column whose values will be used as column headings and a 'value', i.e. name of the column whose values will populate the cells. 'fill' fills in the missing values.  Most importantly, it combines data for one position in to 1 column.
braf.15.pup.f.s <- spread(braf.15.pup.f, pos, count, fill = 0)
braf.15.pup.r.s <- spread(braf.15.pup.r, pos, count, fill = 0)

# Sum the counts in columns of the above dataframes.
braf.15.pup.f.s.sum <- apply(braf.15.pup.f.s[,5:123], 2, sum)
braf.15.pup.r.s.sum <- apply(braf.15.pup.r.s[,5:123], 2, sum)

# rbind results of apply operation to 'consensus matrix'.  This will lose the seqnames, strand and nucleotide columns.
braf.15.pup.2.f.s.2 <- rbind(braf.15.pup.f.s[,5:123], braf.15.pup.f.s.sum)
braf.15.pup.2.r.s.2 <- rbind(braf.15.pup.r.s[,5:123], braf.15.pup.r.s.sum)

# Calculate the proportion of counts for each nucleotide. 
braf.15.a_pc.reads.f <- braf.15.pup.2.f.s.2[1,]/braf.15.pup.2.f.s.2[nrow(braf.15.pup.2.f.s.2),]
braf.15.c_pc.reads.f <- braf.15.pup.2.f.s.2[2,]/braf.15.pup.2.f.s.2[nrow(braf.15.pup.2.f.s.2),]
braf.15.g_pc.reads.f <- braf.15.pup.2.f.s.2[3,]/braf.15.pup.2.f.s.2[nrow(braf.15.pup.2.f.s.2),]
braf.15.t_pc.reads.f <- braf.15.pup.2.f.s.2[4,]/braf.15.pup.2.f.s.2[nrow(braf.15.pup.2.f.s.2),]
braf.15.a_pc.reads.r <- braf.15.pup.2.r.s.2[1,]/braf.15.pup.2.r.s.2[nrow(braf.15.pup.2.r.s.2),]
braf.15.c_pc.reads.r <- braf.15.pup.2.r.s.2[2,]/braf.15.pup.2.r.s.2[nrow(braf.15.pup.2.r.s.2),]
braf.15.g_pc.reads.r <- braf.15.pup.2.r.s.2[3,]/braf.15.pup.2.r.s.2[nrow(braf.15.pup.2.r.s.2),]
braf.15.t_pc.reads.r <- braf.15.pup.2.r.s.2[4,]/braf.15.pup.2.r.s.2[nrow(braf.15.pup.2.r.s.2),]


# Make a single dataframe consisting of proportions for all nucleotides and convert in to percentages.
braf.15.percent.reads.f <- rbind(braf.15.a_pc.reads.f, braf.15.c_pc.reads.f, braf.15.g_pc.reads.f, braf.15.t_pc.reads.f)
braf.15.percent.reads.r <- rbind(braf.15.a_pc.reads.r, braf.15.c_pc.reads.r, braf.15.g_pc.reads.r, braf.15.t_pc.reads.r)

# The column names of 'braf.15.percent.reads.f' and 'braf.15.percent.reads.r' are gDNA co-ordinates.  This needs to be changed to cDNA positions.  This is done by making a dataframe of cDNA positions and assigning that as new column names.
# Make a dataframe of cDNA positions for BRAF exon 15. 
braf.15.cDNA <- as.data.frame(1860:1742)
braf.15.cDNA <- t(braf.15.cDNA)

# Assign the the above column names.
colnames(braf.15.percent.reads.f) <- braf.15.cDNA
colnames(braf.15.percent.reads.r) <- braf.15.cDNA

# Use 'apply' to convert proportions for each nucleotide in to percentages.
braf.15.percent.reads.f <- apply(braf.15.percent.reads.f, 1:2, function(x) round(as.numeric(x*100),2))
braf.15.percent.reads.r <- apply(braf.15.percent.reads.r, 1:2, function(x) round(as.numeric(x*100),2))

# Make a dataframe of 'seqnames', 'strand' and 'nucleotide' and attach to the percent reads dataframe.
# Forward
braf.15.percent.reads.f.info <- braf.15.pup.f.s[1:4,1:3]
braf.15.percent.reads.f.export <- cbind(braf.15.percent.reads.f.info, braf.15.percent.reads.f)
# As BRAF is in the reverse orientation the row names have been reversed too.  However, the 'Nucleotide' column is not reversed.
row.names(braf.15.percent.reads.f.export) <- c('T', 'G', 'C', 'A')
# Reverse
braf.15.percent.reads.r.info <- braf.15.pup.r.s[1:4,1:3]
braf.15.percent.reads.r.export <- cbind(braf.15.percent.reads.r.info, braf.15.percent.reads.r)
# As BRAF is in the reverse orientation the row names have been reversed too.  However, the 'Nucleotide' column is not reversed.
row.names(braf.15.percent.reads.r.export) <- c('T', 'G', 'C', 'A')


# Output percentage reads to a .csv file.
write.csv(braf.15.percent.reads.f.export, paste(path_to_report, BRAF_pc_reads_f, sep = '/'), quote = F)
write.csv(braf.15.percent.reads.r.export, paste(path_to_report, BRAF_pc_reads_r, sep = '/'), quote = F)

# This section rearranges the data into a form suitable for making a stacked bar chart with ggplot2
# melt (reshape2) converts 'wide' data into 'long' data, which is needed for making a stacked bar chart.
# Forward
rownames(braf.15.percent.reads.f) <- c('T', 'G', 'C', 'A')
braf.15.percent.reads.f.melt <- melt(braf.15.percent.reads.f, id.vars = colnames(braf.15.percent.reads.f), value.name = 'Percent')
colnames(braf.15.percent.reads.f.melt) <- c('Nucleotide', 'cDNA_Position', 'Percent')

# Re-order the columns.  This makes the mutation report make more sense.
braf.15.percent.reads.f.melt <- braf.15.percent.reads.f.melt[c('cDNA_Position', 'Nucleotide', 'Percent')]

# Reverse
rownames(braf.15.percent.reads.r) <- c('T', 'G', 'C', 'A')
braf.15.percent.reads.r.melt <- melt(braf.15.percent.reads.r, id.vars = colnames(braf.15.percent.reads.r), value.name = 'Percent')
colnames(braf.15.percent.reads.r.melt) <- c('Nucleotide', 'cDNA_Position', 'Percent')

# Re-order the columns.  This makes the mutation report make more sense.
braf.15.percent.reads.r.melt <- braf.15.percent.reads.r.melt[c('cDNA_Position', 'Nucleotide', 'Percent')]

# Assign some colours to a palette.  See 'http://www.hypergurl.com/rgbcolorchart.html' for hexadecimal codes.
seq.data.colours <- c('#FF0000', '#FF9900', '#0000FF', '#33CC00')

# STACKED BARCHARTS: FORWARD
# Use 'melted' data to make the barcharts. Start by adding a 'c.' to the cDNA number.  Each cDNA number will appear in the column 4 times, making the number of factors in the column = (number of rows)/4.  This makes the data suitable for a stacked barchart.
braf.15.percent.reads.f.melt.1 <- as.data.frame(paste('c.', braf.15.percent.reads.f.melt$cDNA_Position, sep = ''))
braf.15.percent.reads.f.melt.2 <- cbind(braf.15.percent.reads.f.melt.1, braf.15.percent.reads.f.melt)
colnames(braf.15.percent.reads.f.melt.2) <- c('cDNA_Position', 'Position', 'Nucleotide', 'Percent')

# Make a stacked barchart.  The plot has been assigned to an object to prevent it being printed.
braf.15.fwd.plot <- ggplot(braf.15.percent.reads.f.melt.2, aes(x = cDNA_Position, y = Percent, fill = Nucleotide)) + geom_bar(stat = 'identity') + scale_fill_manual(values = seq.data.colours) + labs(title = paste('Sample number', bam_filename, '\n\nStacked bar chart for forward reads in BRAF exon 15', sep = ' '), x = 'cDNA position', y = 'Percentage allele') + theme(axis.text.x = element_text(angle = 90))
ggsave(paste(path_to_report, BRAF_chart_f, sep = '/'))

# STACKED BARCHARTS: REVERSE
# Use 'melted' data to make the barcharts. Start by adding a 'c.' to the cDNA number.  Each cDNA number will appear in the column 4 times, making the number of factors in the column = (number of rows)/4.  This makes the data suitable for a stacked barchart.
braf.15.percent.reads.r.melt.1 <- as.data.frame(paste('c.', braf.15.percent.reads.r.melt$cDNA_Position, sep = ''))
braf.15.percent.reads.r.melt.2 <- cbind(braf.15.percent.reads.r.melt.1, braf.15.percent.reads.r.melt)
colnames(braf.15.percent.reads.r.melt.2) <- c('cDNA_Position', 'Position', 'Nucleotide', 'Percent')

# Make a stacked barchart.  The plot has been assigned to an object to prevent it being printed.
braf.15.rev.plot <- ggplot(braf.15.percent.reads.r.melt.2, aes(x = cDNA_Position, y = Percent, fill = Nucleotide)) + geom_bar(stat = 'identity') + scale_fill_manual(values = seq.data.colours) + labs(title = paste('Sample number', bam_filename, '\n\nStacked bar chart for reverse reads in BRAF exon 15', sep = ' '), x = 'cDNA position', y = 'Percentage allele') + theme(axis.text.x = element_text(angle = 90))
ggsave(paste(path_to_report, BRAF_chart_r, sep = '/'))


# Mutation report: FORWARD 
# Add a 'c.' and the nucleotide to the melted data.  This means that the number of rows in the column = number of factors.  This makes the data suitable for a properly annotated mutation report. 
braf.15.percent.reads.f.melt.3 <- as.data.frame(paste('c.', braf.15.percent.reads.f.melt$cDNA_Position, braf.15.percent.reads.f.melt$Nucleotide, sep = ''))
braf.15.percent.reads.f.melt.4 <- cbind(braf.15.percent.reads.f.melt.3, braf.15.percent.reads.f.melt)
colnames(braf.15.percent.reads.f.melt.4) <- c('cDNA_Position', 'Position', 'Nucleotide', 'Percent')

# Mutation report: REVERSE
# Add a 'c.' and the nucleotide to the melted data.  This means that the number of rows in the column = number of factors.  This makes the data suitable for a properly annotated mutation report. 
braf.15.percent.reads.r.melt.3 <- as.data.frame(paste('c.', braf.15.percent.reads.r.melt$cDNA_Position, braf.15.percent.reads.r.melt$Nucleotide, sep = ''))
braf.15.percent.reads.r.melt.4 <- cbind(braf.15.percent.reads.r.melt.3, braf.15.percent.reads.r.melt)
colnames(braf.15.percent.reads.r.melt.4) <- c('cDNA_Position', 'Position', 'Nucleotide', 'Percent')

# Return positions that have nucleotides >=5% and <=95%.  This should detect nucleotide changes greater than 5%.
braf.15.muts_5_f <- filter(braf.15.percent.reads.f.melt, Percent >= 5.00, Percent <= 95.00); braf.15.muts_5_f
braf.15.muts_5_r <- filter(braf.15.percent.reads.r.melt, Percent >= 5.00, Percent <= 95.00); braf.15.muts_5_r

# Export the data to a mutation report.
write.table(paste('Sample number', bam_samplename, '\n\n', sep = ' '), paste(path_to_report, BRAF_mutations, sep = '/'), row.names = F, col.names = F, quote = F, append = T)
write.table(braf.15.muts_5_f, paste(path_to_report, BRAF_mutations, sep = '/'), row.names = F, quote = F, append = T)
write.table('\n', paste(path_to_report, BRAF_mutations, sep = '/'), row.names = F, col.names = F, quote = F, append = T)
write.table(braf.15.muts_5_r, paste(path_to_report, BRAF_mutations, sep = '/'), row.names = F, quote = F, append = T)


