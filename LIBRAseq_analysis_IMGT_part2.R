#########################################################

# "Copyright 2022 Vanderbilt University. All rights reserved"

#########################################################

#!/usr/bin/env Rscript
library(data.table)
library(Seurat)
library(tidyr)
library(ggplot2)
library(gridExtra)


# read parameters file
args = commandArgs(trailingOnly = T)
if (length(args)==0) {
  stop("Parameters file must be supplied\n", call. = F)
} else if (length(args)==1) {
  source(args[1])
}


# create analysis directory filepath
analysis.dir <- file.path(OUTPUT.DIR, ANALYSIS.ID)
setwd(analysis.dir)


# search for .txz file in IMGT directory
imgt.output <- list.files(path = "IMGT", pattern = ".txz", full.names = T)

if (length(imgt.output) == 0) {
  stop("TXZ file not found in IMGT subdirectory\n", call. = F)
} else if (length(imgt.output) > 1) {
  stop("Multiple TXZ files found in IMGT subdirectory\n", call. = F)
}


# convert IMGT output to CHANGE-O format
changeo.cmd <- paste(CHANGEO.MAKEDB.PATH, "imgt -i", imgt.output, "-s VDJ/filtered_contig.fasta --10x VDJ/filtered_contig_annotations.csv --extended -o IMGT/IMGT_output.tab")


# write CHANGE-O command to logfile
sink(file = "pipeline.log", append = T, split = T)
cat("\nConverting IMGT output to CHANGE-O format:\n", changeo.cmd, "\n\n", sep = "")
sink()


# output CHANGE-O command to system
system(changeo.cmd)


# parse heavy chain data from CHANGE-O table
changeo.cmd <- paste(CHANGEO.PARSEDB.PATH, "select -d IMGT/IMGT_output.tab -f LOCUS -u \"IGH\" --logic all --regex -o IMGT/IMGT_heavy_chains.tab")


# write CHANGE-O command to logfile
sink(file = "pipeline.log", append = T, split = T)
cat("\nParsing heavy chain data from CHANGE-O table:\n", changeo.cmd, "\n\n", sep = "")
sink()


# output CHANGE-O command to system
system(changeo.cmd)

# cluster heavy chain sequences into clonal groups
changeo.cmd <- paste(CHANGEO.DEFINECLONES.PATH, "-d IMGT/IMGT_heavy_chains.tab --mode gene --act first --model ham --dist 0.20 --norm len -o IMGT/IMGT_heavy_chains_with_clonal_groups.tab")


# write CHANGE-O command to logfile
sink(file = "pipeline.log", append = T, split = T)
cat("\nClustering heavy chain sequences into clonal groups:\n", changeo.cmd, "\n\n", sep = "")
sink()


# output CHANGE-O command to system
system(changeo.cmd)

# parse light chain data from CHANGE-O table
changeo.cmd <- paste(CHANGEO.PARSEDB.PATH, "select -d IMGT/IMGT_output.tab -f LOCUS -u \"IG[LK]\" --logic all --regex -o IMGT/IMGT_light_chains.tab")


# write CHANGE-O command to logfile
sink(file = "pipeline.log", append = T, split = T)
cat("\nParsing light chain data from CHANGE-O table:\n", changeo.cmd, "\n\n", sep = "")
sink()


# output CHANGE-O command to system
system(changeo.cmd)


# cluster light chain sequences into clonal groups
changeo.cmd <- paste(CHANGEO.DEFINECLONES.PATH, "-d IMGT/IMGT_light_chains.tab --mode gene --act first --model ham --dist 0.20 --norm len -o IMGT/IMGT_light_chains_with_clonal_groups.tab")


# write CHANGE-O command to logfile
sink(file = "pipeline.log", append = T, split = T)
cat("\nClustering light chain sequences into clonal groups:\n", changeo.cmd, "\n\n", sep = "")
sink()


# output CHANGE-O command to system
system(changeo.cmd)


# read in heavy chain table
heavy.chains <- fread("IMGT/IMGT_heavy_chains_with_clonal_groups.tab", select = c("SEQUENCE_ID", "SEQUENCE_INPUT", "FUNCTIONAL", "IN_FRAME",	"INDELS", "LOCUS", "V_CALL", "D_CALL",	"J_CALL",	"SEQUENCE_VDJ",	"SEQUENCE_IMGT",	"JUNCTION",	"JUNCTION_LENGTH",	"GERMLINE_IMGT",	"V_SCORE",	"V_IDENTITY",	"J_SCORE",	"J_IDENTITY",	"FWR1_IMGT",	"FWR2_IMGT",	"FWR3_IMGT",	"FWR4_IMGT",	"CDR1_IMGT",	"CDR2_IMGT",	"CDR3_IMGT",	"C_CALL",	"CONSCOUNT",	"UMICOUNT",	"JUNCTION_10X_AA", "CLONE"))


# perform table manipulations
heavy.chains <- separate(data = heavy.chains, col = SEQUENCE_ID, into = c("BARCODE", "CONTIG"), sep = "-1_")
barcode.counts <- data.frame(count = length(unique(heavy.chains$BARCODE)), description = "unique heavy chain barcodes")
heavy.chains <- heavy.chains[heavy.chains$FUNCTIONAL == "T"]
barcode.counts <- rbind(barcode.counts, data.frame(count = length(unique(heavy.chains$BARCODE)), description = "unique heavy chain barcodes after filtering out cells with non-functional heavy chains"))
heavy.chains$V_CALL <- sub("Homsap ", "", sub("\\*.*", "", heavy.chains$V_CALL))
heavy.chains$D_CALL <- sub("Homsap ", "", sub("\\*.*", "", heavy.chains$D_CALL))
heavy.chains$J_CALL <- sub("Homsap ", "", sub("\\*.*", "", heavy.chains$J_CALL))


# remove duplicate barcodes with differing V genes
barcode.dups <- heavy.chains[duplicated(heavy.chains$BARCODE)]$BARCODE
v_call.dups <- subset(heavy.chains, subset = BARCODE %in% barcode.dups, select = c(BARCODE, V_CALL))
barcode.dups <- barcode.dups[!barcode.dups %in% v_call.dups[duplicated(v_call.dups)]$BARCODE]
heavy.chains <- subset(heavy.chains, subset = !BARCODE %in% barcode.dups)
barcode.counts <- rbind(barcode.counts, data.frame(count = length(unique(heavy.chains$BARCODE)), description = "unique heavy chain barcodes after filtering out cells with multiple heavy chains"))


# read in light chain table
light.chains <- fread("IMGT/IMGT_light_chains_with_clonal_groups.tab", select = c("SEQUENCE_ID",	"SEQUENCE_INPUT",	"FUNCTIONAL",	"IN_FRAME",	"INDELS",	"LOCUS",	"V_CALL",	"J_CALL",	"SEQUENCE_VDJ",	"SEQUENCE_IMGT", "JUNCTION",	"JUNCTION_LENGTH",	"GERMLINE_IMGT",	"V_SCORE",	"V_IDENTITY",	"J_SCORE",	"J_IDENTITY",	"FWR1_IMGT",	"FWR2_IMGT",	"FWR3_IMGT",	"FWR4_IMGT",	"CDR1_IMGT",	"CDR2_IMGT",	"CDR3_IMGT",	"C_CALL",	"CONSCOUNT", "JUNCTION_10X_AA",	"UMICOUNT"))


# perform table manipulations
light.chains <- separate(data = light.chains, col = SEQUENCE_ID, into = c("BARCODE", "CONTIG"), sep = "-1_")
barcode.counts <- rbind(barcode.counts, data.frame(count = length(unique(light.chains$BARCODE)), description = "unique light chain barcodes"))
light.chains <- light.chains[light.chains$FUNCTIONAL == "T"]
barcode.counts <- rbind(barcode.counts, data.frame(count = length(unique(light.chains$BARCODE)), description = "unique light chain barcodes after filtering out cells with non-functional light chains"))
light.chains$V_CALL <- sub("Homsap ", "", sub("\\*.*", "", light.chains$V_CALL))
light.chains$D_CALL <- sub("Homsap ", "", sub("\\*.*", "", light.chains$D_CALL))
light.chains$J_CALL <- sub("Homsap ", "", sub("\\*.*", "", light.chains$J_CALL))


# merge heavy and light chain tables
heavy.light.chains <- merge(heavy.chains, light.chains, by = "BARCODE", suffixes = c(".H", ".L"))
heavy.light.chains <- merge(heavy.light.chains, data.table(table(BARCODE=heavy.light.chains$BARCODE)), by = "BARCODE")
heavy.light.chains <- unite(heavy.light.chains, col = "CLONE", c(CLONE.H, CLONE.L), remove = F)
barcode.counts <- rbind(barcode.counts, data.frame(count = length(unique(heavy.light.chains$BARCODE)), description = "unique paired heavy-light chain barcodes"))


# read in cellranger antigen count matrix
ag.counts <- t(as.matrix(Read10X("AG_analysis/outs/raw_feature_bc_matrix", )))
ag.counts <- ag.counts[rowSums(ag.counts) > 0, ]
barcode.counts <- rbind(barcode.counts, data.frame(count = nrow(ag.counts), description = "antigen barcodes after filtering out cells with a total UMI count of zero"))


# normalize unfiltered antigen count data
ag.geomean <- apply(ag.counts, 1, function(x) prod(x + 1)^(1 / ncol(ag.counts)))
ag.clr <- log((ag.counts + 1) / ag.geomean)
ag.clr.z <- apply(ag.clr, 2, function(x) (x - mean(x)) / sd(x))
for(i in 1:ncol(ag.clr.z)) {ag.clr.z[ag.counts[, i] == 0, i] <- min(ag.clr.z[, i])}


# create unfiltered antigen count data tables
unfiltered.ag.counts.table <- data.table(BARCODE=rownames(ag.counts), ag.counts, row.names = NULL)
unfiltered.ag.clr.z.table <- data.table(BARCODE=rownames(ag.clr.z), ag.clr.z, row.names = NULL)


# filter antigen count matrix
ag.counts <- ifelse(ag.counts < 4, 0, ag.counts)
ag.counts <- ag.counts[rowSums(ag.counts) > 0, ] 
barcode.counts <- rbind(barcode.counts, data.frame(count = nrow(ag.counts), description = "antigen barcodes after setting UMI counts less than four to zero and re-filtering"))


# determine overlapping barcodes and subset data
barcodes <- intersect(heavy.light.chains$BARCODE, rownames(ag.counts))
heavy.light.chains <- subset(heavy.light.chains, subset = BARCODE %in% barcodes)
ag.counts <- subset(ag.counts, subset = rownames(ag.counts) %in% barcodes)
barcode.counts <- rbind(barcode.counts, data.frame(count = length(barcodes), description = "overlapping barcodes between paired heavy-light chains and antigen counts"))


# normalize antigen count data
ag.geomean <- apply(ag.counts, 1, function(x) prod(x + 1)^(1 / ncol(ag.counts)))
ag.clr <- log((ag.counts + 1) / ag.geomean)
ag.clr.z <- apply(ag.clr, 2, function(x) (x - mean(x)) / sd(x))
for(i in 1:ncol(ag.clr.z)) {ag.clr.z[ag.counts[, i] == 0, i] <- min(ag.clr.z[, i])}


# binarize antigen count data
ag.binary <- ifelse(ag.clr.z >= 1, 1, 0)


# create antigen count data tables
ag.counts.table <- data.table(BARCODE=rownames(ag.counts), ag.counts, row.names = NULL)
ag.clr.z.table <- data.table(BARCODE=rownames(ag.clr.z), ag.clr.z, row.names = NULL)
ag.binary.table <- data.table(BARCODE=rownames(ag.binary), ag.binary, row.names = NULL)


# append antigen counts and LIBRA-seq scores to paired heavy-light chain table
heavy.light.chains <- merge(heavy.light.chains, ag.counts.table, by = "BARCODE")
heavy.light.chains <- merge(heavy.light.chains, ag.clr.z.table, by = "BARCODE", suffixes = c("", ".LSS"))


histograms <- list()
ggdata1 <- melt(ag.counts.table, id.vars = "BARCODE", variable.name = "antigen", value.name = "count")
histograms[[1]] <- ggplot(ggdata1, aes(x = count, fill = antigen)) + geom_histogram(alpha = 0.5, binwidth = 1) +
  facet_wrap(facets = vars(antigen), scales = "free") + xlim(-1,51) + scale_fill_discrete(name = "Antigen") +
  labs(title = "UMI count distributions", x = "UMI count", y = "Frequency") + theme_classic() +
  theme(legend.position="none")
ggdata2 <- melt(unfiltered.ag.counts.table, id.vars = "BARCODE", variable.name = "antigen", value.name = "count")
histograms[[2]] <- ggplot(ggdata2, aes(x = count, fill = antigen)) + geom_histogram(alpha = 0.5, binwidth = 1) +
  facet_wrap(facets = vars(antigen), scales = "free") + xlim(-1,51) + scale_fill_discrete(name = "Antigen") +
  labs(title = "UMI count distributions (unfiltered)", x = "UMI count", y = "Frequency") + theme_classic() +
  theme(legend.position="none")
ggdata3 <- melt(ag.clr.z.table, id.vars = "BARCODE", variable.name = "antigen", value.name = "count")
histograms[[3]] <- ggplot(ggdata3, aes(x = count, fill = antigen)) + geom_histogram(alpha = 0.5, binwidth = 0.1) +
  facet_wrap(facets = vars(antigen), scales = "free") + scale_y_log10() + scale_fill_discrete(name = "Antigen") +
  labs(title = "LIBRA-seq score distributions", x = "LIBRA-seq score", y = "Frequency") + theme_classic() +
  theme(legend.position="none")
ggdata4 <- melt(unfiltered.ag.clr.z.table, id.vars = "BARCODE", variable.name = "antigen", value.name = "count")
histograms[[4]] <- ggplot(ggdata4, aes(x = count, fill = antigen)) + geom_histogram(alpha = 0.5, binwidth = 0.1) +
  facet_wrap(facets = vars(antigen), scales = "free") + scale_y_log10() + scale_fill_discrete(name = "Antigen") +
  labs(title = "LIBRA-seq score distributions (unfiltered)", x = "LIBRA-seq score", y = "Frequency") + theme_classic() +
  theme(legend.position="none")


scatterplots <- list()
antigen.pairs <- t(combn(2:ncol(ag.clr.z.table), 2))
isotypes <- unique(subset(heavy.light.chains, select = c("BARCODE", "C_CALL.H")))
isotypes$C_CALL.H <- gsub('[0-9]+', '', isotypes$C_CALL.H)
isotypes[isotypes$C_CALL.H == ""] <- "NA"
ag.clr.z.table <- merge(ag.clr.z.table, isotypes, by = "BARCODE")
for(i in 1:nrow(antigen.pairs)) {
  scatterplots[[i]] <- ggplot(ag.clr.z.table, aes_string(x = colnames(ag.clr.z.table)[antigen.pairs[i,1]], 
                                                         y = colnames(ag.clr.z.table)[antigen.pairs[i,2]],
                                                         color = "C_CALL.H")) + 
    geom_point(size = 0.5) + geom_hline(yintercept = 1, color = "darkred") + 
    geom_vline(xintercept = 1, color = "darkred") + 
    scale_colour_manual(values = c("seagreen3", "dodgerblue3", "firebrick3", "black")) + 
    theme_classic(base_size = 8) + 
    theme(legend.title = element_blank(), legend.position = "top", legend.margin = margin(0,0,0,0), legend.box.margin = margin(0,-10,-10,-10))
}
      

# generate PDF of plots
pdf(paste0("plots", format(Sys.time(), "_%Y%m%d_%H%M%S"), ".pdf"), width = 8.5, height = 11)
marrangeGrob(histograms, nrow = 2, ncol = 1, top = "", bottom = "", left = "", right = "", layout_matrix = rbind(1,2))
marrangeGrob(scatterplots, nrow = 5, ncol = 3, top = "", bottom = "", left = "", right = "", layout_matrix = matrix(1:15, ncol = 3, byrow = T))
dev.off()


# write data tables to files
timestamp <- format(Sys.time(), "%Y%m%d_%H%M%S")
write.table(heavy.chains, paste0("heavy_chains_", timestamp, ".tab"), quote = F, sep = "\t", row.names = F)
write.table(light.chains, paste0("light_chains_", timestamp, ".tab"), quote = F, sep = "\t", row.names = F)
write.table(heavy.light.chains, paste0("paired_heavy_light_chains_", timestamp, ".tab"), quote = F, sep = "\t", row.names = F)
write.table(ag.counts.table, paste0("AG_umi_counts_", timestamp, ".tab"), quote = F, sep = "\t", row.names = F)
write.table(ag.clr.z.table, paste0("AG_libraseq_scores_", timestamp, ".tab"), quote = F, sep = "\t", row.names = F)
write.table(ag.binary.table, paste0("AG_binary_scores_", timestamp, ".tab"), quote = F, sep = "\t", row.names = F)
write.table(barcode.counts, paste0("barcode_counts_", timestamp, ".tab"), quote = F, sep = " ", row.names = F, col.names = F)

#########################################################

# "Copyright 2022 Vanderbilt University. All rights reserved"

#########################################################

