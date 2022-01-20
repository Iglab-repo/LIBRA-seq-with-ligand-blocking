#########################################################

# "Copyright 2022 Vanderbilt University. All rights reserved"

#########################################################

#!/usr/bin/env Rscript


# read parameters file
args = commandArgs(trailingOnly = T)
if (length(args)==0) {
  stop("Parameters file must be supplied\n", call. = F)
} else if (length(args)==1) {
  source(args[1])
}


# create analysis directory filepath and set as working directory
analysis.dir <- file.path(OUTPUT.DIR, ANALYSIS.ID)
dir.create(analysis.dir, showWarnings = F)
setwd(analysis.dir)


### CELLRANGER VDJ ANALYSIS ###

if(CELLRANGER.VDJ.ANALYSIS) {
  
  # run cellranger VDJ analysis
  cellranger.cmd <- paste0(CELLRANGER.PATH, " vdj",
                           " --id=VDJ_analysis",
                           " --reference=", VDJ.REF,
                           " --fastqs=", VDJ.FASTQ.DIR,
                           " --localmem=", 64,
                           " >> pipeline.log 2>&1")

  # write cellranger command to logfile
  sink(file = "pipeline.log", append = T, split = T)
  cat("Running cellranger VDJ analysis:\n", cellranger.cmd, "\n\n", sep = "")
  sink()

  # output cellranger command to system
  system(cellranger.cmd)
  
  # create IMGT subdirectory in which to place IMGT output
  dir.create("IMGT", showWarnings = F)
  invisible(file.create("IMGT/PLACE_IMGT_OUTPUT_HERE", showWarnings = F))
  invisible(file.copy("VDJ_analysis/outs/filtered_contig.fasta", "IMGT/"))
}


### CELLRANGER AG ANALYSIS ###

if(CELLRANGER.AG.ANALYSIS) {
  
  
  # create antigens CSV file for cellranger
  ANTIGENS <- matrix(unlist(ANTIGENS), nrow = length(ANTIGENS), byrow = T)
  antigens.csv <- data.frame(id=ANTIGENS[,1], name=ANTIGENS[,1], read="R2", pattern="5P(BC)", sequence=ANTIGENS[,2], feature_type="Antibody Capture")
  antigens.path <- file.path(analysis.dir, "antigens.csv")
  write.table(antigens.csv, antigens.path, quote = F, sep = ",", row.names = F)
  
  # create library CSV file for cellranger
  library.csv <- data.frame(fastqs=AG.FASTQ.DIR, sample=AG.FASTQ.PREFIX, library_type="Antibody Capture")
  library.path <- file.path(analysis.dir, "library.csv")
  write.table(library.csv, library.path, quote = F, sep = ",", row.names = F)
  
  # run cellranger antigen counts analysis
  cellranger.cmd <- paste0(CELLRANGER.PATH, " count",
                           " --id=AG_analysis",
                           " --transcriptome=", TRANSCRIPTOME.REF,
                           " --feature-ref=", antigens.path,
                           " --libraries=", library.path,
                           " --localmem=", 64,
                           " >> pipeline.log 2>&1")

  # write cellranger command to logfile
  sink(file = "pipeline.log", append = T, split = T)
  cat("\nRunning cellranger antigen counts analysis:\n", cellranger.cmd, "\n\n", sep = "")
  sink()

  # output cellranger command to system
  system(cellranger.cmd)
}


### CELLRANGER RNA ANALYSIS ###

if(CELLRANGER.RNA.ANALYSIS) {
  
  # run cellranger RNA counts analysis
  cellranger.cmd <- paste0(CELLRANGER.PATH, " count",
                           " --id=RNA_analysis",
                           " --transcriptome=", TRANSCRIPTOME.REF,
                           " --fastqs=", RNA.FASTQ.DIR,
                           " --localmem=", 64,
                           " >> pipeline.log 2>&1") 

  # write cellranger command to logfile
  sink(file = "pipeline.log", append = T, split = T)
  cat("\nRunning cellranger RNA counts analysis:\n", cellranger.cmd, "\n\n", sep = "")
  sink()

  # output cellranger command to system
  system(cellranger.cmd)
}

#########################################################

# "Copyright 2022 Vanderbilt University. All rights reserved"

#########################################################
