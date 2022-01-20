# LIBRA-seq-with-ligand-blocking
Scripts used in LIBRA-seq with ligand blocking data analysis

There are three R scripts used in the LIBRA-seq data analysis that includes
i) parameters_5317
ii) LIBRAseq_analysis_IMGT_part1.R
iii) LIBRAseq_analysis_IMGT_part2.R

Description of the files
-------------------------
parameters.R: In this script, all the filanmes and paths should up added including the paths of Cellranger and Change-O tools. This script should be executed first to set the file and program paths

LIBRAseq_analysis_IMGT_part1.R: This is the main script that process the sequencing data using Cell ranger.
LIBRAseq_analysis_IMGT_part2.R: This script is used to combine the heavy and light chains based on the cell barcodes and to compute LIBRA-seq scores for each antigen.

For further details or help, please contact the corresponding author

