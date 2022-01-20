#########################################################

# "Copyright 2022 Vanderbilt University. All rights reserved"

#########################################################

### ANALYSIS ID ###
ANALYSIS.ID = '5317'

### PIPELINE DIRECTORY PATHS ###
OUTPUT.DIR = '/hdd3/andrea'

### CELLRANGER VDJ ANALYSIS PARAMETERS ###
CELLRANGER.VDJ.ANALYSIS = TRUE
VDJ.FASTQ.DIR = '/hdd3/vantage_data/5317/VDJ'
VDJ.REF = '/hdd3/cellranger/GRCh38_VDJ'

### CELLRANGER AG ANALYSIS PARAMETERS ###
CELLRANGER.AG.ANALYSIS = TRUE
AG.FASTQ.DIR = '/hdd3/vantage_data/5317/AG'
AG.FASTQ.PREFIX = '5317-KK-1_Hash'
TRANSCRIPTOME.REF = '/hdd3/cellranger/GRCh38'
ANTIGENS = list()
ANTIGENS[[1]] = c("CoV_HP6", "GACAAGTGATCTGCA")
ANTIGENS[[2]] = c("FLU_HANC99", "TCATTTCCTCCGATT")
ANTIGENS[[3]] = c("HIV_ZM197", "TACGCCTATAACTTG")
ANTIGENS[[4]] = c("ACE2", "CTTCACTCTGTCAGG")


### CELLRANGER RNA ANALYSIS PARAMETERS ###
CELLRANGER.RNA.ANALYSIS = FALSE
RNA.FASTQ.DIR = '/hdd3/vantage_data/5317/GEX'
TRANSCRIPTOME.REF = '/hdd3/cellranger/GRCh38'

### PROGRAM PATHS ###
CELLRANGER.PATH = '/home/libraseq/.local/opt/cellranger-3.1.0/cellranger-cs/3.1.0/bin/cellranger'
CHANGEO.MAKEDB.PATH = '/home/libraseq/.local/bin/MakeDb.py'
CHANGEO.PARSEDB.PATH = '/home/libraseq/.local/bin/ParseDb.py'
CHANGEO.DEFINECLONES.PATH = '/home/libraseq/.local/bin/DefineClones.py'

#########################################################

# "Copyright 2022 Vanderbilt University. All rights reserved"

#########################################################
